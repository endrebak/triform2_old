import pandas as pd

from natsort import natsorted

from rpy2.robjects import r, pandas2ri, globalenv
from rpy2.robjects.packages import importr
ri2py = pandas2ri.ri2py
pandas2ri.activate()

from triform.helper_functions import df_to_iranges, subset_RS4


def create_bigwig(chip, fdr_table, args):

    locs = fdr_table.index.get_level_values(0).to_series()
    locs = locs.str.split(":", expand=True).iloc[:, 1]
    locs = locs.str.split("-", expand=True).astype(int)
    locs.insert(2, 2, locs.loc[:, 1] - locs.iloc[:, 0] + 1)

    dfs = []
    # TODO: use multiple cores
    for chromosome, d in natsorted(chip.items()):
        chromosome_key = chromosome + ":"
        idx = locs.index.get_level_values(0).str.contains(chromosome_key)
        chromosome_df = locs[idx]

        # if chromosome_df.empty:
        #     continue

        df = _create_bigwig(d, locs[idx], chromosome)
        dfs.append(df)

    df = pd.concat(dfs, axis=0)
    df.to_csv(args.bedgraph, sep=" ", index=False, header=False)
    return df


def _create_bigwig(chip, fdr_table, chromosome):

    center = [v for k, v in chip.items() if "center" == k[2]]
    center_sum = r["Reduce"]("+", center)

    globalenv["df"] = fdr_table
    df = r["df"]
    iranges = df_to_iranges(df)

    overlaps = r["Views"](center_sum, iranges)
    length = r["length"](overlaps)[0]

    if length > 1:
        views = r["viewApply"](overlaps, "as.vector")
        views = [views.rx2(i) for i in range(1, len(views) + 1)]
    else:
        views = r["list"](r["viewApply"](overlaps, "as.vector"))

    rowdicts = []

    for (_, (start, end, _)), v in zip(fdr_table.iterrows(), views):
        v = list(ri2py(v))

        current_count, current_start = v[0], 0
        for pos, count in enumerate(v):

            if count == current_count:
                continue

            rowdict = {"Chromosome": chromosome,
                       "Start": start + current_start,
                       "End": start + pos + 1,
                       "Count": current_count}
            rowdicts.append(rowdict)

            current_start = pos + 1
            current_count = count
            current_end = start + pos + 1

        if count != current_count:
            continue

        if end != current_end:
            rowdict = {"Chromosome": chromosome,
                       "Start": start + current_start,
                       "End": end,
                       "Count": current_count}

            print(rowdict)
            if rowdict["Start"] < rowdict["End"]:
                rowdicts.append(rowdict)

    if rowdicts:
        df = pd.DataFrame.from_dict(rowdicts)[["Chromosome", "Start", "End",
                                               "Count"]]
        df["Count"] = df["Count"].astype(int)
        return df.sort_values("Start").reset_index(drop=True)
    else:
        return pd.DataFrame()
