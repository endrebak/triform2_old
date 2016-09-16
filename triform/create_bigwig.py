import pandas as pd
import numpy as np

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
    for chromosome, d in sorted(chip.items()):
        chromosome_key = chromosome + ":"
        idx = locs.index.get_level_values(0).str.contains(chromosome_key)
        chromosome_df = locs[idx]

        # if chromosome_df.empty:
        #     continue

        df = _create_bigwig(d, locs[idx], chromosome)
        dfs.append(df)

    df = pd.concat(dfs, axis=0)
    df.to_csv(args.bedgraph, sep="\t", index=False, header=False)
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

    dfs = []
    for (_, (start, _, _)), v in zip(fdr_table.iterrows(), views):
        counts = pd.Series(ri2py(v)).reset_index(drop=True)

        counts_shifted = counts.shift(-1)

        counts_changing = counts != counts_shifted

        idx = counts[counts_changing].index.get_level_values(0)

        values = counts[idx].reset_index(drop=True)

        sidx = pd.Series(idx)
        lengths = (sidx.shift(-1) - sidx).shift()
        lengths[0] = sidx[0]

        starts = sidx + start - sidx[0]

        ends = starts.shift(-1)
        ends.iloc[-1] = start + len(counts) - 1

        df = pd.concat([starts, ends, values], axis=1).astype(int)
        df.insert(0, "Chromosome", chromosome)
        df.columns = "Chromosome Start End Count".split()

        dfs.append(df)

    if dfs:
        return pd.concat(dfs,
                         axis=0).sort_values("Start").reset_index(drop=True)
    else:
        return pd.DataFrame()
