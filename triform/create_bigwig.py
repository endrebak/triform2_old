import pandas as pd

from rpy2.robjects import r, pandas2ri, globalenv
from rpy2.robjects.packages import importr
ri2py = pandas2ri.ri2py
pandas2ri.activate()

from triform.helper_functions import df_to_iranges, subset_RS4


def create_bigwig(chip, fdr_table):

    locs = fdr_table.index.get_level_values(0).to_series()
    locs = locs.str.split(":", expand=True).iloc[:, 1]
    locs = locs.str.split("-", expand=True).astype(int)
    locs.insert(2, 2, locs.loc[:, 1] - locs.iloc[:, 0] + 1)

    dfs = []
    for chromosome, d in chip.items():
        chromosome_key = chromosome + ":"
        idx = locs.index.get_level_values(0).str.contains(chromosome_key)
        df = _create_bigwig(d, locs[idx], chromosome)
        dfs.append(df)

    df = pd.concat(dfs, axis=0)
    df.to_csv("chrY.bedgraph", sep=" ", index=False, header=False)


def _create_bigwig(chip, fdr_table, chromosome):

    center = [v for k, v in chip.items() if "center" == k[2]]
    center_sum = r["Reduce"]("+", center)

    globalenv["df"] = fdr_table
    df = r["df"]
    iranges = df_to_iranges(df)

    overlaps = r["Views"](center_sum, iranges)

    views = r["viewApply"](overlaps, "as.vector")
    result_dfs = []
    for (_, (start, end, _)), v in zip(fdr_table.iterrows(), views):
        v = pd.Series(ri2py(v))
        starts = pd.Series(range(start, end + 1))
        ends = starts + 1
        chromosomes = pd.Series([chromosome] * len(v))
        result_df = pd.concat([chromosomes, starts, ends, v], axis=1)
        result_dfs.append(result_df)

    return pd.concat(result_dfs, axis=0)
