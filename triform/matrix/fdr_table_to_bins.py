from collections import defaultdict
import pandas as pd


def fdr_table_to_bins(fdr_table, tilewidth):

    fdr_table = fdr_table.reset_index()
    rowdicts = []
    for chromosome, start, end in zip(fdr_table.CHROMOSOME, fdr_table.START,
                                      fdr_table.END):
        region_start = start - (start % tilewidth) + 1
        region_end = end - (end % tilewidth) + tilewidth
        for start in range(region_start, region_end + 1, tilewidth):
            rowdicts.append({"Chromosome": chromosome,
                             "Start": start,
                             "End": start + tilewidth - 1,
                             "Enriched": 1})

    df = pd.DataFrame.from_dict(rowdicts)
    df = df["Chromosome Start End Enriched".split()].set_index(
        "Chromosome Start End".split())

    d = defaultdict(lambda: pd.DataFrame())

    for chromosome in df.index.get_level_values(0).drop_duplicates():
        d[chromosome] = df.xs(chromosome, drop_level=False)

    return d
