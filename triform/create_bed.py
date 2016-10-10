import logging
from collections import defaultdict
from os.path import join, basename, splitext
from joblib import Parallel, delayed

import pandas as pd
from numpy import log2

from rpy2.robjects import r, pandas2ri
ri2py = pandas2ri.ri2py
from rpy2.robjects.packages import importr
from triform.create_bedgraph import create_bedgraph
from triform.helper_functions import df_to_rle, rle_to_df

importr("GenomicRanges")
importr("rtracklayer")
pandas2ri.activate()


def _enriched_regions_to_dict(enriched_regions):

    enriched_regions.to_csv("enriched_regions.csv", sep=" ")
    d = {}
    df_to_ir = r("function(df) {IRanges(start=df$START, end=df$END)}")
    for chromosome, df in enriched_regions.groupby("CHROMOSOME"):
        df = df["CHROMOSOME START END QVAL".split()]
        ir = df_to_ir(df)
        gr = r["GRanges"](chromosome, ir)

        d[chromosome] = gr

    return d



def overlaps_enriched_regions_and_data(iranges, enriched_dict):

    find_overlaps = r("function(query, subject) findOverlaps(query=query, subject=subject)")
    get_queryhits = r("function(overlaps, query) query[queryHits(overlaps),]")
    get_subjecthits = r("function(overlaps, subject) subject[subjectHits(overlaps),]")

    to_dataframe = r(""" function(subject, query) {
    cdf = cbind(as.data.frame(subject), as.data.frame(query))
    colnames(cdf) = 1:dim(cdf)[2]
    cdf
    }""")

    chip_granges = defaultdict(list)

    for filename, irs in iranges.items():
        for chromosome, enriched in enriched_dict.items():
            data = irs[chromosome]
            gr = r["GRanges"](chromosome, data)

            hits_q_enriched = find_overlaps(enriched, gr)
            subjecthits = get_subjecthits(hits_q_enriched, gr)
            queryhits = get_queryhits(hits_q_enriched, enriched)

            df = to_dataframe(subjecthits, queryhits)

            idxes = "6 7 8".split()

            df = ri2py(df).ix[:, idxes]

            df.columns = "Chromosome Start End".split()

            sizes = df.groupby("Chromosome Start End".split()).size()
            chip_granges[chromosome].append(sizes)


    chip_counts = {}
    for chromosome, granges_list in chip_granges.items():
        counts = pd.concat(granges_list, axis=1).sum(axis=1)

        chip_counts[chromosome] = counts

    return chip_counts


def create_bed(chip_iranges, input_iranges, enriched_regions, args):

    enriched_regions = enriched_regions.reset_index()
    enriched_dict = _enriched_regions_to_dict(enriched_regions)

    chip_counts = overlaps_enriched_regions_and_data(chip_iranges, enriched_dict)
    input_counts = overlaps_enriched_regions_and_data(input_iranges, enriched_dict)

    chromosome_list = []
    for chip, input in zip(chip_counts.values(), input_counts.values()):
        chip.name = "ChIP"
        input.name = "Input"
        df = chip.to_frame().join(input, how="outer")
        chromosome_list.append(df)

    df = pd.concat(chromosome_list)
    enriched_regions.columns = [c.capitalize() for c in enriched_regions.columns]
    enriched_regions = enriched_regions.set_index("Chromosome Start End".split())
    df = df.join(enriched_regions[["Qval", "Loc"]])

    df = df.fillna(0)
    df.loc[df.Input == 0, "Input"] = 0.5
    df.insert(4, "logFC", log2(df.ChIP / (df.Input)) * 100)
    df.loc[df.logFC > 1000, "logFC"] = 1000

    return df["Qval logFC Loc".split()]



def _create_bed(chip_iranges, input_iranges, outfile, args):

    logging.info("Creating bed file " + outfile)
    genomicranges = []

    for chromosome, data in iranges.items():
        gr = r["GRanges"](chromosome, data)
        genomicranges.append(gr)

    genome_coverage = r["coverage"](r["c"](*genomicranges))
    genome_coverage_rpkm = r("function(cv) {1e6 * cv / sum(sum(cv))}")(
        genome_coverage)
    genome_coverage_rpkm = r("function(n, o) {names(n) = names(o); n}")(
        genome_coverage_rpkm, genome_coverage)

    r["export.bed"](genome_coverage_rpkm, outpath)
