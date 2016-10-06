import logging
from os.path import join, basename, splitext
from joblib import Parallel, delayed

from rpy2.robjects import r, pandas2ri
ri2py = pandas2ri.ri2py
from rpy2.robjects.packages import importr
from triform.create_bedgraph import create_bedgraph
from triform.helper_functions import df_to_rle, rle_to_df

importr("GenomicRanges")
importr("rtracklayer")


def create_bigwig(data, outfolder, args):

    matrix_outfiles = Parallel(n_jobs=args.number_cores)(
        delayed(_create_bigwig)(iranges, filename, args.bigwig, args)
        for filename, iranges in data.items())


def _create_bigwig(iranges, filename, outfolder, args):

    outpath = join(outfolder, splitext(basename(filename))[0] + ".bw")

    logging.info("Creating biwgwig " + outpath)
    genomicranges = []

    for chromosome, data in iranges.items():
        gr = r["GRanges"](chromosome, data)
        genomicranges.append(gr)

    genome_coverage = r["coverage"](r["c"](*genomicranges))
    genome_coverage_rpkm = r("function(cv) {1e6 * cv / sum(sum(cv))}")(
        genome_coverage)
    genome_coverage_rpkm = r("function(n, o) {names(n) = names(o); n}")(
        genome_coverage_rpkm, genome_coverage)

    r["export.bw"](genome_coverage_rpkm, outpath)
