from collections import OrderedDict
from joblib import Parallel, delayed

from rpy2.robjects import r, pandas2ri
ri2py = pandas2ri.ri2py
from rpy2.robjects.packages import importr


def make_ranged_data(chromosome_dfs, args):

    ranged_data = Parallel(n_jobs=args.number_cores)(
        delayed(_make_ranged_data)(df) for df in chromosome_dfs.values())

    assert len(ranged_data) == len(chromosome_dfs)
    ranged_data = OrderedDict([(
        c, d) for (c, d) in zip(chromosome_dfs.keys(), ranged_data)])
    return ranged_data


def _make_ranged_data(lines):
    """Write description as command ending in a period."""

    importr("GenomicRanges")

    # TODO: why does this function take a string and not a file?
    make_ranged_data_func = r("""
makeRangedData <- function(lines){
  options(stringsAsFactors=FALSE)
  con <- textConnection(lines)
  dfr <- read.delim(con, colClasses=c("character", rep("integer",2), "character"), header=FALSE)
  close(con)
  # dfr <- read.delim(infile, header=FALSE,
  #                   colClasses=c("character", rep("integer",2), "character"))
  colnames(dfr) <- c("seqnames", "start", "end", "strand")
  rd = makeGRangesFromDataFrame(dfr)
  # rd <- as(dfr, "RangedData")
  rd
  # save(rd, file="ranged.RData")
  # rd
}
""")

    return make_ranged_data_func(lines)
