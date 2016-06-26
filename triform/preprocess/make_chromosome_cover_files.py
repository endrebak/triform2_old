from collections import OrderedDict
from joblib import Parallel, delayed

from rpy2.robjects import r, pandas2ri
ri2py = pandas2ri.ri2py
from rpy2.robjects.packages import importr


def make_chromosome_cover_files(chromosome_granges, args):

    ranged_data = Parallel(n_jobs=args.number_cores)(
        delayed(_make_chromosome_cover_files)(df, args)
        for df in chromosome_granges.values())

    assert len(ranged_data) == len(chromosome_granges)
    ranged_data = OrderedDict([(
        c, d) for (c, d) in zip(chromosome_granges.keys(), ranged_data)])
    return ranged_data


def _make_chromosome_cover_files(granges, args):
    """Write description as command ending in a period."""

    importr("GenomicRanges")

    make_chromosome_cover_files_func = r("""
makeChromosomeCoverFiles <- function(rd, gapped.width){
  options(stringsAsFactors=FALSE)

  strands <- factor(as.vector(strand(rd)), levels=c("+", "-"))
  ly <- split(rd, strands)
  lsize <- lapply(ly, length)
  lstart <- lapply(ly, start)
  lwidth <- lapply(ly, width)
  lgap <- lapply(lwidth, function(w) floor((gapped.width-w)/2))
  irpos <- IRanges(start=lstart[["+"]] - lgap[["+"]], width=lwidth[["+"]] + 2*lgap[["+"]])
  irneg <- IRanges(start=lstart[["-"]] - lgap[["-"]], width=lwidth[["-"]] + 2*lgap[["-"]])
  irl <- IRangesList("-"=irneg, "+"=irpos)
  lcvg <- coverage(irl)
  covers <- list(SIZE=lsize, CVG=lcvg)
  # save(covers, file="covers.RData")
  covers
}
""")

    return make_chromosome_cover_files_func(granges, args.read_width)
