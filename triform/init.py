from itertools import product

from rpy2 import robjects as ro
from rpy2.robjects import r, pandas2ri
ri2py = pandas2ri.ri2py
from rpy2.robjects.packages import importr
importr("GenomicRanges")


def init(covers, is_input, args):

    replicates = set(k[0] for k in covers.keys())

    results = dict()

    flank_delta_pad = r["Rle"](0, args.flank_distance)

    for replicate, direction, location in product(
            replicates, ["reverse", "forward"], ["left", "right", "center"]):

        cvg = covers[replicate, direction]
        if is_input:
            if location == "center":
                results[replicate, direction, location] = cvg
            continue

        if location == "left":
            cvg = r('''function(cvg, flank.delta.pad, flank.delta) {
              c(flank.delta.pad, rev(rev(cvg)[-1:-flank.delta]))
            }''')(cvg, flank_delta_pad, args.flank_distance)
        if location == "right":
            cvg = r('''function(cvg, flank.delta.pad, flank.delta) {
              c(cvg[-1:-flank.delta], flank.delta.pad)
            }''')(cvg, flank_delta_pad, args.flank_distance)

        results[replicate, direction, location] = cvg

    maxlen = r["max"](r["sapply"](results.values(), "length"))
    lapply = r('function(cvg, maxlen) c(cvg,Rle(0,maxlen-length(cvg)))')
    results = {k: lapply(v, maxlen) for (k, v) in results.items()}

    return results
