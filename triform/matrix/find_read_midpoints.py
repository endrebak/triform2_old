from collections import defaultdict

from rpy2.robjects import r
from joblib import Parallel, delayed


def find_read_midpoints(iranges, args):

    read_midpoints = defaultdict(list)
    for filename, irange_dict in iranges.items():

        genomic_ranges = Parallel(n_jobs=args.number_cores)(
            delayed(_find_read_midpoints)(df, chromosome, filename,
                                          args.matrix_bin_size)
            for chromosome, df in irange_dict.items())

        for chromosome, gr in irange_dict.items():
            read_midpoints[chromosome].append((filename, gr))

        # read_midpoints[f] = {chromosome: gr
        #                      for (chromosome, gr) in zip(irange_dict.keys(),
        #                                                  genomic_ranges)}

    return read_midpoints


def _find_read_midpoints(data, chromosome, filename, tilewidth=10):

    iranges_to_midpoint = r("""function(ir) {
    midpoint = (start(ir) + end(ir)) / 2
    coverage(IRanges(midpoint, midpoint))
    }""")

    midpoint = iranges_to_midpoint(data)

    rle_to_granges = r(
        """function(runlengths, chromosome, filename, tilewidth) {
    gr = GRanges(chromosome,IRanges(cumsum(c(0,runLength(runlengths)[-nrun(runlengths)])),
                              width=runLength(runlengths)),
             runlengths = runValue(runlengths))
    gr = gr[elementMetadata(gr)[,1] != 0]
    seqlengths(gr) = suppressWarnings(max(end(gr)))
    tg = tileGenome(seqinfo(gr), tilewidth=tilewidth, cut.last.tile.in.chrom=TRUE)

    runlengthsL = RleList(runlengths)
    names(runlengthsL) = chromosome

    ba = binnedAverage(trim(tg), runlengthsL, varname="runlengths")
    ba = ba[elementMetadata(ba)[,1] != 0]
    elementMetadata(ba)[,1] = elementMetadata(ba)[,1] * tilewidth
    names(values(ba)) = filename
    ba
    }
    """)

    return rle_to_granges(midpoint, chromosome, filename, tilewidth)
