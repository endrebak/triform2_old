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

        for chromosome, gr in zip(irange_dict.keys(), genomic_ranges):
            read_midpoints[chromosome].append((filename, gr))

    return read_midpoints


def _find_read_midpoints(data, chromosome, filename, tilewidth=10):

    iranges_to_midpoint = r("""function(ir) {
    midpoint = (start(ir) + end(ir)) / 2
    cv = coverage(IRanges(midpoint, midpoint))
    cv
    }""")

    midpoint = iranges_to_midpoint(data)

    rle_to_tiles = r("""function(runlengths, chromosome, filename, tilewidth) {

        gr = GRanges(chromosome,IRanges(start(runlengths) - 1,
       width=runLength(runlengths)),
        runlengths = runValue(runlengths))

    gr2 = GRanges(chromosome, IRanges(max(start(runlengths) + tilewidth * 2),
          width=tilewidth), runlengths=1)
    # adding gr2 since cut.last.tile.in.chrom must be true
    gr = c(gr, gr2)
    gr = gr[elementMetadata(gr)[,1] != 0]
    seqlengths(gr) = suppressWarnings(max(end(gr)))
    tg = tileGenome(seqinfo(gr), tilewidth=tilewidth, cut.last.tile.in.chrom=TRUE)
    tg
    }
        """)

    tiles = rle_to_tiles(midpoint, chromosome, filename, tilewidth)

    give_tiles_counts = r(
        """function(runlengths, tg, chromosome, filename, tilewidth){{
    runlengthsL = RleList(runlengths)
    names(runlengthsL) = chromosome

    ba = binnedAverage(tg, runlengthsL, varname="runlengths")
    ba = ba[elementMetadata(ba)[,1] != 0]
    elementMetadata(ba)[,1] = elementMetadata(ba)[,1] * tilewidth
    names(values(ba)) = "counts"
    ba
    }}
    """)

    ct = give_tiles_counts(midpoint, tiles, chromosome, filename, tilewidth)

    return ct
