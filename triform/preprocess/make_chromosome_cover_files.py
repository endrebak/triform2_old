from collections import OrderedDict, defaultdict
from joblib import Parallel, delayed

from rpy2.robjects import r, pandas2ri
ri2py = pandas2ri.ri2py
from rpy2.robjects.packages import importr

from triform.helper_functions import subset_RS4_2


def make_chromosome_cover_files(ranged_data_per_file, args):

    cvgs = defaultdict(dict)
    sizes = defaultdict(dict)
    for f, chromosome_granges in ranged_data_per_file.items():

        ranged_data = Parallel(n_jobs=args.number_cores)(
            delayed(_make_chromosome_cover_files)(df, args)
            for df in chromosome_granges.values())

        assert len(ranged_data) == len(chromosome_granges)

        for c, (cvg, size) in zip(chromosome_granges.keys(), ranged_data):

            cvgs[c][f, "forward"] = cvg["+"]
            cvgs[c][f, "reverse"] = cvg["-"]

            sizes[c][f, "forward"] = size["+"]
            sizes[c][f, "reverse"] = size["-"]

    return cvgs, sizes


def _make_chromosome_cover_files(granges, args):
    """Write description as command ending in a period."""

    importr("GenomicRanges")

    strands = r(
        'function(rd) factor(as.vector(strand(rd)), levels=c("+", "-"))')(
            granges)
    strand_lists = r["split"](granges, strands)
    strand_sizes = r["lapply"](strand_lists, "length")
    strand_starts = r["lapply"](strand_lists, "start")
    strand_widths = r["lapply"](strand_lists, "width")
    strand_gaps = r["lapply"](
        strand_widths,
        r("function(w) floor(({}-w)/2)".format(args.read_width)))

    get_iranges = r('''function(lstart, lgap, lwidth, strand) {
    IRanges(start=lstart[[strand]] - lgap[[strand]], width=lwidth[[strand]] + 2*lgap[[strand]])
    }''')
    forward_iranges = get_iranges(strand_starts, strand_gaps, strand_widths,
                                  "+")
    reverse_iranges = get_iranges(strand_starts, strand_gaps, strand_widths,
                                  "-")

    iranges_list = r(
        'function(irneg, irpos) IRangesList("-"=irneg, "+"=irpos)')(
            reverse_iranges, forward_iranges)
    cvg = r["coverage"](iranges_list)

    cvgs = {"+": subset_RS4_2(cvg, "+"), "-": subset_RS4_2(cvg, "-")}
    sizes = {"+": subset_RS4_2(strand_sizes, "+"),
             "-": subset_RS4_2(strand_sizes, "-")}

    return cvgs, sizes
