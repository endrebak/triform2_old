from __future__ import division
import pandas as pd
import numpy as np

from collections import defaultdict

import pytest

from rpy2.robjects import r, pandas2ri
from rpy2.robjects.packages import importr
importr("GenomicRanges")
importr("S4Vectors")

from triform.init import _init


def df_to_rle(df):

    _df_to_rle = r("""function(df){
    Rle(df$values, df$lengths)
    }
    """)

    return _df_to_rle(df)


def _create_rle_list(run_length_files):

    # run_lengths = r["list"]()
    run_lengths = dict()
    for name, rle in run_length_files.items():
        df = r["read.table"](rle, sep=" ", header=1)
        run_lengths[name] = df_to_rle(df)

    return run_lengths


@pytest.fixture
def input_data(run_length_encodings_full_rep1_backgr,
               run_length_encodings_full_rep2_backgr):

    input = defaultdict(list)

    rep1 = _create_rle_list(run_length_encodings_full_rep1_backgr)
    rep2 = _create_rle_list(run_length_encodings_full_rep2_backgr)

    input["forward"] = [rep1["forward"], rep2["forward"]]
    input["reverse"] = [rep1["reverse"], rep2["reverse"]]

    return input


@pytest.fixture
def chip_data(run_length_encodings_full_rep1, run_length_encodings_full_rep2):

    data = defaultdict(dict)
    data["forward"]["rep1"] = _create_rle_list(run_length_encodings_full_rep1[
        "forward"])
    data["forward"]["rep2"] = _create_rle_list(run_length_encodings_full_rep2[
        "forward"])

    data["reverse"]["rep1"] = _create_rle_list(run_length_encodings_full_rep1[
        "reverse"])
    data["reverse"]["rep2"] = _create_rle_list(run_length_encodings_full_rep2[
        "reverse"])

    return data


@pytest.fixture
def chip_sizes(sizes_rep1, sizes_rep2):
    d = defaultdict(dict)
    d["reverse"]["rep1"] = sizes_rep1["reverse"]
    d["forward"]["rep1"] = sizes_rep1["forward"]
    d["reverse"]["rep2"] = sizes_rep2["reverse"]
    d["forward"]["rep2"] = sizes_rep2["forward"]
    return d


@pytest.fixture
def input_sizes(sizes_rep1_backgr, sizes_rep2_backgr):
    d = defaultdict(dict)
    d["reverse"]["rep1"] = sizes_rep1_backgr["reverse"]
    d["forward"]["rep1"] = sizes_rep1_backgr["forward"]
    d["reverse"]["rep2"] = sizes_rep2_backgr["reverse"]
    d["forward"]["rep2"] = sizes_rep2_backgr["forward"]
    return d


@pytest.fixture
def expected_result():
    return None


@pytest.mark.current
def test_chromosome(chip_data, input_data, chip_sizes, input_sizes, args,
                    expected_result):

    # r["print"](r["sessionInfo"]())
    result = chromosome(chip_data, input_data, chip_sizes, input_sizes, args)
    assert 0  # result == expected_result

# MAX.P <- 0.1
# MIN.Z <- qnorm(MAX.P,low=FALSE)
# MIN.SHIFT <- 10
# MIN.WIDTH <- 10


def collect_key(d, key):
    d2 = {}
    for f, v in d.items():
        d2[f] = v[key]
    return d2


def compute_ok1(d):
    signs = []
    for v in d.values():
        double = r["*"](2, v["center"])
        s = r["sign"](r["Reduce"]("-", [double, v["left"], v["right"]]))
        pos = r["=="](1, s)
        signs.append(pos)
    return r["Reduce"]("*", signs)


def compute_ok23(d, type):
    type = "left" if type == 2 else "right"
    signs = []
    for v in d.values():
        s = r["sign"](r["-"](v["center"], v[type]))
        pos = r["=="](1, s)
        signs.append(pos)
    return r["Reduce"]("*", signs)


def compute_ok4(ratios, center, input):

    assert ratios.keys() == center.keys()

    signs = []
    for k, ratio in ratios.items():
        c = center[k]
        scaled_center = r["*"](c, ratio)
        s = r["sign"](r["-"](scaled_center, input))
        pos = r["=="](1, s)
        signs.append(pos)
    return r["Reduce"]("*", signs)

    # stack_overflow()

    # type = "left" if type == 2 else "right"
    # signs = []
    # for v in d.values():
    #     s = r["sign"](r["-"](v["center"], v[type]))
    #     pos = r["=="](1, s)
    #     signs.append(pos)
    # return r["Reduce"]("*", signs)


def subset_RS4(rs4, subset, drop=False):
    subset_func = r("""function(o, s, d){
    o[s, drop=d]
    }
    """)
    return subset_func(rs4, subset, drop)


def compute_ratios(chip_sizes, input_size):
    ratios = dict()
    for f, d in chip_sizes.items():
        ratios[f] = input_size / d
    return ratios


def zscores(x, y, ratio=1):

    _zscores = r("""
function(x,y,r=1) {  # r = size.y/size.x
  dif <- (r*x-y)
  zs <- dif/sqrt(r*(x+y))
  zs[!dif] <- 0
  zs
}
""")
    return _zscores(x, y, ratio)


def chromosome(chip_covers, input_covers, chip_sizes, input_sizes, args):

    results = defaultdict(dict)

    min_z = r["qnorm"](args.max_p, lower_tail=False)

    reverse_input = input_covers["reverse"]
    reverse_input = r["Reduce"]("+", reverse_input)
    input_size = sum(input_sizes["reverse"].values())
    chip_sizes = chip_sizes["reverse"]

    ratios = compute_ratios(chip_sizes, input_size)

    ratio = input_size / sum(chip_sizes.values())

    reverse = chip_covers["reverse"]

    center = collect_key(reverse, "center")
    left = collect_key(reverse, "left")
    right = collect_key(reverse, "right")

    cvg = r["Reduce"]("+", list(center.values()))
    left = r["Reduce"]("+", list(left.values()))
    right = r["Reduce"]("+", list(right.values()))

    results["cvg"] = cvg
    left_right = r["+"](left, right)

    ok1 = compute_ok1(reverse)
    ok2 = compute_ok23(reverse, 2)
    ok3 = compute_ok23(reverse, 3)
    ok4 = compute_ok4(ratios, center, reverse_input)

    zscores1 = r["*"](ok1, zscores(cvg, left_right, 2))
    zscores2 = r["*"](ok2, zscores(cvg, left))
    zscores3 = r["*"](ok3, zscores(cvg, right))
    zscores4 = r["*"](ok4, zscores(cvg, reverse_input, ratio))

    peaks1 = r["slice"](zscores1, lower=min_z)
    peaks2 = r["slice"](zscores2, lower=min_z)
    peaks3 = r["slice"](zscores3, lower=min_z)
    peaks4 = r["slice"](zscores4, lower=min_z)

    subset1 = r[">"](r["width"](peaks1), args.min_width)
    subset2 = r[">"](r["width"](peaks2), args.min_width)
    subset3 = r[">"](r["width"](peaks3), args.min_width)
    subset4 = r[">"](r["width"](peaks4), args.min_width)

    peaks1 = subset_RS4(peaks1, subset1)
    peaks2 = subset_RS4(peaks2, subset2)
    peaks3 = subset_RS4(peaks3, subset3)
    peaks4 = subset_RS4(peaks4, subset4)

    peaks1 = r["as"](peaks1, "IRanges")
    peaks2 = r["as"](peaks2, "IRanges")
    peaks3 = r["as"](peaks3, "IRanges")
    peaks4 = r["as"](peaks4, "IRanges")

    peaks1 = r["intersect"](peaks1, peaks4)
    peaks2 = r["intersect"](peaks2, peaks4)
    peaks3 = r["intersect"](peaks3, peaks4)

    subset1 = r[">"](r["width"](peaks1), args.min_width)
    subset2 = r[">"](r["width"](peaks2), args.min_width)
    subset3 = r[">"](r["width"](peaks3), args.min_width)

    peaks1 = subset_RS4(peaks1, subset1)
    peaks2 = subset_RS4(peaks2, subset2)
    peaks3 = subset_RS4(peaks3, subset3)

    _peaks = [peaks1, peaks2, peaks3]
    _zscores = [zscores1, zscores2, zscores3]

    views_func = r("function(x,y) Views(x,y)")

    zviews_list = r["mapply"](views_func, x=_zscores, y=_peaks)
    maxz_list = r["lapply"](zviews_list, r["viewMaxs"])

    for peak_type, peaks in enumerate(_peaks, 1):
        maxz = maxz_list.rx2(peak_type)
        peak_nlps = r["-"](r["/"](r["pnorm"](maxz,
                                             lower_tail=False,
                                             log_p=True),
                                  np.log(10)))
        start_peaks = r["start"](peaks)
        end_peaks = r["end"](peaks)
        sum_starts_ends = r["+"](start_peaks, end_peaks)
        peak_locs = r["/"](sum_starts_ends, 2)
        peak_locs = r["round"](peak_locs)

        peak_cvg = subset_RS4(cvg, peak_locs, True)
        peak_ref = subset_RS4(reverse_input, peak_locs, True)

        peak_enrich_cvg = r["+"](1, r["*"](ratio, peak_cvg))
        peak_enrich_ref = r["+"](1, peak_ref)
        peak_enrich = r["/"](peak_enrich_cvg, peak_enrich_ref)

        # min_er will always be defined as it is set in the first loop, but ugh the code
        if peak_type == 1:
            min_er = r["quantile"](peak_enrich, args.min_enrichment)

        ok = r[">"](peak_enrich, min_er)
        nb_ok = r["sum"](ok)

        if not nb_ok:
            continue

        peaks = subset_RS4(peaks, ok)
        peak_locs = subset_RS4(peak_locs, ok)

        peak_nlps = subset_RS4(peak_nlps, ok)
        peak_cvg = subset_RS4(peak_cvg, ok)

        peak_left = subset_RS4(left, peak_locs, True)
        peak_right = subset_RS4(right, peak_locs, True)

        dfr = r["data.frame"](PEAK_LOC=peak_locs,
                              PEAK_CVG=peak_cvg,
                              PEAK_SURL=peak_left,
                              PEAK_SURR=peak_right,
                              PEAK_NLP=r["round"](peak_nlps, 3),
                              PEAK_WIDTH=r["width"](peaks),
                              PEAK_START=r["start"](peaks),
                              PEAK_END=r["end"](peaks))

        # reorder columns since rpy2 does not respect the col order above
        rename_cols = r(
            'function (df) { names(df) = sub("PEAK_", "PEAK.", names(df)); df}')
        dfr = rename_cols(dfr)
        column_order = r(
            'c("PEAK.LOC", "PEAK.CVG", "PEAK.SURL", "PEAK.SURR", "PEAK.NLP", "PEAK.WIDTH", "PEAK.START", "PEAK.END")')
        co = r["match"](column_order, r["names"](dfr))
        dfr = dfr.rx(True, co)

        results["nb_peaks"][peak_type] = r["length"](peaks)
        results["peaks"][peak_type] = peaks
        results["peak_info"][peak_type] = dfr

    r["write.table"](results["peak_info"][1], "dfr_py", sep=" ")
    return results
