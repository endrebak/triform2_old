from __future__ import division

import logging

import numpy as np

from itertools import product
from collections import defaultdict

from joblib import Parallel, delayed

from rpy2.robjects import r, pandas2ri
ri2py = pandas2ri.ri2py
from rpy2.robjects.packages import importr
importr("GenomicRanges")
importr("S4Vectors")

from triform.helper_functions import subset_RS4
import triform.config.logging_settings


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


def compute_ratios(chip_sizes, input_size):
    ratios = dict()
    for f, d in chip_sizes.items():
        ratios[f] = input_size / d

    return ratios


def zscores(x, y, ratio=1):

    _zscores = r("""
function(x,y,r) {  # r = size.y/size.x
  dif <- (r*x-y)
  zs <- dif/sqrt(r*(x+y))
  zs[!dif] <- 0
  zs
}
""")
    return _zscores(x, y, ratio)


def zscores_no_control(x, ratio=1):

    _zscores = r("""
function(x,r) {  # r = size.y/size.x
  dif <- (r*x)
  zs <- dif/sqrt(r*(x))
  zs[!dif] <- 0
  zs
}
""")
    return _zscores(x, ratio)


def compute_peaks_and_zscores(cvg, center, left, right, chip, input, ratios,
                              ratio, args):
    # TODO: rewrite as a loop
    min_z = r["qnorm"](args.max_p, lower_tail=False)
    left_right = r["+"](left, right)

    ok1 = compute_ok1(chip)
    ok2 = compute_ok23(chip, 2)
    ok3 = compute_ok23(chip, 3)

    if input:
        ok4 = compute_ok4(ratios, center, input)

    zscores1 = r["*"](ok1, zscores(cvg, left_right, 2))
    zscores2 = r["*"](ok2, zscores(cvg, left))
    zscores3 = r["*"](ok3, zscores(cvg, right))
    if input:
        zscores4 = r["*"](ok4, zscores(cvg, input, ratio))
    else:
        zscores4 = zscores_no_control(cvg, 1)

    peaks1 = r["slice"](zscores1, lower=min_z)
    peaks2 = r["slice"](zscores2, lower=min_z)
    peaks3 = r["slice"](zscores3, lower=min_z)
    peaks4 = r["slice"](zscores4, lower=min_z)
    # print(peaks1)
    # print(type(peaks1))
    # print(r["typeof"](peaks1))

    subset1 = r[">"](r["width"](peaks1), args.min_width)
    subset2 = r[">"](r["width"](peaks2), args.min_width)
    subset3 = r[">"](r["width"](peaks3), args.min_width)
    subset4 = r[">"](r["width"](peaks4), args.min_width)

    peaks1 = subset_RS4(peaks1, subset1)
    peaks2 = subset_RS4(peaks2, subset2)
    peaks3 = subset_RS4(peaks3, subset3)
    peaks4 = subset_RS4(peaks4, subset4)
    # print(subset1)
    # print(type(subset1))
    # print(r["typeof"](subset1))
    # print(peaks1)
    # print(type(peaks1))
    # print(r["typeof"](peaks1))

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

    return _peaks, _zscores


def transform_input_data(data):

    chromosome_dict = defaultdict(
        lambda: defaultdict(lambda: defaultdict(dict)))
    for chromosome, d in data.items():
        for (infile, direction, location), v in d.items():
            chromosome_dict[chromosome][direction][infile][location] = v

    return chromosome_dict


def transform_input_sizes(sizes):

    chromosome_dict = defaultdict(lambda: defaultdict(dict))
    for chromosome, d in sizes.items():
        for (infile, direction), v in d.items():
            chromosome_dict[chromosome][direction][infile] = v

    return chromosome_dict


def chromosome(chip_data, input_data, chip_sizes, input_sizes, args):

    chip_data = transform_input_data(chip_data)
    chip_sizes = transform_input_sizes(chip_sizes)

    if not input_data:
        keys = list(product(chip_data, ["forward", "reverse"]))

        results = Parallel(n_jobs=args.number_cores)(
            delayed(_chromosome)(chip_data[chromosome][direction], None,
                                 chip_sizes[chromosome][direction], None, args)
            for chromosome, direction in keys)

        results_dict = {k: v for k, v in zip(keys, results)}

        return results_dict

    input_data = transform_input_data(input_data)
    input_sizes = transform_input_sizes(input_sizes)

    assert len(chip_data) == len(input_data) == len(chip_sizes) == len(
        input_sizes), "chip_data, input_data, chip_sizes, and input_sizes do not all have the same length!"

    keys = list(product(chip_data, ["forward", "reverse"]))

    results = Parallel(n_jobs=args.number_cores)(delayed(_chromosome)(
        chip_data[chromosome][direction], input_data[chromosome][direction],
        chip_sizes[chromosome][direction], input_sizes[chromosome][direction],
        args) for chromosome, direction in keys)

    results_dict = {k: v for k, v in zip(keys, results)}

    return results_dict


def _chromosome(chip, input_data, chip_sizes, input_data_sizes, args):

    results = defaultdict(dict)

    if input_data:
        input_data = list(collect_key(input_data, "center").values())
        input_data = r["Reduce"]("+", input_data)

        ratios = compute_ratios(chip_sizes, sum(input_data_sizes.values()))

        ratio = sum(input_data_sizes.values()) / sum(chip_sizes.values())
    else:
        ratio = 1
        ratios = None

    center = collect_key(chip, "center")
    left = collect_key(chip, "left")
    right = collect_key(chip, "right")

    cvg = r["Reduce"]("+", list(center.values()))
    left = r["Reduce"]("+", list(left.values()))
    right = r["Reduce"]("+", list(right.values()))

    results["cvg"] = cvg
    _peaks, _zscores = compute_peaks_and_zscores(
        cvg, center, left, right, chip, input_data, ratios, ratio, args)

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
        if input_data:
            peak_ref = subset_RS4(input_data, peak_locs, True)

        if input_data:

            peak_enrich_cvg = r["+"](1, r["*"](ratio, peak_cvg))
            peak_enrich_ref = r["+"](1, peak_ref)
            peak_enrich = r["/"](peak_enrich_cvg, peak_enrich_ref)
        else:
            peak_enrich = peak_cvg

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

        # perhaps we get away with not returning nb_peaks. Let's see!
        # results["nb_peaks"][peak_type] = r["length"](peaks)
        results["peaks"][peak_type] = peaks
        results["peak_info"][peak_type] = dfr

    return results
