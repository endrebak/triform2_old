from joblib import Parallel, delayed
from rpy2.robjects import r, pandas2ri
ri2py = pandas2ri.ri2py
from rpy2.robjects.packages import importr
importr("GenomicRanges")
importr("S4Vectors")

from triform.helper_functions import (
    subset_RS4, subset_RS4_rows, subset_RS4_cols, df_to_iranges, df_to_rle)


def exclude_redundant_peaks(indata, args):

    return Parallel(n_jobs=args.number_cores)(
        delayed(_exclude_redundant_peaks)(data)
        for chromosome, data in indata.items())


def _exclude_redundant_peaks(indata):

    p1 = indata["peaks"][1]
    p2 = indata["peaks"][2]
    p3 = indata["peaks"][3]
    print("p1")
    print(p1)
    print("p2")
    print(p2)
    print("p3")
    print(p3)

    # # IRanges of length 1
    # #     start      end width
    # # [1] 57420858 57420942    85

    # # p2
    # # IRanges of length 1
    # #     start      end width
    # # [1] 10636518 10636576    59

    # # p3
    # # IRanges of length 2
    # #     start      end width
    # # [1] 10631247 10631396   150
    # # [2] 57420858 57420942    85

    print("indata['peak_info'][1]")
    print(indata['peak_info'][1])
    print("indata['peak_info'][2]")
    print(indata['peak_info'][2])
    print("indata['peak_info'][3]")
    print(indata['peak_info'][3])

    # indata['peak_info'][1]
    #                               NLP  MAX.NLP      LOC WIDTH    START      END CVG
    # chrY:57420858-57420942:1 3.923426 6.317243 57420916    85 57420858 57420942  10
    #                          SURL SURR FORM
    # chrY:57420858-57420942:1    0    2    1

    # indata['peak_info'][2]
    #                               NLP  MAX.NLP      LOC WIDTH    START      END CVG
    # chrY:10636518-10636576:2 2.869699 2.869699 10636545    59 10636518 10636576   9
    #                          SURL SURR FORM
    # chrY:10636518-10636576:2    0   16    2

    # indata['peak_info'][3]
    #                               NLP  MAX.NLP      LOC WIDTH    START      END CVG
    # chrY:10631247-10631396:3 4.039029 4.039029 10631334   150 10631247 10631396  14
    # chrY:57420858-57420942:3 1.980441 3.575114 57420916    85 57420858 57420942  10
    #                          SURL SURR FORM
    # chrY:10631247-10631396:3   24    0    3
    # chrY:57420858-57420942:3    0    2    3

    ov12 = r("function(p1, p2) matrix(as.matrix(findOverlaps(p1,p2)),ncol=2)")(
        p1, p2)
    if r["nrow"](ov12):
        ex2 = r("function(p2, ov12) 1:length(p2) %in% ov12[,2]")(p2, ov12)
        inv_ex2 = r["!"](ex2)
        p2 = subset_RS4_rows(p2, inv_ex2)
        info2 = subset_RS4_rows(indata["peak_info"][2], inv_ex2)

    ov13 = r("function(p1, p3) matrix(as.matrix(findOverlaps(p1,p3)),ncol=2)")(
        p1, p3)
    if r["nrow"](ov13):
        ex3 = r("function(p3, ov13) 1:length(p3) %in% ov13[,2]")(p3, ov13)
        inv_ex3 = r["!"](ex3)
        p3 = subset_RS4_rows(p3, inv_ex3)
        info3 = subset_RS4_rows(indata["peak_info"][3], inv_ex3)

    # TODO: do not use harcoded variable
    ov23 = r(
        "function(p2, p3) matrix(as.matrix(findOverlaps(p2,p3, maxgap=100)),ncol=2)")(
            p1, p3)
    if r["nrow"](ov23):
        ex2 = r("function(p2, ov23) 1:length(p2) %in% ov23[,1]")(p2, ov23)
        inv_ex2 = r["!"](ex2)
        p2 = subset_RS4_rows(p2, inv_ex2)
        info2 = subset_RS4_rows(indata["peak_info"][2], inv_ex2)

        ex2 = r("function(p3, ov23) 1:length(p3) %in% ov23[,2]")(p3, ov23)
        inv_ex3 = r["!"](ex3)
        p3 = subset_RS4_rows(p3, inv_ex3)
        info3 = subset_RS4_rows(indata["peak_info"][3], inv_ex3)

    peak_info = r["rbind"](indata["peak_info"][1], info2, info3)
    nubmer_peaks = r["nrow"](peak_info)[0]

    return peak_info
