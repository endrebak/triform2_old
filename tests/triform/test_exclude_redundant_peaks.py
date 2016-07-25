from collections import defaultdict
from itertools import product

import pandas as pd
import numpy as np

import pytest
from io import StringIO

from rpy2.robjects import r, pandas2ri
ri2py = pandas2ri.ri2py
from rpy2.robjects.packages import importr
importr("GenomicRanges")
importr("S4Vectors")
from triform.chromosome import chromosome
from triform.helper_functions import (
    subset_RS4, subset_RS4_rows, subset_RS4_cols, df_to_iranges, df_to_rle)


@pytest.fixture
def expected_result():
    s = StringIO(
        u"""NLP  MAX.NLP      LOC WIDTH    START      END CVG  SURL SURR FORM
chrY:57420858-57420942:1 3.923426 6.317243 57420916    85 57420858 57420942 10 0 2  1
chrY:10636518-10636576:2 2.869699 2.869699 10636545    59 10636518 10636576 9 0 16 2
chrY:10631247-10631396:3 4.039029 4.039029 10631334   150 10631247 10631396 14 24 0 3""")
    df = pd.read_table(s, sep="\s+", header=0)

    return df


@pytest.fixture
def input_data():

    results = {}
    for peak_type in range(1, 4):

        info = r["read.table"]("tests/test_data/find_peaks_result_%s.csv" %
                               peak_type,
                               sep=" ")
        results["info", peak_type] = info

        peaks = r["read.table"]("tests/test_data/merge_peaks_%s.csv" %
                                peak_type,
                                sep=" ")
        results["peaks", peak_type] = df_to_iranges(peaks)
    return results


@pytest.mark.current
def test_exclude_redundant_peaks(input_data, expected_result):

    result = exclude_redundant_peaks(input_data)

    print("result")
    print(result)
    print("expected_result")
    print(expected_result)
    assert np.allclose(expected_result, ri2py(result))


def exclude_redundant_peaks(indata):

    print(indata.keys())
    p1 = indata["peaks", 1]
    p2 = indata["peaks", 2]
    p3 = indata["peaks", 3]

    ov12 = r("function(p1, p2) matrix(as.matrix(findOverlaps(p1,p2)),ncol=2)")(
        p1, p2)
    if r["nrow"](ov12):
        ex2 = r("function(p2, ov12) 1:length(p2) %in% ov12[,2]")(p2, ov12)
        inv_ex2 = r["!"](ex2)
        p2 = subset_RS4_rows(p2, inv_ex2)
        info2 = subset_RS4_rows(indata["info", 2], inv_ex2)

    ov13 = r("function(p1, p3) matrix(as.matrix(findOverlaps(p1,p3)),ncol=2)")(
        p1, p3)
    if r["nrow"](ov13):
        ex3 = r("function(p3, ov13) 1:length(p3) %in% ov13[,2]")(p3, ov13)
        inv_ex3 = r["!"](ex3)
        p3 = subset_RS4_rows(p3, inv_ex3)
        info3 = subset_RS4_rows(indata["info", 3], inv_ex3)

    ov23 = r(
        "function(p2, p3) matrix(as.matrix(findOverlaps(p2,p3, maxgap=100)),ncol=2)")(
            p1, p3)
    if r["nrow"](ov23):
        ex2 = r("function(p2, ov23) 1:length(p2) %in% ov23[,1]")(p2, ov23)
        inv_ex2 = r["!"](ex2)
        p2 = subset_RS4_rows(p2, inv_ex2)
        info2 = subset_RS4_rows(indata["info", 2], inv_ex2)

        ex2 = r("function(p3, ov23) 1:length(p3) %in% ov23[,2]")(p3, ov23)
        inv_ex3 = r["!"](ex3)
        p3 = subset_RS4_rows(p3, inv_ex3)
        info3 = subset_RS4_rows(indata["info", 3], inv_ex3)

    peak_info = r["rbind"](indata["info", 1], info2, info3)
    nubmer_peaks = r["nrow"](peak_info)[0]

    return peak_info
