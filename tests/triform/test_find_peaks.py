from collections import defaultdict
from itertools import product

import pandas as pd
import numpy as np

import pytest

from rpy2.robjects import r, pandas2ri
ri2py = pandas2ri.ri2py
from rpy2.robjects.packages import importr
importr("GenomicRanges")
importr("S4Vectors")
from triform.find_peaks import _find_peaks
from triform.helper_functions import (
    subset_RS4, subset_RS4_rows, subset_RS4_cols, df_to_iranges, df_to_rle)


@pytest.fixture
def input_data():

    results = defaultdict(lambda: defaultdict(dict))
    for peak_type, direction in product([1, 2, 3], ["reverse", "forward"]):

        df = r["read.table"]("tests/test_data/chromosome_result_%s%s.csv" % (
            direction, peak_type),
                             sep=" ")
        df = r('function(df) df[c("PEAK.START", "PEAK.END", "PEAK.WIDTH")]')(
            df)
        iranges = df_to_iranges(df)

        results["chrY", direction]["peaks"][peak_type] = iranges

    for direction in "reverse forward".split():
        df = r["read.table"]("tests/test_data/chromosome_cvg_%s.csv" %
                             direction)
        results["chrY", direction]["cvg"] = df_to_rle(df)

    for i in range(1, 4):
        results["chrY", "forward"]["peak_info"][i] = r["read.table"](
            "tests/test_data/chromosome_result_forward%s.csv" % i,
            sep=" ")
        results["chrY", "reverse"]["peak_info"][i] = r["read.table"](
            "tests/test_data/chromosome_result_reverse%s.csv" % i,
            sep=" ")

    return results


@pytest.fixture
def expected_result():

    results = {}
    for peak_type in range(1, 4):

        info = pd.read_table("tests/test_data/find_peaks_result_%s.csv" %
                             peak_type,
                             sep=" ")
        results["info", peak_type] = info

        peaks = pd.read_table("tests/test_data/merge_peaks_%s.csv" % peak_type,
                              sep=" ")
        results["peaks", peak_type] = peaks
    return results


@pytest.mark.unit
def test_find_peaks(input_data, expected_result, args):
    peaks, info = _find_peaks(input_data["chrY", "forward"],
                              input_data["chrY", "reverse"], "chrY", args)

    for i in range(1, 4):
        print(i)
        print("expected_result['info', i]")
        expected_result_py = expected_result["info", i]
        print(expected_result_py)
        print("info[i]")
        actual_result_py = ri2py(info[i])
        print(actual_result_py)
        assert np.allclose(expected_result_py, actual_result_py)

    for i in range(1, 4):
        print(i)
        print("expected_result['peaks', i]")
        expected_result_py = expected_result["peaks", i]
        print(expected_result_py)
        print("peaks[i]")
        actual_result_py = ri2py(r["as.data.frame"](peaks[i]))
        print(actual_result_py)
        assert np.allclose(expected_result_py, actual_result_py)
