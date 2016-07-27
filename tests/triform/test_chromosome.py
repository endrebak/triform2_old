from collections import defaultdict
from itertools import izip_longest

import pandas as pd
import numpy as np

import pytest

from rpy2.robjects import r, pandas2ri
ri2py = pandas2ri.ri2py
from rpy2.robjects.packages import importr
importr("GenomicRanges")
importr("S4Vectors")
from triform.chromosome import chromosome
from triform.helper_functions import df_to_rle, rle_to_df


@pytest.fixture
def expected_result():
    results = defaultdict(dict)
    for i in range(1, 4):
        results["chrY", "forward", i] = pd.read_table(
            "tests/test_data/chromosome_result_forward%s.csv" % i,
            sep=" ")
        results["chrY", "reverse", i] = pd.read_table(
            "tests/test_data/chromosome_result_reverse%s.csv" % i,
            sep=" ")
    return results


@pytest.fixture
def expected_result_peaks_zscores():
    _peaks = [pd.read_table(path,
                            sep=" ")
              for path in [
                  "tests/test_results/p1.csv", "tests/test_results/p2.csv",
                  "tests/test_results/p3.csv"
              ]]
    _zscores = [pd.read_table(path,
                              sep=" ")
                for path in [
                    "tests/test_results/z1.csv", "tests/test_results/z2.csv",
                    "tests/test_results/z3.csv"
                ]]
    return _peaks, _zscores


@pytest.fixture
def input_sizes(input_sizes):
    return {"chrY": input_sizes}


@pytest.fixture
def chip_sizes(chip_sizes):
    return {"chrY": chip_sizes}


@pytest.fixture
def input_data(run_length_encodings_full_backgr):

    return {"chrY": {k: df_to_rle(r["read.table"](v,
                                                  sep=" "))
                     for (k, v) in run_length_encodings_full_backgr.items()}}


@pytest.fixture
def chip_data(run_length_encodings_full):

    return {"chrY": {k: df_to_rle(r["read.table"](v,
                                                  sep=" "))
                     for (k, v) in run_length_encodings_full.items()}}


@pytest.mark.unit
def test_chromosome(chip_data, input_data, chip_sizes, input_sizes, args,
                    expected_result):

    results = chromosome(chip_data, input_data, chip_sizes, input_sizes, args)

    # results["chrY", "forward"]["peak_info"][1]
    # results["chrY", "forward"]["peak_info"][2]
    # results["chrY", "forward"]["peak_info"][3]

    asserts = []
    for (_chromosome, direction,
         peak), expected_result in expected_result.items():
        print(_chromosome, direction, peak)
        actual_result = ri2py(results[_chromosome, direction]["peak_info"][
            peak])
        # print(actual_result.head(10), "actual_result")
        # print(expected_result.head(10), "expected_result")

        # print(actual_result.tail(), "actual_result")
        # print(expected_result.tail(), "expected_result")

        # print(actual_result.dtypes, "actual_result")
        # print(expected_result.dtypes, "expected_result")

        # actual_result.to_csv("actual.csv", sep=" ")
        # Using allclose here, since it does not care about 32/64 differences
        assert np.allclose(expected_result, actual_result)

    assert all(asserts)
