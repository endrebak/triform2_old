from collections import defaultdict
import pandas as pd
import numpy as np

import pytest
from io import StringIO

from rpy2.robjects import r, pandas2ri
ri2py = pandas2ri.ri2py
from triform.exclude_redundant_peaks import _exclude_redundant_peaks
from triform.helper_functions import (df_to_iranges)


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

    results = defaultdict(dict)
    for peak_type in range(1, 4):

        info = r["read.table"]("tests/test_data/find_peaks_result_%s.csv" %
                               peak_type,
                               sep=" ",
                               row_names=1)
        results["peak_info"][peak_type] = info

    return results


@pytest.mark.unit
def test_exclude_redundant_peaks(input_data, expected_result, args):

    print(input_data, "input_data")

    result = _exclude_redundant_peaks(input_data, args)

    print("result")
    print(result)
    # print(result.dtypes)
    print("expected_result")
    print(expected_result)
    assert np.allclose(expected_result, ri2py(result))


@pytest.fixture
def expected_result_missing3():
    s = StringIO(
        u"""NLP  MAX.NLP      LOC WIDTH    START      END CVG  SURL SURR FORM
chrY:57420858-57420942:1 3.923426 6.317243 57420916    85 57420858 57420942 10 0 2  1
chrY:10636518-10636576:2 2.869699 2.869699 10636545    59 10636518 10636576 9 0 16 2""")
    df = pd.read_table(s, sep="\s+", header=0)

    return df


@pytest.fixture
def input_data_missing3():

    results = defaultdict(dict)
    for peak_type in range(1, 3):

        info = r["read.table"]("tests/test_data/find_peaks_result_%s.csv" %
                               peak_type,
                               sep=" ",
                               row_names=1)
        results["peak_info"][peak_type] = info

    return results


@pytest.mark.unit
def test_exclude_redundant_peaks_missing3(input_data_missing3,
                                          expected_result_missing3, args):

    print(input_data_missing3, "input_data")

    result = _exclude_redundant_peaks(input_data_missing3, args)

    print("result")
    print(result)
    # print(result.dtypes)
    print("expected_result")
    print(expected_result_missing3)
    assert np.allclose(expected_result_missing3, ri2py(result))


@pytest.fixture
def expected_result_real_error():
    s = StringIO(
        u"""NLP  MAX.NLP      LOC WIDTH    START      END CVG  SURL SURR FORM
chrY:56845739-56845895:1" 4.03902854326486 4.03902854326486 56845812 157 56845739 56845895 7 0 0 1""")
    df = pd.read_table(s, sep="\s+", header=0)

    return df


@pytest.fixture
def input_data_only_real_error():

    results = defaultdict(dict)
    for peak_type in range(1, 4):

        info = r["read.table"]("tests/test_data/real_errors/%s_chrY_pre_exclude_peaks.csv" %
                               peak_type,
                               sep=" ",
                               row_names=1)
        results["peak_info"][peak_type] = info

    return results


@pytest.mark.unit
def test_exclude_redundant_peaks_real_error(input_data_only_real_error, expected_result_real_error, args):

    # print(input_data_only_real_error, "input_data")

    result = _exclude_redundant_peaks(input_data_only_real_error, args)

    print("result")
    print(result)
    # print(result.dtypes)
    print("expected_result")
    print(expected_result_real_error)
    # assert 0
    assert np.allclose(expected_result_real_error, ri2py(result))
