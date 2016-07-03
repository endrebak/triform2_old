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
from triform.init import _init
from triform.chromosome import chromosome
from triform.helper_functions import df_to_rle, rle_to_df


@pytest.fixture
def expected_result():
    path = "tests/test_results/result_df.csv"
    return pd.read_table(path, sep=" ")


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


@pytest.mark.unit
def test_chromosome(chip_data, input_data, chip_sizes, input_sizes, args,
                    expected_result):

    input_data = input_data["reverse"]
    input_size = sum(input_sizes["reverse"].values())
    chip_sizes = chip_sizes["reverse"]

    chip_data = chip_data["reverse"]

    result = chromosome(chip_data, input_data, chip_sizes, input_size, args)

    result_df = result["peak_info"][1]
    result_df = ri2py(r["as.data.frame"](result_df)).reset_index(drop=True)

    for (_, prow), (_, xrow) in izip_longest(result_df.iterrows(),
                                             expected_result.iterrows()):
        if not np.allclose(prow, xrow):
            print(prow, xrow)
            raise

    assert np.allclose(result_df, expected_result)
