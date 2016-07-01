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


def df_to_rle(df):

    _df_to_rle = r("""function(df){
    Rle(df$values, df$lengths)
    }
    """)

    return _df_to_rle(df)


def rle_to_df(rle):

    _rle_to_df = r("""function(rle){
    cbind(runValue(rle), runLength(rle))
    }
    """)

    return _rle_to_df(rle)


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


@pytest.mark.current
def test_chromosome(chip_data, input_data, chip_sizes, input_sizes, args,
                    expected_result):

    input_data = input_data["reverse"]
    input_size = sum(input_sizes["reverse"].values())
    chip_sizes = chip_sizes["reverse"]

    chip_data = chip_data["reverse"]

    result = chromosome(chip_data, input_data, chip_sizes, input_size, args)
    result_df = result["peak_info"][1]
    result_df = ri2py(r["as.data.frame"](result_df)).reset_index(drop=True)
    # result_df.to_csv("tests/test_results/result_df.csv", sep=" ", index=False)

    for (_, prow), (_, xrow) in izip_longest(result_df.iterrows(),
                                             expected_result.iterrows()):
        if not np.allclose(prow, xrow):
            print(prow, xrow)
            raise

    assert np.allclose(result_df, expected_result)

# @pytest.mark.current
# def test_compute_peaks_and_zscores(chip_data, input_data, chip_sizes,
#                                    input_sizes, args,
#                                    expected_result_peaks_zscores):
#     expected_peaks, expected_zscores = expected_result_peaks_zscores

#     input_data = input_data["reverse"]
#     input_size = sum(input_sizes["reverse"].values())
#     chip_sizes = chip_sizes["reverse"]

#     chip = chip_data["reverse"]

#     input = r["Reduce"]("+", input_data)

#     ratios = compute_ratios(chip_sizes, input_size)

#     ratio = input_size / sum(chip_sizes.values())

#     center = collect_key(chip, "center")
#     left = collect_key(chip, "left")
#     right = collect_key(chip, "right")

#     cvg = r["Reduce"]("+", list(center.values()))
#     left = r["Reduce"]("+", list(left.values()))
#     right = r["Reduce"]("+", list(right.values()))

#     _peaks, _zscores = compute_peaks_and_zscores(
#         cvg, center, left, right, chip, input, ratios, ratio, args)

#     results = []
#     for i, (p, x) in enumerate(zip(_peaks, expected_peaks)):
#         print(i, "peaks")
#         p = ri2py(r["as.data.frame"](p)).astype(np.int64)

#         for (_, prow), (_, xrow) in izip_longest(p.iterrows(), x.iterrows()):
#             if not prow.equals(xrow):
#                 print(prow, xrow)

#         results.append(p.equals(x))
#     print(results)
