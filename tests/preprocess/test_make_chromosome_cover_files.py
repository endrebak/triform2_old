import sys
import pytest
from numpy import allclose, int64
import pandas as pd
from io import StringIO
from collections import defaultdict

from rpy2 import robjects as ro
from rpy2.robjects import r, pandas2ri
ri2py = pandas2ri.ri2py
from rpy2.robjects.packages import importr
importr("GenomicRanges")

from triform.preprocess.make_chromosome_cover_files import (
    make_chromosome_cover_files)
from triform.helper_functions import df_to_rle, rle_to_df


@pytest.fixture
def expected_result_control(chromosome_cover_results_input):
    control = {k: pd.read_table(v,
                                sep=" ")
               for (k, v) in chromosome_cover_results_input.items()}
    return control


@pytest.fixture
def expected_result_treatment(chromosome_cover_results_chip):
    treatment = {k: pd.read_table(v,
                                  sep=" ")
                 for (k, v) in chromosome_cover_results_chip.items()}
    return treatment


@pytest.fixture
def indata_chip(indata_make_chromosome_cover_files_chip, args):

    d = defaultdict(dict)
    for k, v in indata_make_chromosome_cover_files_chip.items():
        df = r["read.table"](
            v,
            sep=" ",
            col_names=["chromosome", "start", "end", "strand"])
        rd = r["as"](df, "GRanges")
        d["chrY"][k] = rd

    return d


@pytest.mark.unit
def test_make_chromosome_cover_files_chip(indata_chip,
                                          expected_result_treatment, args):
    cvgs, sizes = make_chromosome_cover_files(indata_chip, args)

    sys.stdout = sys.stderr

    print(expected_result_treatment.keys(), "expected_result_treatment")
    for (f, direction), v in expected_result_treatment.items():
        print(v.head(), "v")
        print(v.tail(), "v")
        print('cvgs[f]["chrY",' + direction + ']')
        actual = ri2py(rle_to_df(cvgs[f]["chrY", direction])).astype(int64)
        print("actual")
        print(actual)
        print("expected")
        print(v)
        assert allclose(actual, v)
