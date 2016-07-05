import pytest
from numpy import allclose, int32
import pandas as pd
from io import StringIO

from rpy2 import robjects as ro
from rpy2.robjects import r, pandas2ri
ri2py = pandas2ri.ri2py
from rpy2.robjects.packages import importr
importr("GenomicRanges")

from triform.preprocess.make_chromosome_cover_files import (
    make_chromosome_cover_files)
from triform.helper_functions import df_to_rle, rle_to_df

# @pytest.fixture
# def input_data():
#     """Creates a GRanges object from the df below."""
#     df = """seqnames start end strand
# chrY 10641380 10641405 +
# chrY 12314848 12314873 +
# chrY 11743048 11743073 -
# chrY 57396702 57396727 +
# chrY 11778685 11778710 +
# chrY 10595070 10595095 +
# chrY 16891566 16891591 +
# chrY 11936613 11936638 -
# chrY 11924977 11925002 +
# chrY 11764681 11764706 -"""

#     command = """
# lines = "{}"
# options(stringsAsFactors=FALSE)
# con <- textConnection(lines)
# df <- read.csv(con, colClasses=c("character", rep("integer",2), "character"), header=1, sep=" ")
# close(con)
# colnames(df) <- c("seqnames", "start", "end", "strand")
#     """.format(df)
#     r(command)
#     granges = r("makeGRangesFromDataFrame(df)")

#     return granges


@pytest.fixture
def expected_result_control(chromosome_cover_results_input):
    control = {k: df_to_rle(r["read.table"](v,
                                            sep=" "))
               for (k, v) in chromosome_cover_results_input.items()}
    return control


@pytest.fixture
def expected_result_treatment(chromosome_cover_results_chip):
    treatment = {k: df_to_rle(r["read.table"](v,
                                              sep=" "))
                 for (k, v) in chromosome_cover_results_chip.items()}
    return treatment


@pytest.fixture
def indata_chip(indata_make_chromosome_cover_files_chip, args):

    d = {}
    for k, v in indata_make_chromosome_cover_files_chip.items():
        df = r["read.table"](v, sep=" ")
        # print(r["head"](df))
        rd = r["as"](df, "GRanges")
        d[k] = rd

    return d


@pytest.mark.current
def test_make_chromosome_cover_files_chip(indata_chip,
                                          expected_result_treatment, args):
    result = make_chromosome_cover_files(indata_chip, args)
    print(result.keys())
    print(expected_result_treatment.keys())

    # print("result")
    # print(result["rep1", "forward"])
    # print("expected_result")
    # print(expected_result_treatment["rep1", "forward"])

    all_equals = r["isTRUE"](r["all.equal"](expected_result_treatment, result))
    assert all_equals
