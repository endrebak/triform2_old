import pytest
from numpy import allclose, int32
import numpy as np

import pandas as pd
from io import StringIO

from rpy2 import robjects as ro
from rpy2.robjects import r, pandas2ri
ri2py = pandas2ri.ri2py
from rpy2.robjects.packages import importr
importr("GenomicRanges")

from triform.preprocess.preprocess import preprocess
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


@pytest.mark.unit
def test_preprocess(args, expected_result_treatment, expected_result_control):
    treatment, control = preprocess(args)
    for rep, l in zip(["rep1", "rep2"], control["chrY"]):
        cvg = l.rx2("CVG")
        get_strand = r("function(d, k) d[k]")

        pos = get_strand(cvg, "+")
        pos = ri2py(rle_to_df(pos)).astype(np.int64)

        neg = get_strand(cvg, "-")
        neg = ri2py(rle_to_df(neg)).astype(np.int64)

        assert expected_result_control[rep, "forward"].equals(pos)
        assert expected_result_control[rep, "reverse"].equals(neg)
