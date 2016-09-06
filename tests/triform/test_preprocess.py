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
    print(control.keys(), 'control.keys()')
    return control


@pytest.fixture
def expected_result_treatment(chromosome_cover_results_chip):
    treatment = {k: pd.read_table(v,
                                  sep=" ")
                 for (k, v) in chromosome_cover_results_chip.items()}
    return treatment


@pytest.mark.unit
def test_preprocess(args, expected_result_treatment, expected_result_control):
    treatment, control, treatment_sizes, control_sizes = preprocess(args)

    t = treatment["chrY"]
    for k, l in t.items():
        l_as_df = ri2py(rle_to_df(l)).astype(np.int64)

        assert np.allclose(expected_result_treatment[k], l_as_df)

    c = control["chrY"]
    print(expected_result_control.keys(), 'expected_result_control.keys()')
    for k, l in c.items():
        f, d = k
        print(f)
        print(d)
        print((f, d) in expected_result_control)
        l_as_df = ri2py(rle_to_df(l)).astype(np.int64)
        print(l_as_df.head(), 'l_as_df.head()')
        print(expected_result_control[k].head(),
              "expected_result_control[k].head()")

        assert np.allclose(expected_result_control[k], l_as_df)
