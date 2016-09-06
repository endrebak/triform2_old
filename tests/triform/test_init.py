## Use function below to create csv data of all the files! Will need it to crate
# test data for the rest of the sotware

# want to convert list of 6 chrcovers to csv so can store for tests

import pandas as pd
import numpy as np
import pytest

from rpy2.robjects import r, pandas2ri
ri2py = pandas2ri.ri2py
from rpy2.robjects.packages import importr
importr("GenomicRanges")
importr("S4Vectors")

from triform.init import _init_background, _init_treatment
from triform.helper_functions import df_to_rle, rle_to_df


@pytest.fixture
def expected_result_chip(init_result_chip):
    return {k: pd.read_table(v,
                             sep=" ",
                             index_col=0)
            for k, v in init_result_chip.items()}


@pytest.fixture
def expected_result_input(init_result_input):
    return {k: pd.read_table(v,
                             sep=" ",
                             index_col=0)
            for k, v in init_result_input.items()}


@pytest.mark.unit
def test_init_treatment(expected_result_chip, input_data_treatment, args):

    asserts = []
    results_treatment = _init_treatment(input_data_treatment, args)
    for k, v in results_treatment.items():
        print(k)
        expected_result = expected_result_chip[k]
        actual_result = ri2py(rle_to_df(v)).astype(np.int64)
        print(expected_result.tail())
        print(actual_result.tail())
        assertion = np.allclose(expected_result, actual_result)
        print(assertion)
        asserts.append(assertion)

    print(asserts)
    assert all(asserts)


@pytest.mark.unit
def test_init_background(expected_result_input, input_data_control, args):

    results_control = _init_background(input_data_control, args)

    asserts = []
    for k, v in results_control.items():
        expected_result = expected_result_input[k]
        actual_result = ri2py(rle_to_df(v)).astype(np.int64)
        asserts.append(np.allclose(expected_result, actual_result))

    print(asserts)

    assert all(asserts)


@pytest.fixture
def input_data_control(chromosome_cover_results_input):
    control = {k: df_to_rle(r["read.table"](v,
                                            sep=" "))
               for (k, v) in chromosome_cover_results_input.items()}
    return control


@pytest.fixture
def input_data_treatment(chromosome_cover_results_chip):
    treatment = {k: df_to_rle(r["read.table"](v,
                                              sep=" "))
                 for (k, v) in chromosome_cover_results_chip.items()}
    return treatment
