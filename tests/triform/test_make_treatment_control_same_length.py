import pytest

import pandas as pd
import numpy as np

from io import StringIO

from triform.helper_functions import df_to_rle, rle_to_pandas_df
from rpy2.robjects import r
from triform.make_treatment_control_same_length import (
    make_treatment_control_same_length)


@pytest.fixture
def input_data(init_result_chip, init_result_input):

    chip_data, input_data = dict(), dict()

    chip_data["chrY"] = {k: df_to_rle(r["read.table"](v,
                                                      sep=" "))
                         for k, v in init_result_chip.items()}

    input_data["chrY"] = {k: df_to_rle(r["read.table"](v,
                                                       sep=" "))
                          for k, v in init_result_input.items()}

    return chip_data, input_data


@pytest.fixture
def expected_result(init_result_chip_plus3, init_result_input_plus3):

    chip_data, input_data = dict(), dict()

    chip_data = {k: pd.read_table(v,
                                  sep=" ",
                                  index_col=0)
                 for k, v in init_result_chip_plus3.items()}

    input_data = {k: pd.read_table(v,
                                   sep=" ",
                                   index_col=0)
                  for k, v in init_result_input_plus3.items()}

    return chip_data, input_data


@pytest.mark.unit
def test_make_treatment_control_same_length(input_data, expected_result):

    treatment, control = input_data
    x_treatment, x_control = expected_result

    treatment_result, control_result = make_treatment_control_same_length(
        treatment, control)

    assertions = []
    for k, v in treatment_result["chrY"].items():
        v = rle_to_pandas_df(v)
        x = x_treatment[k]
        assertions.append(v.equals(x))

    assert all(assertions)
