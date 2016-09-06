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
