import pytest

import pandas as pd
import numpy as np

from io import StringIO

from triform.helper_functions import df_to_rle, df_to_iranges_end, rle_to_df, granges_to_bed_df, bed_df_to_granges
from triform.matrix.find_read_midpoints import _find_read_midpoints
from rpy2.robjects import r


@pytest.fixture
def expected_result():
    df = r["read.table"]("tests/test_data/matrix_result.bed")
    gr = bed_df_to_granges(df)
    gr = r("function(gr) {seqlengths(gr) = 10689210; gr}")(gr)
    return gr


@pytest.fixture
def input_data():

    df = r["read.table"]("tests/test_data/Input_rep1.bed",
                         header=False,
                         sep="\t",
                         nrows=50)
    r["print"](df)
    create_ir_df = r("""function(df) {
    ir_df = data.frame(df[[2]], df[[3]])
    ir_df
    }
    """)

    ir_df = create_ir_df(df)
    print(ir_df)

    return df_to_iranges_end(ir_df)

# @pytest.fixture
# def enriched_bins():
#     return pd.read_table(
#         StringIO("""CHROMOSOME     START       END   ENRICHED
# chrY  10519611 10519620 1"""),
#         sep="\s+",
#         index_col=[0, 1, 2])


@pytest.mark.current
def test_create_matrix(input_data, expected_result):

    granges = _find_read_midpoints(input_data, "chrY", "counts")

    print("granges")
    print(r["head"](granges))
    print(r["seqlengths"](granges))
    print("expected_result")
    print(r["head"](expected_result))

    print(r["all.equal"](granges, expected_result))
    assert r["all.equal"](granges, expected_result)[0]
