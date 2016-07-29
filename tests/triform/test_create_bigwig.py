from collections import defaultdict
from itertools import izip_longest
from io import StringIO

import pandas as pd
import numpy as np

import pytest

from rpy2.robjects import r, pandas2ri
ri2py = pandas2ri.ri2py
from rpy2.robjects.packages import importr
importr("GenomicRanges")
importr("S4Vectors")
from triform.create_bigwig import create_bigwig
from triform.helper_functions import df_to_rle, rle_to_df


@pytest.fixture
def expected_result():
    pass


@pytest.fixture
def fdr_table():
    data = u"""QVAL	NLQ	NLP	MAX.NLP	LOC	WIDTH	START	END	CVG	SURL	SURR	FORM
chrY:57420858-57420942:1	0.000178922590521532	3.74733482254352	3.9234260815992	6.31724273434332	57420916	85	57420858	57420942	10	0	2	1
chrY:10631247-10631396:3	0.000178922590521532	3.74733482254352	4.03902854326486	4.03902854326486	10631334	150	10631247	10631396	14	24	0	3
chrX:10631247-10631396:3	0.000178922590521532	3.74733482254352	4.03902854326486	4.03902854326486	10631334	150	10631247	10631396	14	24	0	3
chrY:10636518-10636576:2	0.00247481305798851	2.60645760115479	2.86969903592937	2.86969903592937	10636545	59	10636518	10636576	9	0	16	2"""

    df = pd.read_table(StringIO(data), index_col=0)

    return df


@pytest.fixture
def chip_data(run_length_encodings_full):

    return {"chrY": {k: df_to_rle(r["read.table"](v,
                                                  sep=" "))
                     for (k, v) in run_length_encodings_full.items()}}


@pytest.mark.current
def test_create_bigwig(chip_data, fdr_table, expected_result, args):

    # print(chip_data.keys())
    # print(chip_data.values())
    results = create_bigwig(chip_data, fdr_table)
    # args)
    assert 0
