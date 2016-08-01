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

    return pd.read_table(
        StringIO(u"""Chromosome     Start       End  Count
0        chrY  10631247  10631295      5
1        chrY  10631295  10631296      8
2        chrY  10631296  10631297      9
3        chrY  10631297  10631298     10
4        chrY  10631298  10631300     13
5        chrY  10631300  10631331     14
6        chrY  10631331  10631332     13
7        chrY  10631332  10631336     12
8        chrY  10631336  10631337     11
9        chrY  10631337  10631339     10
10       chrY  10631339  10631395      9
11       chrY  10631395  10631396      6
12       chrY  10631396  10631397      5
13       chrY  10636518  10636526      8
14       chrY  10636526  10636561      9
15       chrY  10636561  10636563      7
16       chrY  10636563  10636564      6
17       chrY  10636564  10636566      5
18       chrY  10636566  10636567      3
19       chrY  10636567  10636576      2
20       chrY  57420858  57420865      3
21       chrY  57420865  57420925      4
22       chrY  57420925  57420929      6
23       chrY  57420929  57420931      7
24       chrY  57420931  57420932      8
25       chrY  57420932  57420934      9
26       chrY  57420934  57420939     10
27       chrY  57420939  57420942      9"""),
        sep="\s+",
        header=0)


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

    results = create_bigwig(chip_data, fdr_table, args)
    # print(results, "results")
    # print(expected_result, "expected_result")

    assert results.equals(expected_result)


@pytest.fixture
def expected_result2():

    return pd.read_table(
        StringIO(u"""Chromosome     Start       End  Count
0       chrY  10631247  10631396      1
1       chrY  10631396  10631397      2"""),
        sep="\s+",
        header=0)


@pytest.fixture
def fdr_table2():
    data = u"""QVAL	NLQ	NLP	MAX.NLP	LOC	WIDTH	START	END	CVG	SURL	SURR	FORM
chrY:10631247-10631396:3	0.000178922590521532	3.74733482254352	4.03902854326486	4.03902854326486	10631334	150	10631247	10631396	14	24	0	3
"""

    #chrY:10636518-10636576:2	0.00247481305798851	2.60645760115479	2.86969903592937	2.86969903592937	10636545	59	10636518	10636576	9	0	16	2

    df = pd.read_table(StringIO(data), index_col=0)

    return df


@pytest.fixture
def chip_data2(run_length_encodings_full):

    df = """"lengths" "values"
"1" 10631246 0
"2" 148 1
"3" 1 2
"4" 1 1
"""

    # "5" 5122 0
    # "6" 58 1

    string_to_df = r("""function(x) {
    con <- textConnection(x)
    data <- read.csv(con, sep=" ", header=1)
    close(con)
    data
    }""")

    df = string_to_df(df)
    rle = df_to_rle(df)

    return {"chrY": {(1, 2, "center"): rle}}


@pytest.mark.current
def test_create_bigwig2(chip_data2, fdr_table2, expected_result2, args):

    results = create_bigwig(chip_data2, fdr_table2, args)
    print("results")
    print(results)
    print("expected_result2")
    print(expected_result2)

    assert 0, "Result is one too long."
    assert results.equals(expected_result2)
