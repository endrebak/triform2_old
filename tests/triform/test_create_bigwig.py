from collections import defaultdict

try:
    from itertools import izip_longest
except ImportError:
    from itertools import zip_longest

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
0        chrY  10631247  10631248      5
1        chrY  10631248  10631249      8
2        chrY  10631249  10631250      9
3        chrY  10631250  10631252     10
4        chrY  10631252  10631283     13
5        chrY  10631283  10631284     14
6        chrY  10631284  10631288     13
7        chrY  10631288  10631289     12
8        chrY  10631289  10631291     11
9        chrY  10631291  10631347     10
10       chrY  10631347  10631348      9
11       chrY  10631348  10631349      6
12       chrY  10631349  10631350      5
13       chrY  10631350  10631396      4
14       chrY  10636518  10636553      8
15       chrY  10636553  10636555      9
16       chrY  10636555  10636556      7
17       chrY  10636556  10636558      6
18       chrY  10636558  10636559      5
19       chrY  10636559  10636570      3
20       chrY  10636570  10636576      2
21       chrY  57420858  57420918      3
22       chrY  57420918  57420922      4
23       chrY  57420922  57420924      6
24       chrY  57420924  57420925      7
25       chrY  57420925  57420927      8
26       chrY  57420927  57420932      9
27       chrY  57420932  57420937     10
28       chrY  57420937  57420942      9"""),
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


@pytest.mark.unit
def test_create_bigwig(chip_data, fdr_table, expected_result, args):

    results = create_bigwig(chip_data, fdr_table, args)
    print(results, "results")
    print(expected_result, "expected_result")

    assert results.equals(expected_result)

# @pytest.fixture
# def expected_result2():

#     return pd.read_table(
#         StringIO(u"""Chromosome     Start       End  Count
# 0       chrY  10631247  10631395      1
# 1       chrY  10631395  10631396      2"""),
#         sep="\s+",
#         header=0)

# @pytest.fixture
# def fdr_table2():
#     data = u"""QVAL	NLQ	NLP	MAX.NLP	LOC	WIDTH	START	END	CVG	SURL	SURR	FORM
# chrY:10631247-10631396:3	0.000178922590521532	3.74733482254352	4.03902854326486	4.03902854326486	10631334	150	10631247	10631396	14	24	0	3
# """

#     #chrY:10636518-10636576:2	0.00247481305798851	2.60645760115479	2.86969903592937	2.86969903592937	10636545	59	10636518	10636576	9	0	16	2

#     df = pd.read_table(StringIO(data), index_col=0)

#     return df

# @pytest.fixture
# def chip_data2(run_length_encodings_full):

#     df = """"lengths" "values"
# "1" 10631246 0
# "2" 147 1
# "3" 1 2
# "4" 1 1
# """

#     # "5" 5122 0
#     # "6" 58 1

#     string_to_df = r("""function(x) {
#     con <- textConnection(x)
#     data <- read.csv(con, sep=" ", header=1)
#     close(con)
#     data
#     }""")

#     df = string_to_df(df)
#     rle = df_to_rle(df)

#     return {"chrY": {(1, 2, "center"): rle}}

# @pytest.mark.current
# def test_create_bigwig2(chip_data2, fdr_table2, expected_result2, args):

#     results = create_bigwig(chip_data2, fdr_table2, args)
#     print("results")
#     print(results)
#     print("expected_result2")
#     print(expected_result2)

#     assert results.equals(expected_result2)

# @pytest.fixture
# def expected_result3():

#     return pd.read_table(
#         StringIO(u"""Chromosome     Start       End  Count
# 0       chr16  15094395  15094426      2
# 1       chr16  15094426  15094582      1"""),
#         sep="\s+",
#         header=0)

# 15094582 - 15094395

# @pytest.fixture
# def fdr_table3():
#     data = u"""QVAL	NLQ	NLP	MAX.NLP	LOC	WIDTH	START	END	CVG	SURL	SURR	FORM
# chr16:15094426-15094582:1 0.00879484459358 2.0557718302 3.42192088586 22.2407724148 15094509 157 15094426 15094582 27.0 16.0 5.0 1
# # chr16:15094395-15094426:2 0.0090164769225 2.04496312454 3.40389080362 4.26957638035 15094411 32 15094395 15094426 14.0 1.0 16.0 2"""

#     df = pd.read_table(StringIO(data), index_col=0)

#     return df

# @pytest.fixture
# def chip_data3(run_length_encodings_full):

#     df = """"lengths" "values"
# "1" 15094394 0
# "2" 31 2
# "3" 156 1
# """

#     # "5" 5122 0
#     # "6" 58 1

#     string_to_df = r("""function(x) {
#     con <- textConnection(x)
#     data <- read.csv(con, sep=" ", header=1)
#     close(con)
#     data
#     }""")

#     df = string_to_df(df)
#     rle = df_to_rle(df)

#     return {"chr16": {(1, 2, "center"): rle}}

# @pytest.mark.current
# def test_create_bigwig3(chip_data3, fdr_table3, expected_result3, args):

#     results = create_bigwig(chip_data3, fdr_table3, args)
#     print("results")
#     print(results)
#     print("expected_result3")
#     print(expected_result3)

#     assert results.equals(expected_result3)
