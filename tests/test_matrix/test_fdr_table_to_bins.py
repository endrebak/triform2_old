import pytest

import pandas as pd
import numpy as np

from io import StringIO
from collections import defaultdict

from triform.helper_functions import df_to_rle, df_to_iranges_end, rle_to_df, granges_to_bed_df, bed_df_to_granges
from triform.matrix.fdr_table_to_bins import fdr_table_to_bins
from rpy2.robjects import r


@pytest.fixture
def expected_result():
    d = defaultdict(
        lambda: pd.DataFrame(columns="Chromosome     Start       End   Enriched".split()),
        {"chrY": pd.read_table(
            StringIO("""Chromosome     Start       End   Enriched
0        chrY  57420851  57420860 1
1        chrY  57420861  57420870 1
2        chrY  57420871  57420880 1
3        chrY  57420881  57420890 1
4        chrY  57420891  57420900 1
5        chrY  57420901  57420910 1
6        chrY  57420911  57420920 1
7        chrY  57420921  57420930 1
8        chrY  57420931  57420940 1
9        chrY  57420941  57420950 1
10       chrY  10631241  10631250 1
11       chrY  10631251  10631260 1
12       chrY  10631261  10631270 1
13       chrY  10631271  10631280 1
14       chrY  10631281  10631290 1
15       chrY  10631291  10631300 1
16       chrY  10631301  10631310 1
17       chrY  10631311  10631320 1
18       chrY  10631321  10631330 1
19       chrY  10631331  10631340 1
20       chrY  10631341  10631350 1
21       chrY  10631351  10631360 1
22       chrY  10631361  10631370 1
23       chrY  10631371  10631380 1
24       chrY  10631381  10631390 1
25       chrY  10631391  10631400 1"""),
            sep="\s+",
            usecols=[1, 2, 3, 4],
            index_col=[0, 1, 2]),
         "chrX": pd.read_table(
             StringIO("""Chromosome     Start       End   Enriched
26       chrX  10636511  10636520 1
27       chrX  10636521  10636530 1
28       chrX  10636531  10636540 1
29       chrX  10636541  10636550 1
30       chrX  10636551  10636560 1
31       chrX  10636561  10636570 1
32       chrX  10636571  10636580  1"""),
             sep="\s+",
             usecols=[1, 2, 3, 4],
             index_col=[0, 1, 2])})
    return d


@pytest.fixture
def input_data():
    data = u"""CHROMOSOME  QVAL	NLQ	NLP	MAX.NLP	LOC	WIDTH	START	END	CVG	SURL	SURR	FORM
chrY 	0.000178922590521532	3.74733482254352	3.9234260815992	6.31724273434332	57420916	85	57420858	57420942	10	0	2	1
chrY 	0.000178922590521532	3.74733482254352	4.03902854326486	4.03902854326486	10631334	150	10631247	10631396	14	24	0	3
chrX 	0.00247481305798851	2.60645760115479	2.86969903592937	2.86969903592937	10636545	59	10636518	10636576	9	0	16	2"""

    df = pd.read_table(
        StringIO(data),
        index_col=["CHROMOSOME", "START", "END"],
        sep="\s+")

    return df


@pytest.mark.unit
def test_fdr_table_to_bins(input_data, expected_result, args):

    result = fdr_table_to_bins(input_data, 10)
    print(result, "result")
    print(expected_result, "expected_result")
    # print(result.xs("chrY", drop_level=False), "xs " * 100)
    assert result["chrY"].equals(expected_result["chrY"])

    # granges = _find_read_midpoints(input_data, "chrY", "counts")

    # print("granges")
    # print(r["head"](granges))
    # print(r["seqlengths"](granges))
    # print("expected_result")
    # print(r["head"](expected_result))

    # print(r["all.equal"](granges, expected_result))
    # # print(r["identical"](granges, expected_result))
    # # assert r["all.equal"](granges, expected_result)
    # assert r["all.equal"](granges, expected_result)[0]
