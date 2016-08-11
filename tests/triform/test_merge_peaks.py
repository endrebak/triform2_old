import pytest

import pandas as pd
import numpy as np

from io import StringIO


@pytest.fixture
def input_data():
    s = StringIO(
        u"""NLP  MAX.NLP      LOC WIDTH    START      END CVG  SURL SURR FORM
chrY:57420858-57420942:1 3.923426 6.317243 57420916    85 57420858 57420942 10 0 2  1
chrY:10636518-10636576:2 2.869699 2.869699 10636545    59 10636518 10636576 9 0 16 2
chrY:10636618-10636676:3 2.869699 2.869699 10636645    59 10636618 10636676 9 0 16 3
chrY:10631247-10631396:3 4.039029 4.039029 10631334   150 10631247 10631396 14 24 0 3
chrY:10631347-10631496:2 4.039029 4.039029 10631334   150 10631347 10631496 14 24 0 2""")
    df = pd.read_table(s, sep="\s+", header=0)

    return df


@pytest.fixture
def expected_result():
    s = StringIO(
        u"""NLP  MAX.NLP      LOC WIDTH    START      END CVG  SURL SURR FORM
chrY:57420858-57420942:1 3.923426 6.317243 57420916    85 57420858 57420942 10 0 2  1
chrY:10636518-10636676:2 2.869699 2.869699 10636545,10636645    59 10636518 10636676 9 0 16 2
chrY:10631247-10631396:3 4.039029 4.039029 10631334   150 10631247 10631396 14 24 0 3
chrY:10631347-10631496:2 4.039029 4.039029 10631334   150 10631347 10631496 14 24 0 2""")
    df = pd.read_table(s, sep="\s+", header=0)

    return df


def merge_peaks(df):
    """Merge type2 and type3 peaks."""

    from triform.helper_functions import _locs_from_df, _create_intervaltree

    df2 = df.loc[df.FORM == 2]
    df3 = df.loc[df.FORM == 3]
    locs2 = _locs_from_df(df2)
    locs3 = _locs_from_df(df3)

    it2 = _create_intervaltree(locs2)

    pairs_to_merge = []
    for k, (start, end) in locs3.iterrows():

        intervals = it2.find(start, end)
        if intervals:
            pairs_to_merge.append([intervals[0], k])

    for pair in pairs_to_merge:
        indexes = pd.DataFrame(index=pd.Series(pair)).index
        pair_df = df.ix[indexes]
        start, end,


@pytest.mark.current
def test_merge_peaks(input_data, expected_result):

    result = merge_peaks(input_data)
    print(result)
    assert 0
    assert result == expected_result
