import pytest

import pandas as pd
import numpy as np

from io import StringIO


@pytest.fixture
def expected_result():
    pass


def name():
    """Write description as command ending in a period."""


def test_name(args_end_to_end_one_file, expected_result):

    result = name(input_data)
    assert 0  # result == expected_result
