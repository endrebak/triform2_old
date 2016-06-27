import pytest
from glob import glob

from collections import namedtuple

import pandas as pd

__author__ = "Endre Bakken Stovner https://github.com/endrebak/"
__license__ = "MIT"

MockNamespace = namedtuple("MockNamespace",
                           ["number_cores", "genome", "read_width",
                            "flank_distance", "treatment", "control"])

# egs = 2290813547.4  # this is the effective genome size used by the original sicer for hg19

# @pytest.fixture(scope="session")
# def args_200_fast():
#     return MockNamespace(25, "hg19", False, 200, 150, False, 3, 1, egs, False,
#                          ["examples/test.bed"], ["examples/control.bed"])

# @pytest.fixture(scope="session")
# def args_200():
#     return MockNamespace(1, "hg19", False, 200, 150, False, 3, 0.05, egs,
#                          False, ["examples/test.bed"],
#                          ["examples/control.bed"])


@pytest.fixture(scope="session")
def args_tests():
    return MockNamespace(2, "hg19", 100, 150, ["examples/test.bed"],
                         ["examples/control.bed"])


@pytest.fixture(scope="session")
def multiple_chromosomes():
    return pd.read_table("examples/multiple_chromos.bed", header=None)


@pytest.fixture(scope="session")
def run_length_encodings():
    return {"FORWARD.CENTER": "examples/rle/forward_center.csv",
            "FORWARD.LEFT": "examples/rle/forward_left.csv",
            "FORWARD.RIGHT": "examples/rle/forward_right.csv",
            "REVERSE.CENTER": "examples/rle/reverse_center.csv",
            "REVERSE.LEFT": "examples/rle/reverse_left.csv",
            "REVERSE.RIGHT": "examples/rle/reverse_right.csv"}
