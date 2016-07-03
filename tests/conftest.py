import pytest
from glob import glob

from collections import namedtuple

import pandas as pd

__author__ = "Endre Bakken Stovner https://github.com/endrebak/"
__license__ = "MIT"

MockNamespace = namedtuple("MockNamespace",
                           ["number_cores", "genome", "read_width",
                            "flank_distance", "min_enrichment", "max_p",
                            "min_shift", "min_width", "treatment", "control"])
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
def args():
    return MockNamespace(2, "hg19", 100, 150, 0.375, 0.1, 10, 10,
                         ["examples/srf_huds_Gm12878_rep1.bed",
                          "examples/srf_huds_Gm12878_rep2.bed"],
                         ["examples/backgr_huds_gm12878_rep1.bed",
                          "examples/backgr_huds_gm12878_rep2.bed"])


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


@pytest.fixture(scope="session")
def run_length_encodings_full_rep1():
    return {"forward": {"center": "examples/rle_full/rep1_forward_center.csv",
                        "left": "examples/rle_full/rep1_forward_left.csv",
                        "right": "examples/rle_full/rep1_forward_right.csv"},
            "reverse": {"center": "examples/rle_full/rep1_reverse_center.csv",
                        "left": "examples/rle_full/rep1_reverse_left.csv",
                        "right": "examples/rle_full/rep1_reverse_right.csv"}}


@pytest.fixture(scope="session")
def run_length_encodings_full_rep2():
    return {"forward": {"center": "examples/rle_full/rep2_forward_center.csv",
                        "left": "examples/rle_full/rep2_forward_left.csv",
                        "right": "examples/rle_full/rep2_forward_right.csv"},
            "reverse": {"center": "examples/rle_full/rep2_reverse_center.csv",
                        "left": "examples/rle_full/rep2_reverse_left.csv",
                        "right": "examples/rle_full/rep2_reverse_right.csv"}}


@pytest.fixture(scope="session")
def run_length_encodings_full_rep1_backgr():
    return {"forward": "examples/rle_full/backgr_rep1_forward_center.csv",
            "reverse": "examples/rle_full/backgr_rep1_reverse_center.csv"}


@pytest.fixture(scope="session")
def run_length_encodings_full_rep2_backgr():
    return {"forward": "examples/rle_full/backgr_rep2_forward_center.csv",
            "reverse": "examples/rle_full/backgr_rep2_reverse_center.csv"}


@pytest.fixture(scope="session")
def sizes_rep1():
    return {"forward": 8986, "reverse": 8688}


@pytest.fixture(scope="session")
def sizes_rep2():
    return {"forward": 8528, "reverse": 8334}


@pytest.fixture(scope="session")
def sizes_rep1_backgr():
    return {"forward": 7685, "reverse": 7798}


@pytest.fixture(scope="session")
def sizes_rep2_backgr():
    return {"forward": 8414, "reverse": 8214}

#    7685 examples/backgr_huds_Gm12878_rep1.bed:chrY
#    8414 examples/backgr_huds_Gm12878_rep2.bed:chrY
#    8986 examples/srf_huds_Gm12878_rep1.bed:chrY
#    8528 examples/srf_huds_Gm12878_rep2.bed:chrY
# labsenter@ans-180249 ~/h/c/triform> grep "-" examples/*huds* | cut -f 1 | sort | uniq -c
#    7798 examples/backgr_huds_Gm12878_rep1.bed:chrY
#    8214 examples/backgr_huds_Gm12878_rep2.bed:chrY
#    8688 examples/srf_huds_Gm12878_rep1.bed:chrY
#    8334 examples/srf_huds_Gm12878_rep2.bed:chrY


@pytest.fixture(scope="session")
def chromosome_cover_results_input():

    return {
        ("rep1", "reverse"):
        "tests/test_data/chrY_backgr_huds_Gm12878_rep1_neg.csv",
        ("rep1", "forward"):
        "tests/test_data/chrY_backgr_huds_Gm12878_rep1_pos.csv",
        ("rep2", "reverse"):
        "tests/test_data/chrY_backgr_huds_Gm12878_rep2_neg.csv",
        ("rep2", "forward"):
        "tests/test_data/chrY_backgr_huds_Gm12878_rep2_pos.csv",
    }


@pytest.fixture(scope="session")
def chromosome_cover_results_chip():
    return {
        ("rep1", "reverse"):
        "tests/test_data/chrY_srf_huds_Gm12878_rep1_neg.csv",
        ("rep1", "forward"):
        "tests/test_data/chrY_srf_huds_Gm12878_rep1_pos.csv",
        ("rep2", "reverse"):
        "tests/test_data/chrY_srf_huds_Gm12878_rep2_neg.csv",
        ("rep2", "forward"):
        "tests/test_data/chrY_srf_huds_Gm12878_rep2_pos.csv"
    }


@pytest.fixture(scope="session")
def init_result_input():
    return {
        ("rep1", "forward", "center"):
        "tests/test_data/backgr_huds_Gm12878_rep1.FORWARD.CENTER.csv",
        ("rep2", "forward", "center"):
        "tests/test_data/backgr_huds_Gm12878_rep2.FORWARD.CENTER.csv",
        ("rep1", "reverse", "center"):
        "tests/test_data/backgr_huds_Gm12878_rep1.REVERSE.CENTER.csv",
        ("rep2", "reverse", "center"):
        "tests/test_data/backgr_huds_Gm12878_rep2.REVERSE.CENTER.csv",
    }


@pytest.fixture(scope="session")
def init_result_chip():
    return {
        ("rep1", "forward", "center"):
        "tests/test_data/srf_huds_Gm12878_rep1.FORWARD.CENTER.csv",
        ("rep1", "forward", "left"):
        "tests/test_data/srf_huds_Gm12878_rep1.FORWARD.LEFT.csv",
        ("rep1", "forward", "right"):
        "tests/test_data/srf_huds_Gm12878_rep1.FORWARD.RIGHT.csv",
        ("rep2", "forward", "center"):
        "tests/test_data/srf_huds_Gm12878_rep2.FORWARD.CENTER.csv",
        ("rep2", "forward", "left"):
        "tests/test_data/srf_huds_Gm12878_rep2.FORWARD.LEFT.csv",
        ("rep2", "forward", "right"):
        "tests/test_data/srf_huds_Gm12878_rep2.FORWARD.RIGHT.csv",
        ("rep1", "reverse", "center"):
        "tests/test_data/srf_huds_Gm12878_rep1.REVERSE.CENTER.csv",
        ("rep1", "reverse", "left"):
        "tests/test_data/srf_huds_Gm12878_rep1.REVERSE.LEFT.csv",
        ("rep1", "reverse", "right"):
        "tests/test_data/srf_huds_Gm12878_rep1.REVERSE.RIGHT.csv",
        ("rep2", "reverse", "center"):
        "tests/test_data/srf_huds_Gm12878_rep2.REVERSE.CENTER.csv",
        ("rep2", "reverse", "left"):
        "tests/test_data/srf_huds_Gm12878_rep2.REVERSE.LEFT.csv",
        ("rep2", "reverse", "right"):
        "tests/test_data/srf_huds_Gm12878_rep2.REVERSE.RIGHT.csv"
    }
