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


@pytest.fixture(scope="session")
def args():
    return MockNamespace(1, "hg19", 100, 150, 0.375, 0.1, 10, 10,
                         ["examples/srf_huds_Gm12878_rep1.bed",
                          "examples/srf_huds_Gm12878_rep2.bed"],
                         ["examples/backgr_huds_Gm12878_rep1.bed",
                          "examples/backgr_huds_Gm12878_rep2.bed"])


@pytest.fixture(scope="session")
def indata_make_chromosome_cover_files_chip():
    return {"rep1": "tests/test_data/srf_huds_Gm12878_rep1.csv",
            "rep2": "tests/test_data/srf_huds_Gm12878_rep2.csv"}


@pytest.fixture(scope="session")
def indata_make_chromosome_cover_files_input():
    return {"rep1": "tests/test_data/backgr_huds_Gm12878_rep1.csv",
            "rep2": "tests/test_data/backgr_huds_Gm12878_rep2.csv"}


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
        ("examples/backgr_huds_Gm12878_rep1.bed", "reverse"):
        "tests/test_data/chrY_backgr_huds_Gm12878_rep1_neg.csv",
        ("examples/backgr_huds_Gm12878_rep1.bed", "forward"):
        "tests/test_data/chrY_backgr_huds_Gm12878_rep1_pos.csv",
        ("examples/backgr_huds_Gm12878_rep2.bed", "reverse"):
        "tests/test_data/chrY_backgr_huds_Gm12878_rep2_neg.csv",
        ("examples/backgr_huds_Gm12878_rep2.bed", "forward"):
        "tests/test_data/chrY_backgr_huds_Gm12878_rep2_pos.csv",
    }


@pytest.fixture(scope="session")
def chromosome_cover_results_chip():
    return {
        ("examples/srf_huds_Gm12878_rep1.bed", "reverse"):
        "tests/test_data/chrY_srf_huds_Gm12878_rep1_neg.csv",
        ("examples/srf_huds_Gm12878_rep1.bed", "forward"):
        "tests/test_data/chrY_srf_huds_Gm12878_rep1_pos.csv",
        ("examples/srf_huds_Gm12878_rep2.bed", "reverse"):
        "tests/test_data/chrY_srf_huds_Gm12878_rep2_neg.csv",
        ("examples/srf_huds_Gm12878_rep2.bed", "forward"):
        "tests/test_data/chrY_srf_huds_Gm12878_rep2_pos.csv"
    }


@pytest.fixture(scope="session")
def init_result_input():
    return {
        ("examples/srf_huds_Gm12878_rep1.bed", "forward", "center"):
        "tests/test_data/backgr_huds_Gm12878_rep1.FORWARD.CENTER.csv",
        ("rep2", "forward", "center"):
        "tests/test_data/backgr_huds_Gm12878_rep2.FORWARD.CENTER.csv",
        ("examples/srf_huds_Gm12878_rep1.bed", "reverse", "center"):
        "tests/test_data/backgr_huds_Gm12878_rep1.REVERSE.CENTER.csv",
        ("rep2", "reverse", "center"):
        "tests/test_data/backgr_huds_Gm12878_rep2.REVERSE.CENTER.csv",
    }


@pytest.fixture(scope="session")
def init_result_chip():
    return {
        ("examples/srf_huds_Gm12878_rep1.bed", "forward", "center"):
        "tests/test_data/actual_result_rep1_forward_center.csv",
        ("examples/srf_huds_Gm12878_rep1.bed", "forward", "left"):
        "tests/test_data/actual_result_rep1_forward_left.csv",
        ("examples/srf_huds_Gm12878_rep1.bed", "forward", "right"):
        "tests/test_data/actual_result_rep1_forward_right.csv",
        ("examples/srf_huds_Gm12878_rep2.bed", "forward", "center"):
        "tests/test_data/actual_result_rep2_forward_center.csv",
        ("examples/srf_huds_Gm12878_rep2.bed", "forward", "left"):
        "tests/test_data/actual_result_rep2_forward_left.csv",
        ("examples/srf_huds_Gm12878_rep2.bed", "forward", "right"):
        "tests/test_data/actual_result_rep2_forward_right.csv",
        ("examples/srf_huds_Gm12878_rep1.bed", "reverse", "center"):
        "tests/test_data/actual_result_rep1_reverse_center.csv",
        ("examples/srf_huds_Gm12878_rep1.bed", "reverse", "left"):
        "tests/test_data/actual_result_rep1_reverse_left.csv",
        ("examples/srf_huds_Gm12878_rep1.bed", "reverse", "right"):
        "tests/test_data/actual_result_rep1_reverse_right.csv",
        ("examples/srf_huds_Gm12878_rep2.bed", "reverse", "center"):
        "tests/test_data/actual_result_rep2_reverse_center.csv",
        ("examples/srf_huds_Gm12878_rep2.bed", "reverse", "left"):
        "tests/test_data/actual_result_rep2_reverse_left.csv",
        ("examples/srf_huds_Gm12878_rep2.bed", "reverse", "right"):
        "tests/test_data/actual_result_rep2_reverse_right.csv"
    }
