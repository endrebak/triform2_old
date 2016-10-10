from collections import defaultdict

try:
    from itertools import izip_longest
except ImportError:
    from itertools import zip_longest

from io import StringIO
from os import stat

import pandas as pd
import numpy as np

import pytest

from rpy2.robjects import r, pandas2ri
from triform.create_bed import create_bed


@pytest.fixture
def chip():
    bed_file = "tests/test_data/iranges_to_bigwig.bed"
    df = r["read.table"](bed_file, sep=" ")

    treatment = {bed_file: {
        "chrY":
        r("function(df) {IRanges(start=df$start, width=df$width)}")(df),
        "chr1": r("function(df) {IRanges(start=df$start, width=df$width)}")(df)
    }}

    return treatment


@pytest.fixture
def input():
    bed_file = "tests/test_data/input_iranges_to_bigwig.bed"
    df = r["read.table"](bed_file, sep=" ", header=1)

    control = {bed_file: {
        "chrY":
        r("function(df) {IRanges(start=df$start, width=df$width)}")(df),
        "chr1": r("function(df) {IRanges(start=df$start, width=df$width)}")(df)
    }}

    return control


@pytest.fixture
def enriched_regions():

    result = u"""CHROMOSOME START END QVAL NLQ NLP MAX.NLP LOC WIDTH CVG SURL SURR FORM
chrY 27096723 27096953 0.0030190321570660633 2.520132261086883 3.984008558280453 5.638949274037698 27096861 231 19.0 2.0 24.0 2
chr1 2709562 2719659 0.12753262199776602 0.8943787113752741 1.9970405839422427 4.269576380346016 2709612 98 12.0 3.0 9.0 2"""

    return pd.read_table(StringIO(result), index_col=[0, 1, 2], sep=" ", header=0)


@pytest.fixture
def expected_result():
    result = u"""Chromosome      Start        End      Qval  logFC      Loc
chr1  2709562.0  2719659.0  0.127533    0.0  2709612"""

    return pd.read_table(StringIO(result), header=0, index_col=[0, 1, 2], sep="\s+")



@pytest.mark.current
def test_create_bed(chip, input, enriched_regions, expected_result, args):
    print(chip, "chip " * 3)
    print(input, "input " * 3)
    print(expected_result)

    result = create_bed(chip, input, enriched_regions, args)
    print(result, "result")
    print(expected_result, "expected")

    assert np.allclose(result, expected_result)
# /local/home/endrebak/anaconda3/lib/python3.5/site-packages/rpy2/robjects/functions.py:106: UserWarning: Error in .Call2("Rle_constructor", values, lengths, check, 0L, PACKAGE = "S4Vectors") :
#   integer overflow while summing elements in 'lengths'

#   res = super(Function, self).__call__(*new_args, **new_kwargs)
# /local/home/endrebak/anaconda3/lib/python3.5/site-packages/rpy2/robjects/functions.py:106: UserWarning: In addition:
#   res = super(Function, self).__call__(*new_args, **new_kwargs)
# /local/home/endrebak/anaconda3/lib/python3.5/site-packages/rpy2/robjects/functions.py:106: UserWarning: There were 24 warnings (use warnings() to see them)
#   res = super(Function, self).__call__(*new_args, **new_kwargs)
# /local/home/endrebak/anaconda3/lib/python3.5/site-packages/rpy2/robjects/functions.py:106: UserWarning:

#   res = super(Function, self).__call__(*new_args, **new_kwargs)
# Traceback (most recent call last):
#   File "/local/home/endrebak/anaconda3/bin//triform2", line 4, in <module>
#     __import__('pkg_resources').run_script('triform2==0.0.5', 'triform2')
#   File "/local/home/endrebak/anaconda3/lib/python3.5/site-packages/setuptools-20.7.0-py3.5.egg/pkg_resources/__init__.py", line 719, in run_script
#   File "/local/home/endrebak/anaconda3/lib/python3.5/site-packages/setuptools-20.7.0-py3.5.egg/pkg_resources/__init__.py", line 1504, in run_script
#   File "/local/home/endrebak/anaconda3/lib/python3.5/site-packages/triform2-0.0.5-py3.5.egg/EGG-INFO/scripts/triform2", line 157, in <module>
#     run_triform(args)
#   File "/local/home/endrebak/anaconda3/lib/python3.5/site-packages/triform2-0.0.5-py3.5.egg/triform/run_triform.py", line 96, in run_triform
#     create_bed(treatment_iranges, args.bed, args)
#   File "/local/home/endrebak/anaconda3/lib/python3.5/site-packages/triform2-0.0.5-py3.5.egg/triform/create_bed.py", line 20, in create_bed
#     _create_bed(iranges, outpath, args)
#   File "/local/home/endrebak/anaconda3/lib/python3.5/site-packages/triform2-0.0.5-py3.5.egg/triform/create_bed.py", line 37, in _create_bed
#     genome = r["RleList"](genome_coverage_rpkm)
#   File "/local/home/endrebak/anaconda3/lib/python3.5/site-packages/rpy2/robjects/functions.py", line 178, in __call__
#     return super(SignatureTranslatedFunction, self).__call__(*args, **kwargs)
#   File "/local/home/endrebak/anaconda3/lib/python3.5/site-packages/rpy2/robjects/functions.py", line 106, in __call__
#     res = super(Function, self).__call__(*new_args, **new_kwargs)
# rpy2.rinterface.RRuntimeError: Error in .Call2("Rle_constructor", values, lengths, check, 0L, PACKAGE = "S4Vectors") :
#   integer overflow while summing elements in 'lengths'
