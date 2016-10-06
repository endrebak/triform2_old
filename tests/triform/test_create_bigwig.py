from collections import defaultdict

try:
    from itertools import izip_longest
except ImportError:
    from itertools import zip_longest

from io import StringIO
from os import stat
from os.path import dirname

import pandas as pd
import numpy as np

import pytest

from rpy2.robjects import r, pandas2ri
from triform.create_bigwig import _create_bigwig


@pytest.fixture
def input():
    bed_file = "tests/test_data/iranges_to_bigwig.bed"
    df = r["read.table"](bed_file, sep=" ")

    treatment = {bed_file: {
        "chrY":
        r("function(df) {IRanges(start=df$start, width=df$width)}")(df),
        "chr1": r("function(df) {IRanges(start=df$start, width=df$width)}")(df)
    }}

    return treatment


@pytest.fixture
def output_bigwig(tmpdir):
    p = tmpdir.mkdir("sub").join("outfile.bw")
    return str(p)


@pytest.mark.unit
def test_create_bigwig(input, output_bigwig, args):
    print(input)
    print(output_bigwig)

    _create_bigwig(list(input.values())[0], output_bigwig, dirname(output_bigwig), args)
    filesize = stat(output_bigwig).st_size

    print(filesize, "filesize")

    assert filesize > 0
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
#     create_bigwig(treatment_iranges, args.bigwig, args)
#   File "/local/home/endrebak/anaconda3/lib/python3.5/site-packages/triform2-0.0.5-py3.5.egg/triform/create_bigwig.py", line 20, in create_bigwig
#     _create_bigwig(iranges, outpath, args)
#   File "/local/home/endrebak/anaconda3/lib/python3.5/site-packages/triform2-0.0.5-py3.5.egg/triform/create_bigwig.py", line 37, in _create_bigwig
#     genome = r["RleList"](genome_coverage_rpkm)
#   File "/local/home/endrebak/anaconda3/lib/python3.5/site-packages/rpy2/robjects/functions.py", line 178, in __call__
#     return super(SignatureTranslatedFunction, self).__call__(*args, **kwargs)
#   File "/local/home/endrebak/anaconda3/lib/python3.5/site-packages/rpy2/robjects/functions.py", line 106, in __call__
#     res = super(Function, self).__call__(*new_args, **new_kwargs)
# rpy2.rinterface.RRuntimeError: Error in .Call2("Rle_constructor", values, lengths, check, 0L, PACKAGE = "S4Vectors") :
#   integer overflow while summing elements in 'lengths'
