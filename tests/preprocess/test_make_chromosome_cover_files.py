import pytest
from numpy import allclose, int32
import pandas as pd
from io import StringIO

from rpy2 import robjects as ro
from rpy2.robjects import r, pandas2ri
ri2py = pandas2ri.ri2py
from rpy2.robjects.packages import importr
importr("GenomicRanges")

from triform.preprocess.make_chromosome_cover_files import _make_chromosome_cover_files


@pytest.fixture
def input_data():
    """Creates a GRanges object from the df below."""
    df = """seqnames start end strand
chrY 10641380 10641405 +
chrY 12314848 12314873 +
chrY 11743048 11743073 -
chrY 57396702 57396727 +
chrY 11778685 11778710 +
chrY 10595070 10595095 +
chrY 16891566 16891591 +
chrY 11936613 11936638 -
chrY 11924977 11925002 +
chrY 11764681 11764706 -"""

    command = """
lines = "{}"
options(stringsAsFactors=FALSE)
con <- textConnection(lines)
df <- read.csv(con, colClasses=c("character", rep("integer",2), "character"), header=1, sep=" ")
close(con)
colnames(df) <- c("seqnames", "start", "end", "strand")
    """.format(df)
    r(command)
    granges = r("makeGRangesFromDataFrame(df)")

    return granges


@pytest.fixture
def expected_result(input_data):

    create_expected_result = r("""function() {
  minus = "start end width
11743011 11743110 100
11936576 11936675 100
11764644 11764743 100"

  plus = "start end width
10641343 10641442 100
12314811 12314910 100
57396665 57396764 100
11778648 11778747 100
10595033 10595132 100
16891529 16891628 100
11924940 11925039 100"

    pcon <- textConnection(plus)
    mcon <- textConnection(minus)

    pdf = read.delim(pcon, sep=" ", header=1)
    mdf = read.delim(mcon, sep=" ", header=1)

    pi = IRanges(start=pdf$start, end=pdf$end, width=pdf$width)
    mi = IRanges(start=mdf$start, end=mdf$end, width=mdf$width)

    covers_expected_result = coverage(IRangesList("-"=mi, "+"=pi))
    lsize = list("+"=7, "-"=3)
    covers = list(SIZE=lsize, CVG=covers_expected_result)
}
""")
    return create_expected_result()


@pytest.mark.unit
def test_make_chromosome_cover_files(input_data, expected_result, args_tests):
    result = _make_chromosome_cover_files(input_data, args_tests)

    print("result")
    print(result)
    print("expected_result")
    print(expected_result)

    all_equals = r["isTRUE"](r["all.equal"](expected_result, result))
    assert all_equals
