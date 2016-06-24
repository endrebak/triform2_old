import pytest
from numpy import allclose, int32
import pandas as pd
from io import StringIO

from rpy2.robjects import r, pandas2ri
ri2py = pandas2ri.ri2py
from rpy2.robjects.packages import importr


@pytest.fixture
def input_data():

    return 'chrY\t10641380\t10641405\t+\nchrY\t12314848\t12314873\t+\nchrY\t11743048\t11743073\t-\nchrY\t57396702\t57396727\t+\nchrY\t11778685\t11778710\t+\nchrY\t10595070\t10595095\t+\nchrY\t16891566\t16891591\t+\nchrY\t11936613\t11936638\t-\nchrY\t11924977\t11925002\t+\nchrY\t11764681\t11764706\t-\n'


@pytest.fixture
def expected_result():
    input = StringIO(u"""chrY  10641380  10641405  +
chrY  12314848  12314873  +
chrY  11743048  11743073  -
chrY  57396702  57396727  +
chrY  11778685  11778710  +
chrY  10595070  10595095  +
chrY  16891566  16891591  +
chrY  11936613  11936638  -
chrY  11924977  11925002  +
chrY  11764681  11764706  -""")
    return pd.read_table(input,
                         sep="\s+",
                         header=None,
                         dtype={0: str,
                                1: int32,
                                2: int32})


def make_ranged_data(lines):
    """Write description as command ending in a period."""

    importr("GenomicRanges")

    _make_ranged_data = r("""
makeRangedData <- function(lines){
  options(stringsAsFactors=FALSE)
  con <- textConnection(lines)
  dfr <- read.delim(con, colClasses=c("character", rep("integer",2), "character"), header=FALSE)
  close(con)
  # dfr <- read.delim(infile, header=FALSE,
  #                   colClasses=c("character", rep("integer",2), "character"))
  colnames(dfr) <- c("seqnames", "start", "end", "strand")
  rd = makeGRangesFromDataFrame(dfr)
  # rd <- as(dfr, "RangedData")
  # save(rd, file=outfile)
}
""")
    return _make_ranged_data(lines)


@pytest.mark.unit
def test_name(input_data, expected_result):
    result = make_ranged_data(input_data)

    bc = importr("BiocGenerics")
    gidb = importr("GenomeInfoDb")
    base = importr("base")

    seqnames = base.as_vector(gidb.seqnames(result))
    starts = bc.start(result)
    ends = bc.end(result)
    strands = base.as_vector(bc.strand(result))

    bed_cols = [pd.Series(ri2py(c)) for c in [seqnames, starts, ends, strands]]
    bed_dataframe = pd.concat(bed_cols, axis=1)
    print(bed_dataframe)
    print(expected_result)
    print(bed_dataframe.dtypes)
    print(expected_result.dtypes)
    assert bed_dataframe.equals(expected_result)
