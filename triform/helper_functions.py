from quicksect import IntervalTree
from numpy import int64

from rpy2.robjects import r, pandas2ri
ri2py = pandas2ri.ri2py

from rpy2.robjects.packages import importr
# importr("S4Vectors")
importr("GenomicRanges")


def subset_RS4(rs4, subset, drop=False):
    subset_func = r("""function(o, s, d){
    o[s, drop=d]
    }
    """)
    return subset_func(rs4, subset, drop)


def subset_RS4_2(rs4, subset, drop=False):
    subset_func = r("""function(o, s, d){
    o[[s]]
    }
    """)
    return subset_func(rs4, subset, drop)


def subset_RS4_rows(rs4, subset, drop=False):
    subset_func = r("""function(o, s, d){
    o[s, ,drop=d]
    }
    """)
    return subset_func(rs4, subset, drop)


def subset_RS4_cols(rs4, subset, drop=False):
    subset_func = r("""function(o, s, d){
    o[,s,drop=d]
    }
    """)
    return subset_func(rs4, subset, drop)


def df_to_rle(df):

    _df_to_rle = r("""function(df){
    Rle(df$values, df$lengths)
    }
    """)

    return _df_to_rle(df)


def rle_to_df(rle):

    _rle_to_df = r("""function(rle){
    df = data.frame(as.vector(runLength(rle)), as.vector(runValue(rle)))
    colnames(df) = c("lengths", "values")
    df
    }
    """)

    return _rle_to_df(rle)


def rle_to_pandas_df(rle):

    return ri2py(rle_to_df(rle)).astype(int64)


def df_to_iranges_end(df):

    return r("""function(df) {
    IRanges(start=as.integer(df[[1]]), end=as.integer(df[[2]]))}
    """)(df)


def df_to_iranges(df):

    return r("""function(df) {
    IRanges(start=as.integer(df[[1]]), end=as.integer(df[[2]]), width=as.integer(df[[3]]))}
    """)(df)


def _locs_from_df(df):

    locs = df.index.get_level_values(0).to_series()
    locs = locs.astype(str).str.split(":", expand=True).iloc[:, 1]
    locs = locs.astype(str).str.split("-", expand=True).astype(int)
    locs.columns = "Start End".split()

    return locs


def _create_intervaltree(locs):

    it = IntervalTree()

    for k, (start, end) in locs.iterrows():

        intervals = it.find(start, end)
        if intervals:
            continue

        it.add(start, end, k)

    return it


def granges_to_bed_df(gr):

    _granges_to_df = r("""function(gr)
    {df <- data.frame(seqnames=seqnames(gr),
    starts=start(gr),
    ends=end(gr),
    names=c(rep(".", length(gr))),
    scores=elementMetadata(gr)[,1],
    strands=strand(gr))
    df}""")

    return _granges_to_df(gr)


def bed_df_to_granges(df):

    _bed_df_to_granges = r("""function(df){
    GRanges(seqnames=df$seqnames, ranges=IRanges(df$starts, df$ends), strand=df$strands, counts=df$scores)
    }
    """)

    return _bed_df_to_granges(df)
