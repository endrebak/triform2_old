from bx.intervals.intersection import IntervalTree
from numpy import int64

from rpy2.robjects import r, pandas2ri
ri2py = pandas2ri.ri2py


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


def df_to_iranges(df):

    return r("""function(df) {
    IRanges(start=as.integer(df[[1]]), end=as.integer(df[[2]]), width=as.integer(df[[3]]))}
    """)(df)


def _locs_from_df(df):

    locs = df.index.get_level_values(0).to_series()
    locs = locs.str.split(":", expand=True).iloc[:, 1]
    locs = locs.str.split("-", expand=True).astype(int)
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
