from rpy2.robjects import r


def subset_RS4(rs4, subset, drop=False):
    subset_func = r("""function(o, s, d){
    o[s, drop=d]
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
