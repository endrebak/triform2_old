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
