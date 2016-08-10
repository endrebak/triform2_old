import pandas as pd

from joblib import Parallel, delayed
from rpy2.robjects import r, pandas2ri, globalenv
pandas2ri.activate()
ri2py = pandas2ri.ri2py

from bx.intervals.intersection import IntervalTree


def exclude_redundant_peaks(indata, args):

    return Parallel(n_jobs=args.number_cores)(
        delayed(_exclude_redundant_peaks)(data, args)
        for chromosome, data in indata.items())


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


def _exclude_type23_overlapping1(type1_it, type2_locs, type3_locs):

    keys = []
    for locs in [type2_locs, type3_locs]:
        for k, (start, end) in locs.iterrows():

            intervals = type1_it.find(start, end)
            if intervals:
                continue

            keys.append(k)

    return keys


def _exclude_redundant_peaks(indata, args):

    type1_peaks = ri2py(indata["peak_info"][1])
    type1_locs = _locs_from_df(type1_peaks)
    type1_it = _create_intervaltree(type1_locs)

    type1_keys = type1_it.find(0, type1_locs.End.max())

    type2_peaks = ri2py(indata["peak_info"][2])
    type3_peaks = ri2py(indata["peak_info"][3])
    type2_locs = _locs_from_df(type2_peaks)
    type3_locs = _locs_from_df(type3_peaks)

    type23_keys = _exclude_type23_overlapping1(type1_it, type2_locs,
                                               type3_locs)
    keys = type1_keys + type23_keys

    indexes = pd.DataFrame(index=pd.Series(keys))

    df = pd.concat([type1_peaks, type2_peaks, type3_peaks], axis=0)
    df = df.ix[indexes.index]
    df = pandas2ri.DataFrame(df)

    return df
