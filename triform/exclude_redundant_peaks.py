import pandas as pd

from joblib import Parallel, delayed
from rpy2.robjects import r, pandas2ri, globalenv
pandas2ri.activate()
ri2py = pandas2ri.ri2py

from triform.helper_functions import _locs_from_df, _create_intervaltree


def exclude_redundant_peaks(indata, args):

    return Parallel(n_jobs=args.number_cores)(
        delayed(_exclude_redundant_peaks)(data, args)
        for chromosome, data in indata.items())


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

    if not 1 in indata["peak_info"]:
        return []:

    peaks1 = indata["peak_info"][1]
    type1_peaks = ri2py(peaks1)
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
