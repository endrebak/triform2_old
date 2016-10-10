import logging
import pandas as pd

from joblib import Parallel, delayed
from rpy2.robjects import r, pandas2ri, globalenv
pandas2ri.activate()
ri2py = pandas2ri.ri2py

from triform.helper_functions import _locs_from_df, _create_intervaltree



add_rownames = r("""function(df, rnames) {
rownames(df) = rnames
df
}""")


def exclude_redundant_peaks(indata, args):


    ordered_data = list(indata.items())
    chromosomes = [c for c, _ in ordered_data]

    dfs = Parallel(n_jobs=args.number_cores)(
        delayed(_exclude_redundant_peaks)(data, args)
        for chromosome, data in ordered_data)

    return dfs


def _exclude_type23_overlapping1(type1_it, type2_locs, type3_locs):

    keys = []
    for locs in [type2_locs, type3_locs]:

        if locs.empty:
            continue

        for k, (start, end) in locs.iterrows():

            intervals = type1_it.find(start, end)
            if intervals:
                continue

            keys.append(k)

    return keys


def _exclude_redundant_peaks(indata, args):

    peak_types = set()
    for i in range(1, 4):
        if i in indata["peak_info"]:
            peak_types.add(i)
            # print(indata["peak_info"][i])

    # print("peak_types")
    # print(peak_types)
    # peak_types = set(indata["peak_info"].keys())
    if peak_types == set([1, 2, 3]):
        # logging.info("All three exist.")
        return _exclude_redundant_peaks_all_three_exist(indata, args)

    if peak_types == set([1, 2]):
        # logging.info("One and two exist.")
        peaks1 = ri2py(indata["peak_info"][1])
        peaks2 = ri2py(indata["peak_info"][2])
        return _exclude_redundant_peaks_one_and_other_exist(peaks1, peaks2)

    if peak_types == set([1, 3]):
        # logging.info("One and three exist.")
        peaks1 = ri2py(indata["peak_info"][1])
        peaks3 = ri2py(indata["peak_info"][3])
        return _exclude_redundant_peaks_one_and_other_exist(peaks1, peaks3)

    if len(peak_types) == 1:
        # logging.info("Only one exists! " * 3)
        return list(indata["peak_info"].values())[0]

    if len(peak_types) == 0:
        # logging.info("No peaks exist! " * 3)
        return ri2py(pd.DataFrame(
            columns=
            "NLP  MAX.NLP      LOC WIDTH    START      END CVG  SURL SURR FORM".split(
            )))


def _exclude_redundant_peaks_all_three_exist(indata, args):

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
    length = len(df)

    df = pandas2ri.DataFrame(df)

    # ryp2 loses the rownames if there is only one row :(
    if length == 1:
        df = add_rownames(df, list(keys))

    return df


def _exclude_redundant_peaks_one_and_other_exist(peaks1, peaks_other):

    type1_peaks = peaks1
    type1_locs = _locs_from_df(type1_peaks)
    type1_it = _create_intervaltree(type1_locs)

    type1_keys = type1_it.find(0, type1_locs.End.max())

    type_other_locs = _locs_from_df(peaks_other)

    type23_keys = _exclude_type23_overlapping1(type1_it, type_other_locs,
                                               pd.DataFrame())

    keys = type1_keys + type23_keys

    indexes = pd.DataFrame(index=pd.Series(keys))

    df = pd.concat([type1_peaks, peaks_other], axis=0)
    df = df.ix[indexes.index]

    length = len(df)

    df = pandas2ri.DataFrame(df)


    # ryp2 loses the rownames if there is only one row :(
    if length == 1:
        df = add_rownames(df, list(keys))

    return df
