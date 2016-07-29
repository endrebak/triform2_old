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


def _exclude_redundant_peaks(indata, args):

    it1 = IntervalTree()
    it23 = IntervalTree()
    # print("init")

    dfs = []
    for peak_type, v in indata["peak_info"].items():
        # print(peak_type)
        df = ri2py(v)

        locs = df.index.get_level_values(0).to_series()
        locs = locs.str.split(":", expand=True).iloc[:, 1]
        locs = locs.str.split("-", expand=True).astype(int)

        if peak_type == 1:
            keys = []
            for k, (start, end) in locs.iterrows():

                intervals = it1.find(start, end)
                if intervals:
                    continue

                it1.add(start, end, k)
                keys.append(k)

            indexes = pd.DataFrame(index=pd.Series(keys))
            found_rows = df.ix[indexes.index]
            dfs.append(found_rows)

        else:

            keys = []
            for k, (start, end) in locs.iterrows():

                intervals1 = it1.find(start, end)
                intervals23 = it23.find(start - args.read_width,
                                        end + args.read_width)
                if intervals1 or intervals23:
                    continue

                it23.add(start, end, k)
                keys.append(k)

            indexes = pd.DataFrame(index=pd.Series(keys))
            found_rows = df.ix[indexes.index]
            dfs.append(found_rows)

    df = pd.concat(dfs, axis=0)

    df = pandas2ri.DataFrame(df)
    # cols = "LOC WIDTH START END CVG SURL SURR FORM".split()

    globalenv["df"] = df
    df = r["df"]

    # print(type(df))
    # print(df)

    # df[cols] = df[cols].astype(int)

    return df

    # # peak_types = indata["peaks"].keys()
    # # nb_peak_types = len(peak_types)

    # # if nb_peak_types == 1:
    # #     return indata["peak_info"][peak_types[0]]
    # # if nb_peak_types == 2:
    # #     pass

    # p1 = indata["peaks"][1]
    # p2 = indata["peaks"][2]
    # p3 = indata["peaks"][3]

    # ov12 = r("function(p1, p2) matrix(as.matrix(findOverlaps(p1,p2)),ncol=2)")(
    #     p1, p2)

    # if r["nrow"](ov12):
    #     ex2 = r("function(p2, ov12) 1:length(p2) %in% ov12[,2]")(p2, ov12)
    #     inv_ex2 = r["!"](ex2)
    #     print("oo")
    #     print(r["length"](p2))
    #     print(r["length"](inv_ex2))
    #     p2 = subset_RS4_rows(p2, inv_ex2)
    #     info2 = subset_RS4_rows(indata["peak_info"][2], inv_ex2)

    # ov13 = r("function(p1, p3) matrix(as.matrix(findOverlaps(p1,p3)),ncol=2)")(
    #     p1, p3)
    # if r["nrow"](ov13):
    #     ex3 = r("function(p3, ov13) 1:length(p3) %in% ov13[,2]")(p3, ov13)
    #     inv_ex3 = r["!"](ex3)

    #     print(r["length"](p3))
    #     print(r["length"](inv_ex3))
    #     p3 = subset_RS4_rows(p3, inv_ex3)
    #     info3 = subset_RS4_rows(indata["peak_info"][3], inv_ex3)

    # # TODO: do not use harcoded variable
    # ov23 = r(
    #     "function(p2, p3) matrix(as.matrix(findOverlaps(p2,p3, maxgap=100)),ncol=2)")(
    #         p2, p3)
    # print("ov23")
    # print(ov23)

    # if r["nrow"](ov23):
    #     ex2 = r("function(p2, ov23) 1:length(p2) %in% ov23[,1]")(p2, ov23)
    #     inv_ex2 = r["!"](ex2)
    #     p2 = subset_RS4_rows(p2, inv_ex2)
    #     info2 = subset_RS4_rows(indata["peak_info"][2], inv_ex2)

    #     ex3 = r("function(p3, ov23) 1:length(p3) %in% ov23[,2]")(p3, ov23)
    #     inv_ex3 = r["!"](ex3)
    #     print(r["length"](p3))
    #     print(r["length"](inv_ex3))
    #     p3 = subset_RS4_rows(p3, inv_ex3)
    #     info3 = subset_RS4_rows(indata["peak_info"][3], inv_ex3)

    # peak_info = r["rbind"](indata["peak_info"][1], info2, info3)
    # # nubmer_peaks = r["nrow"](peak_info)[0]

    # return peak_info
