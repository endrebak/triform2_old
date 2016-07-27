from collections import defaultdict

from natsort import natsorted
from joblib import Parallel, delayed
from rpy2.robjects import r

from triform.helper_functions import (
    subset_RS4, subset_RS4_rows, subset_RS4_cols, df_to_iranges, df_to_rle)


def find_peaks(indata, args):

    chromosomes = natsorted(set([k[0] for k in indata.keys()]))

    results = Parallel(n_jobs=args.number_cores)(delayed(_find_peaks)(
        indata[chromosome, "forward"], indata[chromosome, "reverse"],
        chromosome, args) for chromosome in chromosomes)

    d = defaultdict(lambda: defaultdict(dict))
    for chromosome, (peaks, peak_info) in zip(chromosomes, results):

        for k, v in peaks.items():
            d[chromosome]["peaks"][k] = v

        for k, v in peak_info.items():
            d[chromosome]["peak_info"][k] = v

    return d


r('''
zscore <- function(x,y,r=1) {  # r = size.y/size.x
  dif <- (r*x-y)
  zs <- dif/sqrt(r*(x+y))
  zs[!dif] <- 0
  zs
}
''')


def _find_peaks(forward, reverse, chromosome, args):

    merged_peaks = {}
    merged_info = {}

    pos_cvg = forward["cvg"]
    neg_cvg = reverse["cvg"]

    # print("pos_cvg")
    # print(pos_cvg)
    # print("neg_cvg")
    # print(neg_cvg)

    # pos_cvg
    # integer-Rle of length 57442693 with 22249 runs
    #   Lengths: 2709647      29      71      29 ...     100      20     100      75
    #   Values :       0       1       2       1 ...       1       0       1       0

    # neg_cvg
    # integer-Rle of length 57442693 with 21667 runs
    #   Lengths: 2709676     100   31062     100 ...       5      89     100       3
    #   Values :       0       1       0       1 ...       1       0       1       0

    number_peaks = 0
    for i in range(1, 4):
        p1 = reverse["peaks"][i]
        p2 = forward["peaks"][i]
        # print("p1")
        # print(p1)
        # print("p2")
        # print(p2)

        # p1
        # IRanges of length 138
        #          start      end width
        # [1]    3778332  3778354    23
        # [2]    3779776  3779843    68
        # [3]    7589496  7589586    91
        # [4]   10534276 10534287    12
        # [5]   10568307 10568406   100
        # ...        ...      ...   ...
        # [134] 57415256 57415288    33
        # [135] 57420930 57420942    13
        # [136] 57427449 57427541    93
        # [137] 57428262 57428328    67
        # [138] 57442124 57442137    14

        # p2
        # IRanges of length 150
        #          start      end width
        # [1]    2928874  2928970    97
        # [2]    3772923  3773021    99
        # [3]    4012212  4012307    96
        # [4]    7087921  7088026   106
        # [5]    7184247  7184346   100
        # ...        ...      ...   ...
        # [146] 57414596 57414627    32
        # [147] 57414655 57414665    11
        # [148] 57420858 57420934    77
        # [149] 57434349 57434371    23
        # [150] 57435799 57435896    98

        # TODO:
        # if(!length(p1) | !length(p2)) next

        find_overlaps = r(
            "function(p1, p2) matrix(as.matrix(findOverlaps(p1,p2)),ncol=2)")
        ov = find_overlaps(p1, p2)
        # TODO:
        # if(!nrow(ov)) next

        _find_duplicates = r(
            "function(ov, col) (ov[,col] %in% ov[duplicated(ov[,col]),col])")
        dup1 = _find_duplicates(ov, 1)
        dup2 = _find_duplicates(ov, 2)
        is_multi = r["!"](r["|"](dup1, dup2))
        # TODO:
        ov = subset_RS4_rows(ov, is_multi)

        subset_peaks = r("function(p, ov, col) p[1:length(p) %in% ov[,col]]")
        p1 = subset_peaks(p1, ov, 1)
        p2 = subset_peaks(p2, ov, 2)

        __find_peaks = r(
            """function(p1, p2) IRanges(start=pmin(start(p1),start(p2)), end=pmax(end(p1),end(p2)))""")
        peaks = __find_peaks(p1, p2)

        _shift_peaks = r("""function (peaks, i, FLANK.DELTA) {
          ranges = switch(i,
            ranges <- IRanges(start=start(peaks)-FLANK.DELTA, end=end(peaks)+FLANK.DELTA),
            ranges <- IRanges(start=start(peaks)-FLANK.DELTA, end=end(peaks)),
            ranges <- IRanges(start=start(peaks), end=end(peaks)+FLANK.DELTA))
        }""")

        ranges = _shift_peaks(peaks, i, args.flank_distance)

        _viewApply = r(
            "function(cvg, ranges) viewApply(Views(cvg, ranges), as.numeric)")

        neg_peak_cvg = _viewApply(neg_cvg, ranges)
        pos_peak_cvg = _viewApply(pos_cvg, ranges)

        lags = r("""function (x, y) mapply(function(x,y) {
        							 cc=ccf(x,y,lag.max=100,plot=FALSE)
        							 with(cc,lag[which.max(acf)])}, x=x, y=y)

        """)(x=neg_peak_cvg,
             y=pos_peak_cvg)

        ok = r[">"](lags, args.min_shift)

        # TODO:
        #if(!any(ok)) next

        ov = subset_RS4_rows(ov, ok)
        peaks = subset_RS4(peaks, ok)
        n_peaks = r["length"](peaks)[0]

        merged_peaks[i] = peaks

        # not really used anywhere?
        # merged_info[i, 1] = subset_RS4_cols(ov, 1)
        # merged_info[i, 2] = subset_RS4_cols(ov, 2)

        # print('reverse["peak_info"][i]')
        # print(r["head"](reverse["peak_info"][i]))

        # print('forward["peak_info"][i]')
        # print(r["head"](forward["peak_info"][i]))

        # indata["reverse", i, "peak_info"]
        #   PEAK.LOC PEAK.CVG PEAK.SURL PEAK.SURR PEAK.NLP PEAK.WIDTH PEAK.START PEAK.END
        # 1  3778343        3         0         0    2.146         23    3778332  3778354
        # 2  3779810        4         0         0    2.631         68    3779776  3779843
        # 3  7589541        4         0         0    2.631         91    7589496  7589586
        # 4 10534282        3         0         0    2.146         12   10534276 10534287
        # 5 10568356        3         0         0    2.146        100   10568307 10568406
        # 6 10591945        4         0         0    2.631         73   10591909 10591981

        # indata["forward", i, "peak_info"]
        #   PEAK.LOC PEAK.CVG PEAK.SURL PEAK.SURR PEAK.NLP PEAK.WIDTH PEAK.START PEAK.END
        # 1  2928922        5         0         0    3.106         97    2928874  2928970
        # 2  3772972        3         0         0    2.631         99    3772923  3773021
        # 3  4012260        3         0         0    2.146         96    4012212  4012307
        # 4  7087974       14         0         0    7.217        106    7087921  7088026
        # 5  7184296        7         0         0    4.039        100    7184247  7184346
        # 6  8361022        5         0         0    3.106        102    8360972  8361073

        info1 = subset_RS4_rows(reverse["peak_info"][i],
                                subset_RS4_cols(ov, 1))
        info2 = subset_RS4_rows(forward["peak_info"][i],
                                subset_RS4_cols(ov, 2))

        peak_locs = r(
            "function(info1, info2) as.integer(round((info1$PEAK.LOC + info2$PEAK.LOC)/2))")(
                info1, info2)
        peak_cvg = r("function(info1, info2) info1$PEAK.CVG + info2$PEAK.CVG")(
            info1, info2)
        peak_surL = r(
            "function(info1, info2) info1$PEAK.SURL + info2$PEAK.SURL")(info1,
                                                                        info2)
        peak_surR = r(
            "function(info1, info2) info1$PEAK.SURR + info2$PEAK.SURR")(info1,
                                                                        info2)

        _compute_z = r('''function(i, peak.cvg, peak.surL, peak.surR) {
        switch(i,
               {zscores <- zscore(peak.cvg, peak.surL+peak.surR,2)
                max.z <- zscore(peak.cvg+peak.surL+peak.surR,0,2)},
               {zscores <- zscore(peak.cvg, peak.surL)
                max.z <- zscore(peak.cvg+peak.surL,0)},
               {zscores <- zscore(peak.cvg, peak.surR)
                max.z <- zscore(peak.cvg+peak.surR,0)})
        list(zscores, max.z)
        }
        ''')

        l = _compute_z(i, peak_cvg, peak_surL, peak_surR)
        zscores, maxz = l

        peak_nlps = r(
            "function(zscores) -pnorm(zscores, low=FALSE, log=TRUE)/log(10)")(
                zscores)
        max_nlps = r(
            "function(max.z) -pnorm(max.z, low=FALSE, log=TRUE)/log(10)")(maxz)

        # for x in "peak_nlps, max_nlps, peak_locs, peaks, peak_cvg, peak_surL, peak_surR, i, chromosome".split(
        #         ", "):
        #     var = vars()[x]
        #     print(x)
        #     print(r["length"](var))

        dfr = r(
            """function(peak.nlps, max.nlps, peak.locs, peaks, peak.cvg, peak.surL, peak.surR, i, CHR) {

            dfr = data.frame(NLP=peak.nlps, MAX.NLP=max.nlps, LOC=peak.locs,
                                                WIDTH=width(peaks), START=start(peaks), END=end(peaks),
                                                CVG=peak.cvg, SURL=peak.surL, SURR=peak.surR, FORM=i)
            rownames(dfr) <- with(dfr,sprintf("%s:%d-%d:%d",CHR,START,END,FORM))
            dfr
            }""")(peak_nlps, max_nlps, peak_locs, peaks, peak_cvg, peak_surL,
                  peak_surR, i, chromosome)
        number_peaks += n_peaks
        merged_info[i] = dfr

    return merged_peaks, merged_info
