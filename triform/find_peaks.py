from collections import defaultdict

from natsort import natsorted
from joblib import Parallel, delayed
from rpy2.robjects import r, pandas2ri
ri2py = pandas2ri.ri2py

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

    # print(chromosome)
    merged_peaks = {}
    merged_info = {}

    pos_cvg = forward["cvg"]
    neg_cvg = reverse["cvg"]

    number_peaks = 0
    for i in range(1, 4):
        # print(i)
        p1 = reverse["peaks"][i]
        p2 = forward["peaks"][i]

        # if(!length(p1) | !length(p2)) next
        either_empty = ri2py(r("function(p1, p2) (!length(p1) | !length(p2))")(
            p1, p2))[0]
        if either_empty:
            # print("either_empty")
            continue

        find_overlaps = r(
            "function(p1, p2) matrix(as.matrix(findOverlaps(p1,p2)),ncol=2)")
        ov = find_overlaps(p1, p2)

        # if(!nrow(ov)) next
        overlap = ri2py(r["nrow"](ov))[0]
        if not overlap:
            # print("not overlap")
            continue

        _find_duplicates = r(
            "function(ov, col) (ov[,col] %in% ov[duplicated(ov[,col]),col])")
        dup1 = _find_duplicates(ov, 1)
        dup2 = _find_duplicates(ov, 2)

        is_multi = r["|"](dup1, dup2)
        all_multi = ri2py(r["all"](is_multi))[0]
        if all_multi:
            # print("all_multi")
            continue

        rows = r["!"](is_multi)
        ov = subset_RS4_rows(ov, rows)

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

        #if(!any(ok)) next
        if not ri2py(r["any"](ok))[0]:
            # print('not ri2py(r["any"](ok))[0]')
            continue

        ov = subset_RS4_rows(ov, ok)
        peaks = subset_RS4(peaks, ok)
        n_peaks = r["length"](peaks)[0]

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
        assertions = []
        for var in [peak_nlps, max_nlps, peak_locs, peaks, peak_cvg, peak_surL,
                    peak_surR]:
            length = ri2py(r["length"](var))[0]
            assertions.append(length >= 1)

        if not all(assertions):
            # print("not all(assertions)")
            continue

        merged_peaks[i] = peaks

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
