from collections import defaultdict
from itertools import product

import pandas as pd
import numpy as np

import pytest

from rpy2.robjects import r, pandas2ri
ri2py = pandas2ri.ri2py
from rpy2.robjects.packages import importr
importr("GenomicRanges")
importr("S4Vectors")
from triform.init import _init
from triform.chromosome import chromosome
from triform.helper_functions import (
    subset_RS4, subset_RS4_rows, subset_RS4_cols, df_to_iranges, df_to_rle)


@pytest.fixture
def input_data():

    results = defaultdict(dict)
    for peak_type, direction in product([1, 2, 3], ["reverse", "forward"]):

        df = r["read.table"]("tests/test_data/chromosome_result_%s%s.csv" % (
            direction, peak_type),
                             sep=" ")
        df = r('function(df) df[c("PEAK.START", "PEAK.END", "PEAK.WIDTH")]')(
            df)
        iranges = df_to_iranges(df)

        results["chrY"][direction, peak_type] = iranges

    for direction in "reverse forward".split():
        df = r["read.table"]("tests/test_data/chromosome_cvg_%s.csv" %
                             direction)
        results["chrY"][direction, "cvg"] = df_to_rle(df)

    return results


@pytest.fixture
def expected_result():

    results = {}
    for peak_type in range(1, 4):

        df = r["read.table"]("tests/test_data/find_peaks_result_%s.csv" %
                             peak_type,
                             sep=" ")
        results[peak_type] = df

    return df


@pytest.mark.current
def test_find_peaks(input_data, expected_result, args):
    peaks = _find_peaks(input_data["chrY"], args)

    assert expected_result[1] == peaks[1]
    assert 0


def _find_peaks(indata, args):

    merged_peaks = {}
    merged_info = {}

    pos_cvg = indata["forward", "cvg"]
    neg_cvg = indata["reverse", "cvg"]

    for i in range(1, 4):
        p1 = indata["reverse", i]
        p2 = indata["forward", i]
        print(i)
        print("p1")
        print(p1)
        print("p2")
        print(p2)

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

        _find_peaks = r(
            """function(p1, p2) IRanges(start=pmin(start(p1),start(p2)), end=pmax(end(p1),end(p2)))""")
        peaks = _find_peaks(p1, p2)

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
        n_peaks = r["length"](peaks)

        merged_peaks[i] = peaks

        merged_info["reverse"] = subset_RS4_cols(ov, 1)
        merged_info["forward"] = subset_RS4_cols(ov, 2)

        # print(subset_RS4_cols(ov, 1), 'subset_RS4_cols(ov, 1)')
        # print(subset_RS4_cols(ov, 2), 'subset_RS4_cols(ov, 2)')
        info1 = subset_RS4_rows(indata["peak_info"]["reverse"],
                                subset_RS4_cols(ov, 1))
        info2 = subset_RS4_rows(indata["peak_info"]["forward"],
                                subset_RS4_cols(ov, 2))
        # print("info1 " * 100)
        # print(info1)
        # print("info2 " * 100)
        # print(info2)
        # raise
        #                         )
        # print(info1)

        #info1 <- PEAK.INFO[[type]][[1]][[i]][ov[,1],]
        #info2 <- PEAK.INFO[[type]][[2]][[i]][ov[,2],]
        #peak.locs <- as.integer(round((info1$PEAK.LOC + info2$PEAK.LOC)/2))

        #peak.cvg <- info1$PEAK.CVG + info2$PEAK.CVG
        #peak.surL <- info1$PEAK.SURL + info2$PEAK.SURL
        #peak.surR <- info1$PEAK.SURR + info2$PEAK.SURR

        #switch(i,
        #			 {zscores <- zscore(peak.cvg, peak.surL+peak.surR,2)
        #				max.z <- zscore(peak.cvg+peak.surL+peak.surR,0,2)},
        #			 {zscores <- zscore(peak.cvg, peak.surL)
        #				max.z <- zscore(peak.cvg+peak.surL,0)},
        #			 {zscores <- zscore(peak.cvg, peak.surR)
        #				max.z <- zscore(peak.cvg+peak.surR,0)})

        #peak.nlps <- -pnorm(zscores, low=FALSE, log=TRUE)/log(10)
        #max.nlps <- -pnorm(max.z, low=FALSE, log=TRUE)/log(10)

        #dfr <- data.frame(NLP=peak.nlps, MAX.NLP=max.nlps, LOC=peak.locs,
        #									WIDTH=width(peaks), START=start(peaks), END=end(peaks),
        #									CVG=peak.cvg, SURL=peak.surL, SURR=peak.surR, FORM=i)

        #rownames(dfr) <- with(dfr,sprintf("%s:%d-%d:%d",CHR,START,END,FORM))
        #PEAK.INFO[[type]][[direction]][[i]] <<- dfr
        #N.PEAKS <<- N.PEAKS + n.peaks
        #
        # }
