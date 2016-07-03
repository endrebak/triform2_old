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


@pytest.fixture
def input_data():

    results = defaultdict(dict)
    for peak_type, direction in product([1, 2, 3], ["reverse", "forward"]):
        df = r["read.table"]("tests/test_data/{}_peaks_{}.csv".format(
            peak_type, direction),
                             sep=" ")
        gr = r["IRanges"](start=df[0], end=df[1])
        results[direction][peak_type] = gr

    for direction in ["reverse", "forward"]:
        df = r["read.table"]("tests/test_data/cvg_{}.csv".format(direction),
                             sep=" ")
        results["cvg"][direction] = r["Rle"](df[0], df[1])

    return results


@pytest.mark.current
def test_find_peaks(input_data):
    find_peaks(input_data)

    assert 0


def find_peaks(result):
    print(result)
    print(result["cvg"]["forward"])

    pass
# 	direction <- "merged"
# 	PEAKS[[type]][[direction]] <<- list(IRanges(),IRanges(),IRanges())
# 	PEAK.INFO[[type]][[direction]] <<- list(NULL,NULL,NULL)
# 	PEAK.INFO[[type]][["regions"]] <<- list(NULL,NULL,NULL)

# 	neg.cvg <- CENTER.CVG[[type]][[1]]
# 	pos.cvg <- CENTER.CVG[[type]][[2]]

# 	for (i in 1:3) {
# 		p1 <- PEAKS[[type]][[1]][[i]]
# 		p2 <- PEAKS[[type]][[2]][[i]]
# 		if(!length(p1) | !length(p2)) next

# 		ov <- matrix(as.matrix(findOverlaps(p1,p2)),ncol=2)
# 		if(!nrow(ov)) next

# 		dup1 <- (ov[,1] %in% ov[duplicated(ov[,1]),1])
# 		dup2 <- (ov[,2] %in% ov[duplicated(ov[,2]),2])
# 		is.multi <- dup1 | dup2
# 		if(all(is.multi)) next
# 		ov <- ov[!is.multi,,drop=FALSE]

# 		p1 <- p1[1:length(p1) %in% ov[,1]]
# 		p2 <- p2[1:length(p2) %in% ov[,2]]
# 		peaks <- IRanges(start=pmin(start(p1),start(p2)),
# 										 end=pmax(end(p1),end(p2)))

# 		switch(i,
# 					 ranges <- IRanges(start=start(peaks)-FLANK.DELTA,
# 														 end=end(peaks)+FLANK.DELTA),
# 					 ranges <- IRanges(start=start(peaks)-FLANK.DELTA,
# 														 end=end(peaks)),
# 					 ranges <- IRanges(start=start(peaks),
# 														 end=end(peaks)+FLANK.DELTA))

# 		neg.peak.cvg <- viewApply(Views(neg.cvg,ranges),as.numeric)
# 		pos.peak.cvg <- viewApply(Views(pos.cvg,ranges),as.numeric)

# 		lags <- mapply(function(x,y) {
# 									 cc=ccf(x,y,lag.max=100,plot=FALSE)
# 									 with(cc,lag[which.max(acf)])
# 		}, x=neg.peak.cvg, y=pos.peak.cvg)

# 		ok <- (lags > min.shift)
# 		if(!any(ok)) next
# 		ov <- ov[ok,,drop=FALSE]
# 		peaks <- peaks[ok]
# 		n.peaks <- length(peaks)
# 		PEAKS[[type]][[direction]][[i]] <<- peaks

# 		info1 <- PEAK.INFO[[type]][[1]][[i]][ov[,1],]
# 		info2 <- PEAK.INFO[[type]][[2]][[i]][ov[,2],]
# 		peak.locs <- as.integer(round((info1$PEAK.LOC + info2$PEAK.LOC)/2))

# 		peak.cvg <- info1$PEAK.CVG + info2$PEAK.CVG
# 		peak.surL <- info1$PEAK.SURL + info2$PEAK.SURL
# 		peak.surR <- info1$PEAK.SURR + info2$PEAK.SURR

# 		switch(i,
# 					 {zscores <- zscore(peak.cvg, peak.surL+peak.surR,2)
# 						max.z <- zscore(peak.cvg+peak.surL+peak.surR,0,2)},
# 					 {zscores <- zscore(peak.cvg, peak.surL)
# 						max.z <- zscore(peak.cvg+peak.surL,0)},
# 					 {zscores <- zscore(peak.cvg, peak.surR)
# 						max.z <- zscore(peak.cvg+peak.surR,0)})

# 		peak.nlps <- -pnorm(zscores, low=FALSE, log=TRUE)/log(10)
# 		max.nlps <- -pnorm(max.z, low=FALSE, log=TRUE)/log(10)

# 		dfr <- data.frame(NLP=peak.nlps, MAX.NLP=max.nlps, LOC=peak.locs,
# 											WIDTH=width(peaks), START=start(peaks), END=end(peaks),
# 											CVG=peak.cvg, SURL=peak.surL, SURR=peak.surR, FORM=i)

# 		rownames(dfr) <- with(dfr,sprintf("%s:%d-%d:%d",CHR,START,END,FORM))
# 		PEAK.INFO[[type]][[direction]][[i]] <<- dfr
# 		N.PEAKS <<- N.PEAKS + n.peaks
# 	}
# }
