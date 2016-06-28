import pandas as pd
import pytest
from rpy2.robjects import r, pandas2ri
from rpy2.robjects.packages import importr
importr("GenomicRanges")
importr("S4Vectors")

from triform.init import _init


@pytest.fixture
def input_data(run_length_encodings):

    df_to_rle = r("""function(df){
    Rle(df$values, df$lengths)
    }
    """)

    run_lengths = r["list"]()
    for name, rle in run_length_encodings.items():
        df = r["read.table"](rle, sep=" ", header=1)
        run_lengths.rx2[name] = df_to_rle(df)

    return run_lengths


@pytest.fixture
def expected_result():
    return None


@pytest.mark.current
def test_chromosome(input_data, args, expected_result):

    result = chromosome(input_data, args)
    assert result == expected_result


"""
# Hoel test
zscore <- function(x,y,r=1) {  # r = size.y/size.x
  dif <- (r*x-y)
  zs <- dif/sqrt(r*(x+y))
  zs[!dif] <- 0
  zs
}
"""

# MAX.P <- 0.1
# MIN.Z <- qnorm(MAX.P,low=FALSE)
# MIN.SHIFT <- 10
# MIN.WIDTH <- 10


def chromosome(chromosome_covers, args):
    """produces PEAKS, PEAK.INFO and CENTER.CVG
    PEAK.INFO: Six dataframes, reverse/forward * type peak (1, 2, 3)
    CENTER.CVG: Rle forward backward for the center coverage

    The original loop ran once and included both replicates
    We want to add another loop for replicates

    Idea: create nec scaffolding to run once, store result

    To run function below must have list with regular and input.
    Must create necessary list.

    Have one input list, one chip list.

    Idea: create p1, p2, p3 and zscores1,zscores2,zscores3
      Then do "for(i in 1:3) {  # separate analyses of different peak forms"
      as own function

    """

    chromosome_function = r(
        """function(chr, min.z=MIN.Z, min.shift=MIN.SHIFT, min.width=MIN.WIDTH) {
  # test.init(chr)

  DIRECTIONS <- factor(1:2, lab=c("REVERSE","FORWARD"))
  LOCATIONS <- factor(1:3, lab=c("LEFT","RIGHT","CENTER"))
  PEAKS <<- list()
  PEAK.INFO <<- list()
  CENTER.CVG <<- list()
  N.PEAKS <<- 0

    PEAKS <<- list()
    PEAK.INFO <<- list()
    CENTER.CVG <<- list()

    for(direction in DIRECTIONS){
      PEAKS[[direction]] <<- list(IRanges(),IRanges(),IRanges())
      PEAK.INFO[[direction]] <<- list(NULL,NULL,NULL)

      is.dir <- grepl(direction, CVG.NAMES)
      # ref <- (CVG[[CVG.NAMES[IS.INPUT & is.dir & IS.CENTER & IS.REP1]]] +
      #         CVG[[CVG.NAMES[IS.INPUT & is.dir & IS.CENTER & IS.REP2]]])

      # surL1 <- CVG[[CVG.NAMES[is.type & is.dir & IS.LEFT & IS.REP1]]]
      # surR1 <- CVG[[CVG.NAMES[is.type & is.dir & IS.RIGHT & IS.REP1]]]
      # cvg1 <- CVG[[CVG.NAMES[is.type & is.dir & IS.CENTER & IS.REP1]]]

      # surL2 <- CVG[[CVG.NAMES[is.type & is.dir & IS.LEFT & IS.REP2]]]
      # surR2 <- CVG[[CVG.NAMES[is.type & is.dir & IS.RIGHT & IS.REP2]]]
      # cvg2 <- CVG[[CVG.NAMES[is.type & is.dir & IS.CENTER & IS.REP2]]]

      surL <- surL1 + surL2
      surR <- surR1 + surR2
      cvg <- cvg1 + cvg2
      CENTER.CVG[[direction]] <<- cvg

      # Form-1 test with consistency check
      signs1 <- sign(2*cvg1-surL1-surR1)
      signs2 <- sign(2*cvg2-surL2-surR2)
      ok <- (signs1==1)*(signs2==1)
      zscores1 <- zscore(cvg,surL+surR,2) * ok
      peaks1 <- slice(zscores1,lower=min.z)
      peaks1 <- peaks1[width(peaks1)>min.width]
      p1 <- as(peaks1,"IRanges")

      # Form-2 test with consistency check
      signs1 <- sign(cvg1-surL1)
      signs2 <- sign(cvg2-surL2)
      ok <- (signs1==1)*(signs2==1)
      zscores2 <- zscore(cvg,surL) * ok
      peaks2 <- slice(zscores2,lower=min.z)
      peaks2 <- peaks2[width(peaks2)>min.width]
      p2 <- as(peaks2,"IRanges")

      # Form-3 test with consistency check
      signs1 <- sign(cvg1-surR1)
      signs2 <- sign(cvg2-surR2)
      ok <- (signs1==1)*(signs2==1)
      zscores3 <- zscore(cvg,surR) * ok
      peaks3 <- slice(zscores3,lower=min.z)
      peaks3 <- peaks3[width(peaks3)>min.width]
      p3 <- as(peaks3,"IRanges")

      # enrichment test with consistency check
      ref.size <- (SIZES[IS.INPUT & is.dir & IS.CENTER & IS.REP1] +
                   SIZES[IS.INPUT & is.dir & IS.CENTER & IS.REP2])
      cvg1.size <- SIZES[is.type & is.dir & IS.CENTER & IS.REP1]
      cvg2.size <- SIZES[is.type & is.dir & IS.CENTER & IS.REP2]
      cvg.size <- cvg1.size + cvg2.size
      ratio1 <- ref.size/cvg1.size
      ratio2 <- ref.size/cvg2.size
      ratio <- ref.size/cvg.size

      signs1 <- sign(ratio1*cvg1-ref)
      signs2 <- sign(ratio2*cvg2-ref)
      ok <- (signs1==1)*(signs2==1)
      zscores4 <- zscore(cvg,ref,ratio) * ok
      peaks4 <- slice(zscores4,lower=min.z)
      peaks4 <- peaks4[width(peaks4)>min.width]
      p4 <- as(peaks4,"IRanges")

      p1 <- intersect(p1,p4)
      ok <- (width(p1)>min.width)
      p1 <- p1[ok]

      p2 <- intersect(p2,p4)
      ok <- (width(p2)>min.width)
      p2 <- p2[ok]

      p3 <- intersect(p3,p4)
      ok <- (width(p3)>min.width)
      p3 <- p3[ok]

      peaks.list <- list(p1,p2,p3)
      zscores.list <- list(zscores1,zscores2,zscores3)
      zviews.list <- mapply(function(x,y) Views(x,y),
                            x=zscores.list, y=peaks.list)
      maxz.list <- lapply(zviews.list, viewMaxs)

      for(i in 1:3) {  # separate analyses of different peak forms
        peaks <- peaks.list[[i]]
        if(!length(peaks)) next

        maxz <- maxz.list[[i]]
        peak.nlps <- -pnorm(maxz, low=FALSE,log=TRUE)/log(10)

        peak.locs <- round((start(peaks)+end(peaks))/2)
        peak.cvg <- cvg[peak.locs,drop=TRUE]
        peak.ref <- ref[peak.locs,drop=TRUE]
        peak.enrich <- (1+ratio*peak.cvg)/(1+peak.ref)

        if(i==1) min.er <<- quantile(peak.enrich,MIN.QUANT)
        ok <- (peak.enrich>min.er)
        if(!any(ok)) next

        peaks <- peaks[ok]
        peak.locs <- peak.locs[ok]
        peak.nlps <- peak.nlps[ok]

        peak.cvg <- peak.cvg[ok]
        peak.surL <- surL[peak.locs,drop=TRUE]
        peak.surR <- surR[peak.locs,drop=TRUE]

        PEAKS[[direction]][[i]] <<- peaks
        n.peaks <- length(peaks)
        dfr <- data.frame(PEAK.LOC=peak.locs, PEAK.CVG=peak.cvg,
                          PEAK.SURL=peak.surL, PEAK.SURR=peak.surR,
                          PEAK.NLP=round(peak.nlps,3), PEAK.WIDTH=width(peaks),
                          PEAK.START=start(peaks), PEAK.END=end(peaks))
        PEAK.INFO[[direction]][[i]] <<- dfr
      }
    }
}""")

    chromosome_function(chromosome_covers, args.min_z, args.min_shift,
                        args.min_width)
