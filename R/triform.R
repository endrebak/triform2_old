##' Implementation of Triform peak detection

options(stringsAsFactors=FALSE)
library(IRanges)
library(yaml)

load("temp/cvg.RData")
load("temp/sizes.RData")
load("temp/chr.RData")

##' Runs Triform according to configuration file or given parameters.
##'
##' If configPath is NULL, all the other arguments must be supplied.
##' @title triform
##' @param configPath Path to a configuration file in YAML format, or NULL.
##' @param COVER.PATH Path for coverage files (from preprocessing).
##' @param TARGETS Filenames for TFs. Must include .bed ending and _rep1 to indicate replicate number.
##' @param CONTROLS Filenames for control signal. Must include .bed ending and _rep1 to indicate replicate number.
##' @param OUTPUT.PATH Path for output file.
##' @param MAX.P Minimum p-value, used to calculate min.z
##' @param MIN.WIDTH Minimum peak width (min.n)
##' @param MIN.QUANT Quantile of enrichment ratios, used to calculate min.er
##' @param MIN.SHIFT Minimum inter-strand lag between peak coverage distributions.
##' @param FLANK.DELTA Fixed spacing between central and flanking locations (δ).
##' @param CHRS A list of chromosomes to be processed.
##' @return
triform <- function(configPath="./config.yml", params=list()){
  if (! is.null(configPath)){
    message("Using config file ", configPath)
    config <- yaml.load_file(configPath)
    message("here!")
    for (i in seq_along(config)){
      assign(names(config)[i], config[[i]], pos=".GlobalEnv")
    }
  }


  # reads if not config given on the command line
  for (i in seq_along(params)){
    assign(names(params)[i], params[[i]], pos=".GlobalEnv")
  }

  ## Make SUMCVG.NAMES and TARGET.NAMES contain these two strings
  if(!(all(grepl(".bed", c(TARGETS, CONTROLS))) && all(grepl("_rep.", c(TARGETS,CONTROLS)))))
    stop("Error: Make sure filenames in TARGETS and CONTROLS are correct (i.e. contains .bed ending and _rep1 or _rep2 to indicate replicate number.")
  TARGETS <- sub(".bed", "", TARGETS)
  CONTROLS <- sub(".bed", "", CONTROLS)
  SUMCVG.NAMES <<- c(TARGETS, CONTROLS)
  TARGET.NAMES <<- unique(sub("_rep.", "", TARGETS))

  # Below code just sets a bunch of names

  MIN.Z <<- qnorm(MAX.P, lower.tail=FALSE)
  FLANK.DELTA.PAD <<- Rle(0, FLANK.DELTA)

  N.TYPES <<- length(SUMCVG.NAMES)
  TYPES <<- factor(1:N.TYPES, labels=SUMCVG.NAMES)

  DIRECTIONS <<- factor(1:2, labels=c("REVERSE","FORWARD"))
  N.DIRS <<- length(DIRECTIONS)

  LOCATIONS <<- factor(1:3, labels=c("LEFT","RIGHT","CENTER"))
  N.LOCS <<- length(LOCATIONS)

  N.DIRLOCS <<- N.DIRS*N.LOCS

  DIRECTION <<- rep(rep(DIRECTIONS,ea=N.LOCS),N.TYPES)
  LOCATION <<- rep(LOCATIONS,N.DIRS*N.TYPES)
  TYPE <<- rep(TYPES,ea=N.DIRLOCS)

  CVG.NAMES <<- paste(TYPE,DIRECTION,LOCATION,sep=".")
## > CVG.NAMES
##  [1] "srf_huds_Gm12878_rep1.REVERSE.LEFT"
##  [2] "srf_huds_Gm12878_rep1.REVERSE.RIGHT"
##  [3] "srf_huds_Gm12878_rep1.REVERSE.CENTER"
##  [4] "srf_huds_Gm12878_rep1.FORWARD.LEFT"
##  [5] "srf_huds_Gm12878_rep1.FORWARD.RIGHT"
##  [6] "srf_huds_Gm12878_rep1.FORWARD.CENTER"
##  [7] "srf_huds_Gm12878_rep2.REVERSE.LEFT"
##  [8] "srf_huds_Gm12878_rep2.REVERSE.RIGHT"
##  [9] "srf_huds_Gm12878_rep2.REVERSE.CENTER"
## [10] "srf_huds_Gm12878_rep2.FORWARD.LEFT"
## [11] "srf_huds_Gm12878_rep2.FORWARD.RIGHT"
## [12] "srf_huds_Gm12878_rep2.FORWARD.CENTER"
## [13] "backgr_huds_Gm12878_rep1.REVERSE.LEFT"
## [14] "backgr_huds_Gm12878_rep1.REVERSE.RIGHT"
## [15] "backgr_huds_Gm12878_rep1.REVERSE.CENTER"
## [16] "backgr_huds_Gm12878_rep1.FORWARD.LEFT"
## [17] "backgr_huds_Gm12878_rep1.FORWARD.RIGHT"
## [18] "backgr_huds_Gm12878_rep1.FORWARD.CENTER"
## [19] "backgr_huds_Gm12878_rep2.REVERSE.LEFT"
## [20] "backgr_huds_Gm12878_rep2.REVERSE.RIGHT"
## [21] "backgr_huds_Gm12878_rep2.REVERSE.CENTER"
## [22] "backgr_huds_Gm12878_rep2.FORWARD.LEFT"
## [23] "backgr_huds_Gm12878_rep2.FORWARD.RIGHT"
## [24] "backgr_huds_Gm12878_rep2.FORWARD.CENTER"

  IS.LEFT <<- grepl("LEFT",CVG.NAMES)
  IS.RIGHT <<- grepl("RIGHT",CVG.NAMES)
  IS.CENTER <<- grepl("CENTER",CVG.NAMES)
  IS.REP1 <<- grepl("_rep1",CVG.NAMES)
  IS.REP2 <<- grepl("_rep2",CVG.NAMES)
  ## IS.CONTROL is True if dataset is control/background data (not in TARGET.NAMES)
  IS.CONTROL <<- !apply(sapply(TARGET.NAMES, function(t) {grepl(t,CVG.NAMES)}), 1, any)

  test.genome(min.z=MIN.Z,
              min.shift=MIN.SHIFT,
              min.width=MIN.WIDTH,
              chromosomes=CHRS,
              chrcoversPath=COVER.PATH,
              outputFilePath=OUTPUT.PATH)

}


##' Finds peaks for all given chromosomes and outputs results
##'
##'
##' @title test.genome
##' @param min.z Minimum z.value
##' @param min.shift Minimum inter-strand lag between peak coverage distributions
##' @param min.width Minimum peak width
##' @param chromosomes The chromosomes to search in
##' @param chrcoversPath Path to chromosome coverage files
##' @param outputFilePath Path for output file with peak predictions
##' @return
test.genome <- function(min.z=MIN.Z,
                        min.shift=MIN.SHIFT,
                        min.width=MIN.WIDTH,
                        chromosomes=CHRS,
                        chrcoversPath="./chrcovers",
                        outputFilePath="./Triform_output.csv") {
  INFO <<- NULL
  for(chr in chromosomes) {
    message("Triform processing ", chr)
    flush.console()
    INFO <<- rbind(INFO, test.chr(
                                  chr=chr,
                                  min.z=min.z,
                                  min.shift=min.shift,
                                  min.width=min.width,
                                  filePath=chrcoversPath))
    message("Found ", as.character(N.PEAKS), " peaks")
  }
  message("\n\nSaving results to path: ", outputFilePath)
  flush.console()
  write.table(INFO, file=outputFilePath, col.names=NA, quote=FALSE, sep="\t")
  message("Finished.")
}


##' Finds peaks for a given chromosome
##'
##'
##' @title test.chr
##' @param chr Chromosome
##' @param min.z Minimum z.value
##' @param min.shift Minimum inter-strand lag between peak coverage distributions
##' @param min.width Minimum peak width
##' @return A list of peak loci
test.chr <- function(chr,
                     min.z=MIN.Z,
                     min.shift=MIN.SHIFT,
                     min.width=MIN.WIDTH,
                     filePath="./chrcovers") {
  ## test.init(chr, filePath)
  ## save(CVG, file="temp/cvg.RData")
  ## save(SIZES, file="temp/sizes.RData")
  ## save(CHR, file="temp/chr.RData")

  PEAKS <<- list()
  PEAK.INFO <<- list()
  CENTER.CVG <<- list()
  N.PEAKS <<- 0

  # target names only [1] "srf_huds_Gm12878" - why the loop?


  print(TARGET.NAMES)
	for(type in TARGET.NAMES)  {
    print("type")
    print(type)
    PEAKS[[type]] <<- list()
    PEAK.INFO[[type]] <<- list()
    CENTER.CVG[[type]] <<- list()
		is.type <- grepl(type, CVG.NAMES)

		for(direction in DIRECTIONS){
      print("direction")
      print(direction)
			PEAKS[[type]][[direction]] <<- list(IRanges(),IRanges(),IRanges())
			PEAK.INFO[[type]][[direction]] <<- list(NULL,NULL,NULL)

			is.dir <- grepl(direction, CVG.NAMES)
			ref <- (CVG[[CVG.NAMES[IS.CONTROL & is.dir & IS.CENTER & IS.REP1]]] +
							CVG[[CVG.NAMES[IS.CONTROL & is.dir & IS.CENTER & IS.REP2]]])

			surL1 <- CVG[[CVG.NAMES[is.type & is.dir & IS.LEFT & IS.REP1]]]
      surR1 <- CVG[[CVG.NAMES[is.type & is.dir & IS.RIGHT & IS.REP1]]]
      cvg1 <- CVG[[CVG.NAMES[is.type & is.dir & IS.CENTER & IS.REP1]]]

      surL2 <- CVG[[CVG.NAMES[is.type & is.dir & IS.LEFT & IS.REP2]]]
      surR2 <- CVG[[CVG.NAMES[is.type & is.dir & IS.RIGHT & IS.REP2]]]
      cvg2 <- CVG[[CVG.NAMES[is.type & is.dir & IS.CENTER & IS.REP2]]]

      surL <- surL1 + surL2
      surR <- surR1 + surR2
      cvg <- cvg1 + cvg2
			CENTER.CVG[[type]][[direction]] <<- cvg
      print(cvg)

			# Form-1 test with consistency check
			signs1 <- sign(2*cvg1-surL1-surR1)
			signs2 <- sign(2*cvg2-surL2-surR2)
			ok <- (signs1==1)*(signs2==1)
      zscores1 <- zscore(cvg,surL+surR,2) * ok
      peaks1 <- slice(zscores1,lower=min.z)
      ## write.table(peaks1, "peaks1_r")
      subset1 = width(peaks1)>min.width
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
			ref.size <- (SIZES[IS.CONTROL & is.dir & IS.CENTER & IS.REP1] +
									 SIZES[IS.CONTROL & is.dir & IS.CENTER & IS.REP2])
			cvg1.size <- SIZES[is.type & is.dir & IS.CENTER & IS.REP1]
			cvg2.size <- SIZES[is.type & is.dir & IS.CENTER & IS.REP2]
			cvg.size <- cvg1.size + cvg2.size
			ratio1 <- ref.size/cvg1.size
			ratio2 <- ref.size/cvg2.size
      print(ref.size)
      print(cvg.size)
			ratio <- ref.size/cvg.size
      print(ratio)

			signs1 <- sign(ratio1*cvg1-ref)
			signs2 <- sign(ratio2*cvg2-ref)
			ok <- (signs1==1)*(signs2==1)

			zscores4 <- zscore(cvg,ref,ratio) * ok
      peaks4 <- slice(zscores4,lower=min.z)
      ## write.table(peaks4, "peaks4_r")
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

      ## rle_to_df = function(rle){
      ##   cbind(runValue(rle), runLength(rle))
      ## }
      ## z1 = rle_to_df(zscores1)
      ## z2 = rle_to_df(zscores2)
      ## z3 = rle_to_df(zscores3)
      ## write.table(z1, "tests/test_results/z1.csv", sep=" ")
      ## write.table(z2, "tests/test_results/z2.csv", sep=" ")
      ## write.table(z3, "tests/test_results/z3.csv", sep=" ")

      ## write.table(p1, "tests/test_results/p1.csv", sep=" ")
      ## write.table(p2, "tests/test_results/p2.csv", sep=" ")
      ## write.table(p3, "tests/test_results/p3.csv", sep=" ")
			peaks.list <- list(p1,p2,p3)
			zscores.list <- list(zscores1,zscores2,zscores3)
			zviews.list <- mapply(function(x,y) Views(x,y),
														x=zscores.list, y=peaks.list)
      print(typeof(zviews.list))
			maxz.list <- lapply(zviews.list, viewMaxs)

			for(i in 1:3) {	# separate analyses of different peak forms
        print("i")
        print(i)
				peaks <- peaks.list[[i]]

				if(!length(peaks)) next

				maxz <- maxz.list[[i]]
				peak.nlps <- -pnorm(maxz, low=FALSE,log=TRUE)/log(10)

        start.peaks = start(peaks)
        end.peaks = end(peaks)
				peak.locs <- round((start(peaks)+end(peaks))/2)
				peak.cvg <- cvg[peak.locs,drop=TRUE]
				peak.ref <- ref[peak.locs,drop=TRUE]

				peak.enrich <- (1+ratio*peak.cvg)/(1+peak.ref)

				if(i==1) min.er <<- quantile(peak.enrich,MIN.QUANT)
				ok <- (peak.enrich>min.er)
				if(!any(ok)) next

				peaks <- peaks[ok]
				peak.locs <- peak.locs[ok]
        write.table(peak.locs, "peak_locs_r")
				peak.nlps <- peak.nlps[ok]

        print("length(ok)")
        print(length(ok))

        print("length(peak.locs)")
        print(length(peak.locs))


        print("length(peak.nlps)")
        print(length(peak.nlps))

				peak.cvg <- peak.cvg[ok]
        print("surL")
        print(surL)
        print("length(surL)")
        print(length(surL))
				peak.surL <- surL[peak.locs,drop=TRUE]
        print("peak.surL")
        print(peak.surL)
        print("length(peak.surL)")
        print(length(peak.surL))
				peak.surR <- surR[peak.locs,drop=TRUE]

				PEAKS[[type]][[direction]][[i]] <<- peaks
				n.peaks <- length(peaks)
				dfr <- data.frame(PEAK.LOC=peak.locs, PEAK.CVG=peak.cvg,
													PEAK.SURL=peak.surL, PEAK.SURR=peak.surR,
													PEAK.NLP=round(peak.nlps,3), PEAK.WIDTH=width(peaks),
													PEAK.START=start(peaks), PEAK.END=end(peaks))
				PEAK.INFO[[type]][[direction]][[i]] <<- dfr
        write.table(dfr, "dfr_r", sep=" ")
			}

    }


		direction <- "merged"
		PEAKS[[type]][[direction]] <<- list(IRanges(),IRanges(),IRanges())
		PEAK.INFO[[type]][[direction]] <<- list(NULL,NULL,NULL)
		## PEAK.INFO[[type]][["regions"]] <<- list(NULL,NULL,NULL)

		neg.cvg <- CENTER.CVG[[type]][[1]]
		pos.cvg <- CENTER.CVG[[type]][[2]]

		for (i in 1:3) {
			p1 <- PEAKS[[type]][[1]][[i]]
			p2 <- PEAKS[[type]][[2]][[i]]
			if(!length(p1) | !length(p2)) next

			ov <- matrix(as.matrix(findOverlaps(p1,p2)),ncol=2)
			if(!nrow(ov)) next

			dup1 <- (ov[,1] %in% ov[duplicated(ov[,1]),1])
			dup2 <- (ov[,2] %in% ov[duplicated(ov[,2]),2])
			is.multi <- dup1 | dup2
			if(all(is.multi)) next
			ov <- ov[!is.multi,,drop=FALSE]

			p1 <- p1[1:length(p1) %in% ov[,1]]
			p2 <- p2[1:length(p2) %in% ov[,2]]
      print(p1)
      print(p2)
      stop()
			peaks <- IRanges(start=pmin(start(p1),start(p2)),
											 end=pmax(end(p1),end(p2)))

			switch(i,
						 ranges <- IRanges(start=start(peaks)-FLANK.DELTA,
															 end=end(peaks)+FLANK.DELTA),
						 ranges <- IRanges(start=start(peaks)-FLANK.DELTA,
															 end=end(peaks)),
						 ranges <- IRanges(start=start(peaks),
															 end=end(peaks)+FLANK.DELTA))

			neg.peak.cvg <- viewApply(Views(neg.cvg,ranges),as.numeric)
			pos.peak.cvg <- viewApply(Views(pos.cvg,ranges),as.numeric)

			lags <- mapply(function(x,y) {
										 cc=ccf(x,y,lag.max=100,plot=FALSE)
										 with(cc,lag[which.max(acf)])
			}, x=neg.peak.cvg, y=pos.peak.cvg)

			ok <- (lags > min.shift)
			if(!any(ok)) next
			ov <- ov[ok,,drop=FALSE]
			peaks <- peaks[ok]
			n.peaks <- length(peaks)
			PEAKS[[type]][[direction]][[i]] <<- peaks

      print("ov[,1]")
      print(ov[,1])
      print("ov[,2]")
      print(ov[,2])

			info1 <- PEAK.INFO[[type]][[1]][[i]][ov[,1],]
			info2 <- PEAK.INFO[[type]][[2]][[i]][ov[,2],]
      print("info1")
      print(info1)
      print("info2")
      print(info2)
      stop()
      stop("Endre Endre Endre")
			peak.locs <- as.integer(round((info1$PEAK.LOC + info2$PEAK.LOC)/2))

			peak.cvg <- info1$PEAK.CVG + info2$PEAK.CVG
			peak.surL <- info1$PEAK.SURL + info2$PEAK.SURL
			peak.surR <- info1$PEAK.SURR + info2$PEAK.SURR

			switch(i,
						 {zscores <- zscore(peak.cvg, peak.surL+peak.surR,2)
							max.z <- zscore(peak.cvg+peak.surL+peak.surR,0,2)},
						 {zscores <- zscore(peak.cvg, peak.surL)
							max.z <- zscore(peak.cvg+peak.surL,0)},
						 {zscores <- zscore(peak.cvg, peak.surR)
							max.z <- zscore(peak.cvg+peak.surR,0)})

			peak.nlps <- -pnorm(zscores, low=FALSE, log=TRUE)/log(10)
			max.nlps <- -pnorm(max.z, low=FALSE, log=TRUE)/log(10)

			dfr <- data.frame(NLP=peak.nlps, MAX.NLP=max.nlps, LOC=peak.locs,
												WIDTH=width(peaks), START=start(peaks), END=end(peaks),
												CVG=peak.cvg, SURL=peak.surL, SURR=peak.surR, FORM=i)

			rownames(dfr) <- with(dfr,sprintf("%s:%d-%d:%d",CHR,START,END,FORM))
			PEAK.INFO[[type]][[direction]][[i]] <<- dfr
			N.PEAKS <<- N.PEAKS + n.peaks
		}
	}
  print(PEAK.INFO)
  stop("Endre Endre Endre")
}
  ## for(type in TARGET.NAMES)  {
  ##   PEAKS[[type]] <<- list()
  ##   PEAK.INFO[[type]] <<- list()
  ##   CENTER.CVG[[type]] <<- list()
  ##   is.type <- grepl(type, CVG.NAMES)



  ##   print(DIRECTIONS)
  ##   for(direction in DIRECTIONS){
  ##     # why three? type 1, type 2, type 3?
  ##     PEAKS[[type]][[direction]] <<- list(IRanges(),IRanges(),IRanges())
  ##     PEAK.INFO[[type]][[direction]] <<- list(NULL,NULL,NULL)

  ##     if (any(IS.REP2)){
  ##       ## Use 2 replicates
  ##       formsList <- findForms2Replicates(direction, is.type, type,
  ##                                         min.z, min.width)
  ##     } else{
  ##       ## Use 1 replicate
  ##       formsList <- findForms1Replicate(direction, is.type, type,
  ##                                        min.z, min.width)
  ##     }
  ##     p1 <- formsList$p1
  ##     p2 <- formsList$p2
  ##     p3 <- formsList$p3
  ##     p4 <- formsList$p4
  ##     ratio <- formsList$ratio
  ##     ref <- formsList$ref
  ##     zscores.list <- formsList$zscores.list
  ##     cvg <- formsList$cvg
  ##     rm(formsList)


  ##     p1 <- intersect(p1,p4)
  ##     ok <- (width(p1)>min.width)
  ##     p1 <- p1[ok]

  ##     p2 <- intersect(p2,p4)
  ##     ok <- (width(p2)>min.width)
  ##     p2 <- p2[ok]

  ##     p3 <- intersect(p3,p4)
  ##     ok <- (width(p3)>min.width)
  ##     p3 <- p3[ok]

  ##     peaks.list <- list(p1,p2,p3)
  ##     zviews.list <- mapply(function(x,y) Views(x,y),
  ##                           x=zscores.list, y=peaks.list)
  ##     maxz.list <- lapply(zviews.list, viewMaxs)

  ##     for(i in 1:3) {	# separate analyses of different peak forms
  ##       peaks <- peaks.list[[i]]
  ##       if(!length(peaks)) next

  ##       maxz <- maxz.list[[i]]
  ##       peak.nlps <- -pnorm(maxz, lower.tail=FALSE,log.p=TRUE)/log(10)

  ##       peak.locs <- round((start(peaks)+end(peaks))/2)
  ##       peak.cvg <- cvg[peak.locs,drop=TRUE]
  ##       peak.ref <- ref[peak.locs,drop=TRUE]
  ##       peak.enrich <- (1+ratio*peak.cvg)/(1+peak.ref)

  ##       if(i==1) min.er <<- quantile(peak.enrich,MIN.QUANT)
  ##       ok <- (peak.enrich>min.er)
  ##       if(!any(ok)) next

  ##       peaks <- peaks[ok]
  ##       peak.locs <- peak.locs[ok]
  ##       peak.nlps <- peak.nlps[ok]

  ##       PEAKS[[type]][[direction]][[i]] <<- peaks
  ##       n.peaks <- length(peaks)
  ##       dfr <- data.frame(PEAK.LOC=peak.locs,
  ##                         PEAK.NLP=round(peak.nlps,3), PEAK.WIDTH=width(peaks),
  ##                         PEAK.START=start(peaks), PEAK.END=end(peaks))
  ##       PEAK.INFO[[type]][[direction]][[i]] <<- dfr
  ##     }
  ##   }
##     direction <- "merged"
##     PEAKS[[type]][[direction]] <<- list(IRanges(),IRanges(),IRanges())
##     PEAK.INFO[[type]][[direction]] <<- list(NULL,NULL,NULL)
##     PEAK.INFO[[type]][["regions"]] <<- list(NULL,NULL,NULL)

##     neg.cvg <- CENTER.CVG[[type]][[1]]
##     pos.cvg <- CENTER.CVG[[type]][[2]]

##     for (i in 1:3) {
##       p1 <- PEAKS[[type]][[1]][[i]]
##       p2 <- PEAKS[[type]][[2]][[i]]
##       if(!length(p1) | !length(p2)) next

##       ov <- matrix(as.matrix(findOverlaps(p1,p2)),ncol=2)
##       if(!nrow(ov)) next

##       dup1 <- (ov[,1] %in% ov[duplicated(ov[,1]),1])
##       dup2 <- (ov[,2] %in% ov[duplicated(ov[,2]),2])
##       is.multi <- dup1 | dup2
##       if(all(is.multi)) next
##       ov <- ov[!is.multi,,drop=FALSE]

##       p1 <- p1[1:length(p1) %in% ov[,1]]
##       p2 <- p2[1:length(p2) %in% ov[,2]]
##       peaks <- IRanges(start=pmin(start(p1),start(p2)),
##                        end=pmax(end(p1),end(p2)))

##       switch(i,
##              ranges <- IRanges(start=start(peaks)-FLANK.DELTA,
##                                end=end(peaks)+FLANK.DELTA),
##              ranges <- IRanges(start=start(peaks)-FLANK.DELTA,
##                                end=end(peaks)),
##              ranges <- IRanges(start=start(peaks),
##                                end=end(peaks)+FLANK.DELTA))

##       neg.peak.cvg <- viewApply(Views(neg.cvg,ranges),as.numeric)
##       pos.peak.cvg <- viewApply(Views(pos.cvg,ranges),as.numeric)

##       lags <- mapply(function(x,y) {
##         cc=ccf(x,y,lag.max=100,plot=FALSE)
##         with(cc,lag[which.max(acf)])
##       }, x=neg.peak.cvg, y=pos.peak.cvg)

##       ok <- (lags > min.shift)
##       if(!any(ok)) next
##       ov <- ov[ok,,drop=FALSE]
##       peaks <- peaks[ok]

##       if(i==1) type.delta <<- round(median(lags[ok]))

##       info1 <- PEAK.INFO[[type]][[1]][[i]][ov[,1],]
##       info2 <- PEAK.INFO[[type]][[2]][[i]][ov[,2],]
##       peak.locs <- round((info1$PEAK.LOC + info2$PEAK.LOC)/2)
##       peak.nlps <- info1$PEAK.NLP + info2$PEAK.NLP

##       PEAKS[[type]][[direction]][[i]] <<- peaks
##       n.peaks <- length(peaks)

##       dfr <- data.frame(PEAK.FORM=i, PEAK.NLP=peak.nlps,
##                         PEAK.WIDTH=width(peaks), PEAK.LOC=peak.locs,
##                         PEAK.START=start(peaks), PEAK.END=end(peaks))

##       rownames(dfr) <- with(dfr,sprintf("%s:%d-%d:%d",CHR,PEAK.START,PEAK.END,PEAK.FORM))
##       PEAK.INFO[[type]][[direction]][[i]] <<- dfr
##       N.PEAKS <<- N.PEAKS + n.peaks
##     }
##   }
##   if(!N.PEAKS) return()

##                                         # exclude redundant Form-2 and Form-3 peaks
##   p1 <- PEAKS[[type]][[direction]][[1]]
##   p2 <- PEAKS[[type]][[direction]][[2]]
##   p3 <- PEAKS[[type]][[direction]][[3]]
##   peak.info <- PEAK.INFO[[type]][[direction]][[1]]

##   ov12 <- matrix(as.matrix(findOverlaps(p1,p2)),ncol=2)
##   if(!!nrow(ov12)) {
##     ex2 <- (1:length(p2) %in% ov12[,2])
##     p2 <- p2[!ex2]
##     PEAKS[[type]][[direction]][[2]] <<- p2
##     info <- PEAK.INFO[[type]][[direction]][[2]][!ex2,,drop=FALSE]
##     PEAK.INFO[[type]][[direction]][[2]] <<- info
##     peak.info <- rbind(peak.info,info)
##   }
##   ov13 <- matrix(as.matrix(findOverlaps(p1,p3)),ncol=2)
##   if(!!nrow(ov13)) {
##     ex3 <- (1:length(p3) %in% ov13[,2])
##     p3 <- p3[!ex3]
##     PEAKS[[type]][[direction]][[3]] <<- p3
##     info <- PEAK.INFO[[type]][[direction]][[3]][!ex3,,drop=FALSE]
##     PEAK.INFO[[type]][[direction]][[3]] <<- info
##     peak.info <- rbind(peak.info,info)
##   }

##                                         # merge overlapping Form-2 and Form-3 peaks into Form-1 peaks
##   peak.info <- with(peak.info,peak.info[order(PEAK.START,PEAK.END),])
##   rng <- with(peak.info,IRanges(start=PEAK.START,end=PEAK.END))
##   ov <- matrix(as.matrix(findOverlaps(rng,maxgap=1,drop.self=TRUE,drop.redundant=TRUE)),ncol=2)
##   if(!!nrow(ov)) {
##     peak.info[ov[,1],"PEAK.FORM"] <- 1
##     peak.info[ov[,1],"PEAK.LOC"] <- round((peak.info[ov[,1],"PEAK.LOC"] + peak.info[ov[,2],"PEAK.LOC"])/2)
##     peak.info[ov[,1],"PEAK.NLP"] <- peak.info[ov[,1],"PEAK.NLP"] + peak.info[ov[,2],"PEAK.NLP"]
##     peak.info[ov[,1],"PEAK.START"] <- pmin(peak.info[ov[,1],"PEAK.START"],peak.info[ov[,2],"PEAK.START"])
##     peak.info[ov[,1],"PEAK.END"] <- pmax(peak.info[ov[,1],"PEAK.END"],peak.info[ov[,2],"PEAK.END"])
##     peak.info[ov[,1],"PEAK.WIDTH"] <- 1 + peak.info[ov[,1],"PEAK.END"] - peak.info[ov[,1],"PEAK.START"]
##     peak.info <- peak.info[-ov[,2],]
##   }

##   N.PEAKS <<- nrow(peak.info)
##   peak.info
## }



## ##' Creates sample/direction/location-specific coverage data for a chromosome
## ##'
## ##' Stores result in global variable CVG
## ##' @title test.init
## ##' @param chr The chromosome
## ##' @param filePath The path to the chromosome coverage file
## ##' @return
## test.init <- function(chr, filePath="./chrcovers") {
##   # Checking if CHR exists in the environment
##   if(!exists("CHR",inherits=TRUE)) CHR <<- "none"
##   if(chr==CHR) return()

##   ## print(file.path(filePath, paste(chr,".RData",sep=""))) # "./chrcovers/chrY.RData"
##   load(file.path(filePath, paste(chr,".RData",sep="")), .GlobalEnv) # load chrcovers

##   CVG <<- list()
##   SIZES <<- NULL
##   # N.TYPES is the number of files
##   for(h in 1:N.TYPES) {       			# TYPE index
##     type <- SUMCVG.NAMES[h] # "srf_huds_Gm12878_rep1"
##     # what is c(SIZES, rep(unlist(chrcovers[[type]]$SIZE),ea=N.LOCS)), the below?
##     #    -    -    -    +    +    +
##     # 8688 8688 8688 8986 8986 8986
##     # each file has chromosome coverage - what does it consists of?
##     # length 100 regions and the regions inbetween them
##     print(chrcovers[[type]])
##     SIZES <<- c(SIZES, rep(unlist(chrcovers[[type]]$SIZE),ea=N.LOCS))

##     ## RleList of length 2
##     ## $`-`
##     ## integer-Rle of length 57442690 with 12472 runs
##     ##   Lengths: 2709676     100   31062     100 ...     232     100      94     100
##     ##   Values :       0       1       0       1 ...       0       1       0       1
##     ## ## $`+`
##     ## ## integer-Rle of length 57442498 with 12854 runs
##     ## ##   Lengths: 2709647      29      71      29 ...      56     100     161     100
##     ## ##   Values :       0       1       2       1 ...       0       1       0       1
##     cvgs <- chrcovers[[type]]$CVG		# coverage on each strand

##     for (i in 1:N.DIRLOCS) {   			# DIRECTON.LOCATION index
##       print("i")
##       print(i)
##       n <- i+N.DIRLOCS*(h-1)   			# SAMPLE.DIRECTON.LOCATION index
##       print("n")
##       print(n)
##       j <- ceiling(i/N.LOCS)    		# DIRECTION index
##       print("j")
##       print(j)
##       cvg <- cvgs[[j]]        			# strand-specific coverage
##       print(cvg)

##       if(IS.CONTROL[n]) {
##         if(IS.CENTER[n]) {
##           CVG[[n]] <<- cvg
##         }
##         next                                    # no need for flanking input coverage
##       }

##       # This decides whether left, centre or right flank
##       # N.
##       switch(1 + (i-1)%%N.LOCS,			# LOCATION index
##              CVG[[n]] <<- c(FLANK.DELTA.PAD, rev(rev(cvg)[-1:-FLANK.DELTA])), # strand-specific coverage on left flank
##              CVG[[n]] <<- c(cvg[-1:-FLANK.DELTA], FLANK.DELTA.PAD), # strand-specific coverage on right flank
##              CVG[[n]] <<- cvg     # strand-specific coverage on center
##              )
##       print(CVG[[1 + (i-1)%%N.LOCS]])
##     }
##   }
##   # loop over "srf_huds_Gm12878_rep1" etc done
##   # have collected coverage data for each strand
##   # what is it now used for?

##   names(CVG) <<- CVG.NAMES
##   maxlen <- max(sapply(CVG,length))
##   print(CVG)
##   CVG <<- lapply(CVG,function(cvg) c(cvg,Rle(0,maxlen-length(cvg))))
##   names(SIZES) <<- CVG.NAMES

##   CHR <<- chr
## }


##' Calculates Z score
##'
##' @title zscore
##' @param x Signal
##' @param y Background
##' @param r Ratio (background size / signal size)
##' @return z-score
zscore <- function(x,y,r=1) {  # r = size.y/size.x
  dif <- (r*x-y)
  zs <- dif/sqrt(r*(x+y))
  zs[!dif] <- 0
  zs
}






findForms1Replicate <- function(direction, is.type, type,
                                min.z, min.width){
  is.dir <- grepl(direction, CVG.NAMES)
  ref <- (CVG[[CVG.NAMES[IS.CONTROL & is.dir & IS.CENTER & IS.REP1]]])


  surL <- CVG[[CVG.NAMES[is.type & is.dir & IS.LEFT & IS.REP1]]]
  surR <- CVG[[CVG.NAMES[is.type & is.dir & IS.RIGHT & IS.REP1]]]
  cvg <- CVG[[CVG.NAMES[is.type & is.dir & IS.CENTER & IS.REP1]]]

  CENTER.CVG[[type]][[direction]] <<- cvg


  signs <- sign(2*cvg-surL-surR)
  ok <- (signs==1)
  zscores1 <- zscore(cvg,surL+surR,2) * ok
  peaks1 <- slice(zscores1,lower=min.z)
  peaks1 <- peaks1[width(peaks1)>min.width]
  p1 <- as(peaks1,"IRanges")


  signs <- sign(cvg-surL)
  ok <- (signs==1)
  zscores2 <- zscore(cvg,surL) * ok
  peaks2 <- slice(zscores2,lower=min.z)
  peaks2 <- peaks2[width(peaks2)>min.width]
  p2 <- as(peaks2,"IRanges")


  signs <- sign(cvg-surR)
  ok <- (signs==1)
  zscores3 <- zscore(cvg,surR) * ok
  peaks3 <- slice(zscores3,lower=min.z)
  peaks3 <- peaks3[width(peaks3)>min.width]
  p3 <- as(peaks3,"IRanges")


  ref.size <- (SIZES[IS.CONTROL & is.dir & IS.CENTER & IS.REP1])
  cvg.size <- SIZES[is.type & is.dir & IS.CENTER & IS.REP1]
  ratio <- ref.size/cvg.size

  signs <- sign(ratio*cvg-ref)
  ok <- (signs==1)
  zscores4 <- zscore(cvg,ref,ratio) * ok
  peaks4 <- slice(zscores4,lower=min.z)
  peaks4 <- peaks4[width(peaks4)>min.width]
  p4 <- as(peaks4,"IRanges")

  zscores.list <- list(zscores1,zscores2,zscores3)

  return(list(p1=p1, p2=p2, p3=p3, p4=p4, ratio=ratio, ref=ref, zscores.list=zscores.list, cvg=cvg))

}


findForms2Replicates <- function(direction, is.type, type,
                                 min.z, min.width){
  ## print(direction) # FORWARD or REVERSE
  is.dir <- grepl(direction, CVG.NAMES)
  ## print(CVG.NAMES)
  ## ref <- (CVG[[CVG.NAMES[IS.CONTROL & is.dir & IS.CENTER & IS.REP1]]] +
  ##         CVG[[CVG.NAMES[IS.CONTROL & is.dir & IS.CENTER & IS.REP2]]])
  print("CVG[which(IS.CONTROL & is.dir & IS.CENTER)]")
  print(CVG[which(IS.CONTROL & is.dir & IS.CENTER)])

  ref <- Reduce("+",CVG[which(IS.CONTROL & is.dir & IS.CENTER)])	# no need for replicate inputs
  # does this ^^ mean they average the inputs?
  print("ref")
  print(ref)

  # what does sur stand for?
  surL1 <- CVG[[CVG.NAMES[is.type & is.dir & IS.LEFT & IS.REP1]]]
  surR1 <- CVG[[CVG.NAMES[is.type & is.dir & IS.RIGHT & IS.REP1]]]
  cvg1 <- CVG[[CVG.NAMES[is.type & is.dir & IS.CENTER & IS.REP1]]]

  surL2 <- CVG[[CVG.NAMES[is.type & is.dir & IS.LEFT & IS.REP2]]]
  surR2 <- CVG[[CVG.NAMES[is.type & is.dir & IS.RIGHT & IS.REP2]]]
  cvg2 <- CVG[[CVG.NAMES[is.type & is.dir & IS.CENTER & IS.REP2]]]

  # adding the two datasets together
  surL <- surL1 + surL2
  surR <- surR1 + surR2
  cvg <- cvg1 + cvg2
  # why only storing the center?
  CENTER.CVG[[type]][[direction]] <<- cvg

                                        # Form-1 test with consistency check
  # w
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
  ## ref.size <- (SIZES[IS.CONTROL & is.dir & IS.CENTER & IS.REP1] +
  ##              SIZES[IS.CONTROL & is.dir & IS.CENTER & IS.REP2])
  ref.size <- sum(SIZES[IS.CONTROL & is.dir & IS.CENTER])				# no need for replicate inputs
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

  zscores.list <- list(zscores1,zscores2,zscores3)

  return(list(p1=p1, p2=p2, p3=p3, p4=p4, ratio=ratio, ref=ref, zscores.list=zscores.list, cvg=cvg))
}
