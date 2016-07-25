
options(stringsAsFactors=FALSE)
library(IRanges)

MAX.P <- 0.1
MIN.Z <- qnorm(MAX.P,low=FALSE)

MIN.SHIFT <- 10
MIN.WIDTH <- 10
MIN.QUANT <- 3/8
FLANK.DELTA <- 150
FLANK.SHIFT <- 2*FLANK.DELTA
FLANK.DELTA.PAD <- Rle(0,FLANK.DELTA)

CHRS <- paste("chr",c(1:22,"X","Y"),sep="")

SUMCVG.NAMES <- c("backgr_huds_Gm12878_rep1","backgr_huds_Gm12878_rep2",
									"srf_huds_Gm12878_rep1","srf_huds_Gm12878_rep2")
TARGET.NAMES <- "srf_huds_Gm12878"

N.TYPES <- length(SUMCVG.NAMES)
TYPES <- factor(1:N.TYPES, lab=SUMCVG.NAMES)

DIRECTIONS <- factor(1:2, lab=c("REVERSE","FORWARD"))
N.DIRS <- length(DIRECTIONS)

LOCATIONS <- factor(1:3, lab=c("LEFT","RIGHT","CENTER"))
N.LOCS <- length(LOCATIONS)

N.DIRLOCS <- N.DIRS*N.LOCS

DIRECTION <- rep(rep(DIRECTIONS,ea=N.LOCS),N.TYPES)
LOCATION <- rep(LOCATIONS,N.DIRS*N.TYPES)
TYPE <- rep(TYPES,ea=N.DIRLOCS)

CVG.NAMES <- paste(TYPE,DIRECTION,LOCATION,sep=".")
IS.LEFT <- grepl("LEFT",CVG.NAMES)
IS.RIGHT <- grepl("RIGHT",CVG.NAMES)
IS.CENTER <- grepl("CENTER",CVG.NAMES)
IS.REP1 <- grepl("_rep1",CVG.NAMES)
IS.REP2 <- grepl("_rep2",CVG.NAMES)
IS.INPUT <- !grepl("srf",CVG.NAMES)

test.init <- function(chr) {
  if(!exists("CHR",inherits=TRUE)) CHR <<- "none"
  if(chr==CHR) return()

  load(paste("ChrCovers/",chr,".RData",sep=""), .GlobalEnv)

  CVG <<- list()
  SIZES <<- NULL
  for(h in 1:N.TYPES) {       			# TYPE index
    type <- SUMCVG.NAMES[h]
    SIZES <<- c(SIZES, rep(unlist(chrcovers[[type]]$SIZE),ea=N.LOCS))
    cvgs <- chrcovers[[type]]$CVG		# coverage on each strand

    for (i in 1:N.DIRLOCS) {   			# DIRECTON.LOCATION index
      n <- i+N.DIRLOCS*(h-1)   			# SAMPLE.DIRECTON.LOCATION index
      j <- ceiling(i/N.LOCS)    		# DIRECTION index
      cvg <- cvgs[[j]]        			# strand-specific coverage

			if(IS.INPUT[n]) {
				if(IS.CENTER[n]) {
					CVG[[n]] <<- cvg
				}
				next												# no need for flanking input coverage
			}
      switch(1 + (i-1)%%N.LOCS,			# LOCATION index
             CVG[[n]] <<- c(FLANK.DELTA.PAD, rev(rev(cvg)[-1:-FLANK.DELTA])),		# strand-specific coverage on left flank
             CVG[[n]] <<- c(cvg[-1:-FLANK.DELTA], FLANK.DELTA.PAD),							#	strand-specific coverage on right flank
             CVG[[n]] <<- cvg																										# strand-specific coverage on center
      )
    }
  }
  names(CVG) <<- CVG.NAMES
  maxlen <- max(sapply(CVG,length))
  CVG <<- lapply(CVG,function(cvg) c(cvg,Rle(0,maxlen-length(cvg))))
  names(SIZES) <<- CVG.NAMES

  CHR <<- chr
}

# Hoel test
zscore <- function(x,y,r=1) {  # r = size.y/size.x
  dif <- (r*x-y)
  zs <- dif/sqrt(r*(x+y))
  zs[!dif] <- 0
  zs
}

test.chr <- function(chr,
                     min.z=MIN.Z,
										 min.shift=MIN.SHIFT,
                     min.width=MIN.WIDTH) {
  test.init(chr)

  PEAKS <<- list()
  PEAK.INFO <<- list()
	CENTER.CVG <<- list()
	N.PEAKS <<- 0

	for(type in TARGET.NAMES)  {
    PEAKS[[type]] <<- list()
    PEAK.INFO[[type]] <<- list()
    CENTER.CVG[[type]] <<- list()
		is.type <- grepl(type, CVG.NAMES)

		for(direction in DIRECTIONS){
			PEAKS[[type]][[direction]] <<- list(IRanges(),IRanges(),IRanges())
			PEAK.INFO[[type]][[direction]] <<- list(NULL,NULL,NULL)

			is.dir <- grepl(direction, CVG.NAMES)
			ref <- (CVG[[CVG.NAMES[IS.INPUT & is.dir & IS.CENTER & IS.REP1]]] +
							CVG[[CVG.NAMES[IS.INPUT & is.dir & IS.CENTER & IS.REP2]]])

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

			for(i in 1:3) {	# separate analyses of different peak forms
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

				PEAKS[[type]][[direction]][[i]] <<- peaks
				n.peaks <- length(peaks)
				dfr <- data.frame(PEAK.LOC=peak.locs, PEAK.CVG=peak.cvg,
													PEAK.SURL=peak.surL, PEAK.SURR=peak.surR,
													PEAK.NLP=round(peak.nlps,3), PEAK.WIDTH=width(peaks),
													PEAK.START=start(peaks), PEAK.END=end(peaks))
				PEAK.INFO[[type]][[direction]][[i]] <<- dfr
			}
    }
    # CENTER.CVG and PEAK.INFO is used below
		direction <- "merged"
		PEAKS[[type]][[direction]] <<- list(IRanges(),IRanges(),IRanges())
		PEAK.INFO[[type]][[direction]] <<- list(NULL,NULL,NULL)
		PEAK.INFO[[type]][["regions"]] <<- list(NULL,NULL,NULL)

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

			info1 <- PEAK.INFO[[type]][[1]][[i]][ov[,1],]
			info2 <- PEAK.INFO[[type]][[2]][[i]][ov[,2],]
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
	if(!N.PEAKS) return()

# exclude redundant Form-2 and Form-3 peaks
	p1 <- PEAKS[[type]][[direction]][[1]]
	p2 <- PEAKS[[type]][[direction]][[2]]
	p3 <- PEAKS[[type]][[direction]][[3]]

	ov12 <- matrix(as.matrix(findOverlaps(p1,p2)),ncol=2)
	if(!!nrow(ov12)) {
		ex2 <- (1:length(p2) %in% ov12[,2])
		p2 <- p2[!ex2]
		PEAKS[[type]][[direction]][[2]] <<- p2
		info <- PEAK.INFO[[type]][[direction]][[2]][!ex2,,drop=FALSE]
		PEAK.INFO[[type]][[direction]][[2]] <<- info
	}
	ov13 <- matrix(as.matrix(findOverlaps(p1,p3)),ncol=2)
	if(!!nrow(ov13)) {
		ex3 <- (1:length(p3) %in% ov13[,2])
		p3 <- p3[!ex3]
		PEAKS[[type]][[direction]][[3]] <<- p3
		info <- PEAK.INFO[[type]][[direction]][[3]][!ex3,,drop=FALSE]
		PEAK.INFO[[type]][[direction]][[3]] <<- info
	}

# exclude overlapping Form-2 and Form-3 peaks
	ov23 <- matrix(as.matrix(findOverlaps(p2,p3,maxgap=100)),ncol=2)
	if(!!nrow(ov23)) {
		ex2 <- (1:length(p2) %in% ov23[,1])
		p2 <- p2[!ex2]
		PEAKS[[type]][[direction]][[2]] <<- p2
		info <- PEAK.INFO[[type]][[direction]][[2]][!ex2,,drop=FALSE]
		PEAK.INFO[[type]][[direction]][[2]] <<- info

		ex3 <- (1:length(p3) %in% ov23[,2])
		p3 <- p3[!ex3]
		PEAKS[[type]][[direction]][[3]] <<- p3
		info <- PEAK.INFO[[type]][[direction]][[3]][!ex3,,drop=FALSE]
		PEAK.INFO[[type]][[direction]][[3]] <<- info
	}
	peak.info <- rbind(PEAK.INFO[[type]][[direction]][[1]],
										 PEAK.INFO[[type]][[direction]][[2]],
										 PEAK.INFO[[type]][[direction]][[3]])
	N.PEAKS <<- nrow(peak.info)
	peak.info
}

test.genome <- function(min.z=MIN.Z,
												min.shift=MIN.SHIFT,
												min.width=MIN.WIDTH) {
  INFO <<- NULL
  for(chr in CHRS) {
    cat("\n\tProcessing",chr,"... "); flush.console()
    # calls chromosome function for each chr, adds to list
    INFO <<- rbind(INFO, test.chr(
									 chr=chr,
									 min.z=min.z,
									 min.shift=min.shift,
									 min.width=min.width))
    cat("done\n\t\tfound",N.PEAKS,"peaks")
  }
  print(INFO)
	INFO <<- INFO[order(-INFO$NLP),]
	nlps <- unique(INFO$NLP)
	sizes <- sapply(nlps,function(nlp) sum(INFO$NLP == nlp))
	indices <- sapply(nlps,function(nlp) sum(INFO$NLP >= nlp))

	nlrs <- mapply(function(nlp, j) {
		m <- sum(INFO$MAX.NLP >= nlp)	# Tarone modification for discrete nlp
		b.y <- log10(sum((1:m)^-1))		# discrete Benjamini-Yekutieli offset
		nls <- nlp + log10(j/m)				# discrete Benjamini-Hochberg adjustment
		max(nls-b.y,0)								# discrete Benjamini-Yekutieli adjustment
	}, nlp=nlps, j=indices)

	M <- length(nlrs)
	nlqs <- numeric(M)
	for(i in 1:M) nlqs[i] <- max(nlrs[i:M])		# step-up procedure
	nlqss <- unlist(mapply(function(nlq,size) rep(nlq,size),
									nlq=nlqs, size=sizes))

	INFO <<- cbind(QVAL=10^-nlqss, NLQ=nlqss, INFO)

  cat("\n\nSaving results ... "); flush.console()
  write.table(INFO,file="Triform.srf.xls",col.names=NA,quote=FALSE,sep="\t")
	save(INFO,file="Triform.srf.info.RData")
  cat("Finished.\n\n")
}
