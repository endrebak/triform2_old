
options(stringsAsFactors=FALSE)
suppressMessages(library(GenomicRanges))
args <- commandArgs(TRUE)

infile = args[1]
outfile = args[2]
gapped.width = strtoi(args[3])

files = strsplit(outfile, ",")[[1]]
forward.outfile = files[[1]]
reverse.outfile = files[[2]]

makeChromosomeCoverFiles <- function(infile, forward.outfile, reverse.outfile, gapped.width){
    rd <- get(load(infile))

    # hack to avoid * levels
    strands <- factor(as.vector(strand(rd)), levels=c("+", "-"))

    ly <- split(rd, strands)
    lsize <- lapply(ly, length)
    lstart <- lapply(ly, start)
    lwidth <- lapply(ly, width)

    lgap <- lapply(lwidth, function(w) floor((gapped.width-w)/2))
    irpos <- IRanges(start=lstart[["+"]] - lgap[["+"]], width=lwidth[["+"]] + 2*lgap[["+"]])
    irneg <- IRanges(start=lstart[["-"]] - lgap[["-"]], width=lwidth[["-"]] + 2*lgap[["-"]])
    irl <- IRangesList("-"=irneg, "+"=irpos)
    lcvg <- coverage(irl)
    covers <- list(SIZE=lsize, CVG=lcvg)

    test.init(covers)
    ## forward = covers$CVG"+"
    ## reverse = covers$CVG"-"
    ## save(forward, file=forward.outfile)
    ## save(reverse, file=reverse.outfile)
}

makeChromosomeCoverFiles(infile, forward.outfile, reverse.outfile, gapped.width)


test.init <- function(chrcovers) {
  ## chrcovers <- get(load(infile))

  CVG <<- list()

  # what this loop does:
  # create coverage for product(["+", "-"], ["left", "centre", "right"))
  # can replace with for
  # for (i in c("+", "-")){ for (j in c("f", "c", "b")) {print(paste(i, j))}}
  SIZES <<- list()
  for (i in 1:N.DIRLOCS) {   			# DIRECTON.LOCATION index
    n <- i	# SAMPLE.DIRECTON.LOCATION index
    j <- ceiling(i/N.LOCS)    		# DIRECTION index
    cvg <- chrcovers$CVG[[j]]        			# strand-specific coverage
    SIZES <<- rep(chrcovers$SIZE, N.LOCS)

    if(is.input) {
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
  names(CVG) <<- CVG.NAMES

  maxlen <- max(sapply(CVG,length))
  CVG <<- lapply(CVG,function(cvg) c(cvg,Rle(0,maxlen-length(cvg))))
  names(SIZES) <<- CVG.NAMES
  save(CVG, file="CVG.RData")
  save(SIZES, file="SIZES.RData")
}
