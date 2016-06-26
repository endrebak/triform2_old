from rpy2 import robjects as ro
from rpy2.robjects import r, pandas2ri
ri2py = pandas2ri.ri2py
from rpy2.robjects.packages import importr
importr("GenomicRanges")


def _init(covers, is_input, args):

    init_function = r("""test.init <- function(chr, is.input, flank.delta) {

    flank.delta.pad <- Rle(0, flank.delta)

    DIRECTIONS <- factor(1:2, lab=c("REVERSE","FORWARD"))
    N.DIRS <- length(DIRECTIONS)

    LOCATIONS <- factor(1:3, lab=c("LEFT","RIGHT","CENTER"))
    N.LOCS <- length(LOCATIONS)

    N.DIRLOCS <- N.DIRS*N.LOCS

    DIRECTION <- rep(DIRECTIONS, 3)
    LOCATION <- rep(LOCATIONS, 2)

    CVG <<- list()
    SIZES <<- NULL
    SIZES <<- c(SIZES, rep(unlist(chr$SIZE),ea=N.LOCS))
    print(SIZES)
    cvgs <- chr$CVG		# coverage on each strand

    h = 1
    for (i in 1:N.DIRLOCS) {   			# DIRECTON.LOCATION index
      n <- i+N.DIRLOCS*(h-1)   			# SAMPLE.DIRECTON.LOCATION index
      j <- ceiling(i/N.LOCS)    		# DIRECTION index
      cvg <- cvgs[[j]]        			# strand-specific coverage

      if(is.input) {
          if(IS.CENTER[n]) {
              CVG[[n]] <<- cvg
          }
          next												# no need for flanking input coverage
      }
      switch(1 + (i-1)%%N.LOCS,			# LOCATION index
              CVG[[n]] <<- c(flank.delta.pad, rev(rev(cvg)[-1:-flank.delta])),		# strand-specific coverage on left flank
              CVG[[n]] <<- c(cvg[-1:-flank.delta], flank.delta.pad),							#	strand-specific coverage on right flank
              CVG[[n]] <<- cvg																										# strand-specific coverage on center
      )
    }
    names(CVG) <<- paste(DIRECTION,LOCATION,sep=".")
    maxlen <- max(sapply(CVG,length))

    names(SIZES) <<- names(CVG)

    CVG <<- lapply(CVG,function(cvg) c(cvg,Rle(0,maxlen-length(cvg))))
    }""")

    return init_function(covers, is_input, args.flank_distance)
