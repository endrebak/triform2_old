suppressMessages(library(GenomicRanges))
args <- commandArgs(TRUE)

infile = args[1]
outfile = args[2]

makeRangedData <- function(infile){
  options(stringsAsFactors=FALSE)
  dfr <- read.delim(infile, header=FALSE,
                    colClasses=c("character", rep("integer",2), rep("NULL",2), "character"))
  colnames(dfr) <- c("seqnames", "start", "end", "strand")
  print(head(dfr))
  rd = makeGRangesFromDataFrame(dfr)
  save(rd, file=outfile)
}

makeRangedData(infile, outfile)
