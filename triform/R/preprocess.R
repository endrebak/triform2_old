


## makeChromosomeCoverFiles <- function(filePath=".", filePattern="RData$", gapped.width=100, outputPath="./chrcovers"){
##   rdfiles <- list.files(path=filePath, pattern=filePattern)
##   for (file in rdfiles){
##     message("Making chromosome coverage files for ", file)
##     load(file.path(filePath, file))
##     for (chr in names(rd)){
##       if (grepl("[_M]", chr)) next       # Ignore chromosome M
##       strand <- values(rd)[[chr]]$strand
##       y <- ranges(rd)[[chr]]
##       ly <- split(y, strand)
##       lsize <- lapply(ly, length)
##       lstart <- lapply(ly, start)
##       lwidth <- lapply(ly, width)
##       lgap <- lapply(lwidth, function(w) floor((gapped.width-w)/2))
##       irpos <- IRanges(start=lstart[["+"]] - lgap[["+"]], width=lwidth[["+"]] + 2*lgap[["+"]])
##       irneg <- IRanges(start=lstart[["-"]] - lgap[["-"]], width=lwidth[["-"]] + 2*lgap[["-"]])
##       irl <- IRangesList("-"=irneg, "+"=irpos)
##       lcvg <- coverage(irl)
##       covers <- list(SIZE=lsize, CVG=lcvg)

##       print(outputPath)
##       print(file)
##       save(covers, file=file.path(outputPath, paste(chr, "_", file, sep="")))
##     }
##   }
## }
