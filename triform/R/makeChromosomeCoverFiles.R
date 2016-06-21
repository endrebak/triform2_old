
options(stringsAsFactors=FALSE)
suppressMessages(library(IRanges))
args <- commandArgs(TRUE)

infile = args[1]
outfile = args[2]
chr = args[3]
gapped.width = strtoi(args[4])

makeChromosomeCoverFiles <- function(infile, outfile, chr, gapped.width){
    if (grepl("[_M]", chr)) stop()       # Ignore chromosome M
    rd <- get(load(infile))

    strands <- strand(rd)
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

    save(covers, file=outfile)
}

makeChromosomeCoverFiles(infile, outfile, chr, gapped.width)
