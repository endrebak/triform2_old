
options(stringsAsFactors=FALSE)
suppressMessages(library(GenomicRanges))
args <- commandArgs(TRUE)

infile = args[1]
outfile = args[2]
gapped.width = strtoi(args[3])

makeChromosomeCoverFiles <- function(infile, outfile, gapped.width){
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

    save(covers, file=outfile)
}

makeChromosomeCoverFiles(infile, outfile, gapped.width)
