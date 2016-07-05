##' Pre-processes BED files with tag data
##' Makes coverage files (1 pr chromosome) describing tag counts along the chromosomes.
##'

options(stringsAsFactors=FALSE)
library(IRanges)
library(yaml)


##' Runs preprocessing according to configuration file.
##'
##' If configPath is NULL, the readPath, coverPath and readWidth are taken from the function arguments.
##' @title preprocess
##' @param configPath Path to a configuration file in YAML format (see config.yml), or NULL.
##' @param readPath Path to BED files with reads.
##' @param coverPath Path for coverage files (preprocessing output).
##' @param readWidth Read width (w)
##' @return
preprocess <- function(configPath="./config.yml", params=list()){
  if (! is.null(configPath)){
    message("Using config file ", configPath)
    config <- yaml.load_file(configPath)
    READ.PATH <- config$READ.PATH
    COVER.PATH <- config$COVER.PATH
    READ.WIDTH <- config$READ.WIDTH
  }

  for (i in seq_along(params)){
    assign(names(params)[i], params[[i]])
  }

  if (! file.exists(READ.PATH))
    stop("Error! READ.PATH does not exist: ", READ.PATH)
  if (! file.exists(COVER.PATH)){
    warning("COVER.PATH: '", COVER.PATH, "' does not exist. Creating directory.")
    dir.create(COVER.PATH, recursive=TRUE)
  }

  makeRangedData(READ.PATH)
  makeChromosomeCoverFiles(filePath=READ.PATH, gapped.width=READ.WIDTH,
                           outputPath=COVER.PATH)
  mergeChromosomeCoverFiles(filePath=COVER.PATH)
}



##' Converts BED files with position of mapped reads to IRanges RangedData objects
##'
##'
##' @title makeRangedData
##' @param filePath The file path for the BED files
##' @return
makeRangedData <- function(filePath="."){
  files <- list.files(path = filePath, pattern = "bed$")
  if (length(files)==0)
    warning("Cannot find any BED files in path: ", filePath)
  for (file in files){
    message("Converting BED to RangedData for file ", file.path(filePath, file))
    dfr <- read.delim(file.path(filePath, file), header=FALSE,
      colClasses=c("character", rep("integer",2), rep("NULL",2), "character"))
    colnames(dfr) <- c("space", "start", "end", "strand")
    dfr <- with(dfr, dfr[order(space, start, end), ])
    write.table(dfr, file=file.path(filePath, sub("bed$", "csv", file)), sep=" ")
    rd <- as(dfr, "RangedData")
    save(rd, file=file.path(filePath, sub("bed$", "RData", file)))
  }
}

## filePath = "inst/extdata/srf_huds_Gm12878_rep1.bed"
## dfr <- read.delim(file.path(filePath, file), header=FALSE,
##       colClasses=c("character", rep("integer",2), rep("NULL",2), "character"))
## colnames(dfr) <- c("space", "start", "end", "strand")
## rd <- as(dfr, "RangedData")
## RangedData with 17674 rows and 1 value column across 1 space
##          space               ranges   |      strand
##       <factor>            <IRanges>   | <character>
## 1         chrY [10592905, 10592930]   |           -
## 2         chrY [12968115, 12968140]   |           -
## 3         chrY [11926812, 11926837]   |           -

##' Makes files describing tag coverage for each chromosome
##'
##' Each file is split on chromosomes and the coverage pr strand is calculated in intervals
##' @title makeChromosomeCoverFiles
##' @param filePath The file path for the files
##' @param filePattern The name pattern for RangedData files, ending with RData
##' @param gapped.width The width of each coverage interval
##' @return
makeChromosomeCoverFiles <- function(filePath=".", filePattern="RData$", gapped.width=100, outputPath="./chrcovers"){
  rdfiles <- list.files(path=filePath, pattern=filePattern)
  for (file in rdfiles){
    message("Making chromosome coverage files for ", file)
    load(file.path(filePath, file))
    for (chr in names(rd)){
      if (grepl("[_M]", chr)) next       # Ignore chromosome M
      strand <- values(rd)[[chr]]$strand
      y <- ranges(rd)[[chr]]
      ly <- split(y, strand)
      lsize <- lapply(ly, length)
      lstart <- lapply(ly, start)
      lwidth <- lapply(ly, width)
      lgap <- lapply(lwidth, function(w) floor((gapped.width-w)/2))
      irpos <- IRanges(start=lstart[["+"]] - lgap[["+"]], width=lwidth[["+"]] + 2*lgap[["+"]])
      irneg <- IRanges(start=lstart[["-"]] - lgap[["-"]], width=lwidth[["-"]] + 2*lgap[["-"]])
      irl <- IRangesList("-"=irneg, "+"=irpos)
      lcvg <- coverage(irl)
      covers <- list(SIZE=lsize, CVG=lcvg)
      print(outputPath)
      print(covers)
      save(covers, file=file.path(outputPath, paste(chr, "_", file, sep="")))
    }
  }
}


##' Merges/pools chromosome coverage files pr chromosome
##'
##' Makes one file pr chromosome, containing a named list of all files merged.
##' Each list item has the coverage for each strand, according to that file.
##' @title mergeChromosomeCoverFiles
##' @param filePath Path to coverage files
##' @param filePattern File name pattern of coverage files
##' @return
mergeChromosomeCoverFiles <- function(filePath="./chrcovers", filePattern="RData$"){
  allfiles <- list.files(filePath, pattern=filePattern, full.names=TRUE)
  chrfiles <- allfiles[grepl(".+/chr([^_]+)_", allfiles)]
  chrs <- unique(sub(".+/([^_]+)_.+", "\\1", chrfiles))

  for (chr in chrs){
    message("Merging coverage files for chromosome ", chr)
    chrcovers <- list()
    files <- chrfiles[grep(paste(chr, "_", sep=""), chrfiles)]
    for (file in files){
      cat(file, '\n')
      load(file)
      cat("Loaded", file, '\n')
      target <- sub("[^_]+_(.+).RData", "\\1", file)
      print(file)
      print(target)
      chrcovers[[target]] <- covers
    }
    save(chrcovers, file=file.path(filePath, paste(chr, ".RData", sep="")))
  }
}
