pkgname <- "triform"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('triform')

assign(".oldSearch", search(), pos = 'CheckExEnv')
cleanEx()
nameEx("preprocess")
### * preprocess

flush(stderr()); flush(stdout())

### Name: preprocess
### Title: Preprocess BED-files with reads to calculate chromosome coverage
### Aliases: preprocess
### Keywords: manip

### ** Examples

## Not run:
# Run preprocessing using a configuration file in YAML format:
## print(getwd())
## source("/Users/labsenter/havpryd/code/triform/R/preprocess.R")
## preprocess("/Users/labsenter/havpryd/code/triform/inst/extdata/config.yml")

# Run preprocessing without using a configuration file:
## preprocess(configPath = NULL, params=list(READ.PATH="./tmp",
##            COVER.PATH="./chrcovers", READ.WIDTH=100))
## End(Not run)



cleanEx()
nameEx("triform-package")
### * triform-package

flush(stderr()); flush(stdout())

### Name: triform-package
### Title: Triform finds peaks in ChIP-sequencing data.
### Aliases: triform-package
### Keywords: package

### ** Examples

## Not run:
##D preprocess("./config.yml")
source("/Users/labsenter/havpryd/code/triform/R/triform.R")
triform("/Users/labsenter/havpryd/code/triform/inst/extdata/config.yml", params=list(hoo="hoo"))
## End(Not run)



cleanEx()
nameEx("triform")
### * triform

flush(stderr()); flush(stdout())

### Name: triform
### Title: Run Triform peak detection algorithm.
### Aliases: triform
### Keywords: model

### ** Examples


## Not run:
##D # Run Triform with configuration file:
##D triform(configPath = "./config.yml")
##D
##D # Run Triform with arguments instead of configuration file:
##D triform(configPath=NULL, params=list(COVER.PATH = "./chrcovers",
##D        OUTPUT.PATH = "./Triform_output.csv",
##D        CONTROLS = c("backgr_huds_Gm12878_rep1", "backgr_huds_Gm12878_rep2"),
##D        TARGETS = c("srf_huds_Gm12878_rep1", "srf_huds_Gm12878_rep2"),
##D        MAX.P = 0.1, MIN.WIDTH = 10, MIN.QUANT = 0.375, MIN.SHIFT = 10,
##D        FLANK.DELTA = 150, CHRS = c("chrY")))
## End(Not run)



### * <FOOTER>
###
cat("Time elapsed: ", proc.time() - get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
# quit('no')
