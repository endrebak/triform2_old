from subprocess import call
"""
PROBABLY best to have R files that are part of distribution,
get with pkg_data, call from here
"""

make_ranged_data_template = """
library(IRanges)
dfr <- read.delim("{}", header=FALSE,
colClasses=c("character", rep("integer",2), rep("NULL",2), "character"))
colnames(dfr) <- c("space", "start", "end", "strand")
dfr <- with(dfr, dfr[order(space, start, end), ])
rd <- as(dfr, "RangedData")
save(rd, file="{}")"""


def make_ranged_data(infile, outfile):

    cmd = make_ranged_data_template.format(infile, outfile)
    call("")
