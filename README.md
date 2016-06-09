# triform

Attempt at reviving the triform algorithm, perhaps the best punctuate peak finder (think transcription factors).

This is done in my free time (which I do not have much of), so progress will be leisurely.

#### TODO

The list will grow as I discover more stuff that needs to be done. Feel free to add your own.

- Take input/chip files as command line arguments (easy)
- Do not require input files to follow a naming convention (easy)
- Do not use temp-files, keep data in memory (to speed up the algorithm and avoid file paths in the code) (medium)
- Run chromosomes in parallel (unknown difficulty)
- Accept an arbitrary number of chip and input-files (hard?)

#### Install

Download the latest version of R. Then install bioconductor with:

```R
source("https://bioconductor.org/biocLite.R")
biocLite()
```

Download the package sources for yaml [here](https://cran.rstudio.com/web/packages/yaml/index.html).

Then install it using:

```R
install.packages("yaml_<version>.tgz", repos=NULL, type="source")
```

Now you should be able to run the triform example with

```
source("triform-Ex.R")
```

#### Paper

[The Triform algorithm: improved sensitivity and specificity in ChIP-Seq peak finding](http://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-13-176)
