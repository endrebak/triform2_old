# triform2

ChIP-Seq caller for transcription factor binding sites.

Very alpha edition, so please post questions, issues, docrequests etc at [the triform2 google groups forum.](https://groups.google.com/forum/#!forum/triform2)

#### Install

Clone repo and run with `python triform2.py -h`

#### Papers

* [The Triform algorithm: improved sensitivity and specificity in ChIP-Seq peak finding](http://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-13-176)
* [ChIP analysis unravels an exceptionally wide distribution of DNA binding sites for the NtcA transcription factor in a heterocyst-forming cyanobacterium](http://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-15-22)

#### CLI

```
usage: triform2.py [-h] --treatment TREATMENT [TREATMENT ...] --control
                   CONTROL [CONTROL ...] [--number-cores NUMBER_CORES]
                   [--genome GENOME] [--bedgraph BEDGRAPH] [--max-p MAX_P]
                   [--min-shift MIN_SHIFT] [--min-width MIN_WIDTH]
                   [--read-width READ_WIDTH] [--flank-distance FLANK_DISTANCE]
                   [--min-enrichment MIN_ENRICHMENT] [--version]

Improved sensitivity, specificity and control of false discovery rates in
ChIP-Seq peak finding. (Visit github.com/endrebak/triform for examples and
help.)

optional arguments:
  -h, --help            show this help message and exit
  --treatment TREATMENT [TREATMENT ...], -t TREATMENT [TREATMENT ...]
                        Treatment (pull-down) file(s) in
                        bam/bed/bed.gz/bed.bz2 format.
  --control CONTROL [CONTROL ...], -c CONTROL [CONTROL ...]
                        Control (input) file(s) in bam/bed/bed.gz/bed.bz2
                        format.
  --number-cores NUMBER_CORES, -cpu NUMBER_CORES
                        Number of cpus to use. Can use at most one per
                        chromosome. Default: 1.
  --genome GENOME, -g GENOME
                        Genome version to use.
  --bedgraph BEDGRAPH, -b BEDGRAPH
                        Path to write bedgraph file to, if desired.
  --max-p MAX_P, -mp MAX_P
                        Used to calculate minimum upper-tail z-value (default
                        corresponds to standard normal p = 0.1)
  --min-shift MIN_SHIFT, -ms MIN_SHIFT
                        Minimum inter-strand shift (lag) between peak coverage
                        distributions (default 10 bp).
  --min-width MIN_WIDTH, -mw MIN_WIDTH
                        Minimum number of bp (peak width) in peak-like region
                        (default 10 bp).
  --read-width READ_WIDTH, -rw READ_WIDTH
                        Read width w, symmetrically extended to a fixed value.
                        Must be larger than the flank size. Default: 100 bp.
  --flank-distance FLANK_DISTANCE, -fd FLANK_DISTANCE
                        Fixed spacing between central and flanking locations
                        (must be > w). Default: 150 bp.
  --min-enrichment MIN_ENRICHMENT, -mr MIN_ENRICHMENT
                        Minimum local enrichment ratio (default 3/8 quantile
                        of the enrichment ratio)
  --version, -v         show program's version number and exit
```

#### TODO

- only require a certain percentage of ChIP files to show a peak?
- ensure chromosomes operated on in correct order
- merge type 2/3 peaks
- create setup.py
- add cl-arg for merge distance
- add to travis ci
- add to bioconda
- create wrapper for asciigenome to print peaks
- create galaxy wrapper
- add to coveralls
- also make bedgraph of peaks before merging on strand? (ask advisors if this will be valuable)

#### Done

- Make possible to use without input data
- Suite of unit-tests
- Implemented FDR calculations
- Take input/chip files as command line arguments (easy)
- Do not require input files to follow a naming convention (easy)
- Do not use temp-files, keep data in memory (to speed up the algorithm and avoid file paths in the code) (medium)
- Run chromosomes in parallel (unknown difficulty)
- Accept an arbitrary number of chip and input-files (hard?)
