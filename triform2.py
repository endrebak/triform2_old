import pkg_resources

from collections import defaultdict, OrderedDict
import os
import argparse
from sys import argv
import pandas as pd
from joblib import Parallel, delayed
import rpy2
from rpy2.robjects.packages import importr

from triform.version import __version__
from triform.preprocess.make_ranged_data import make_ranged_data
from triform.preprocess.bed_to_chromosome_dfs import bed_to_chromosome_dfs
from triform.preprocess.make_chromosome_cover_files import make_chromosome_cover_files

parser = argparse.ArgumentParser(
    description=
    """Improved sensitivity, specificity and control of false discovery rates in ChIP-Seq peak finding.

(Visit github.com/endrebak/triform for examples and help.)

    """,
    prog=os.path.basename(__file__))

parser.add_argument(
    '--treatment',
    '-t',
    required=True,
    type=str,
    nargs='+',
    help='''Treatment (pull-down) file(s) in bam/bed/bed.gz/bed.bz2 format.
                   ''')

parser.add_argument(
    '--control',
    '-c',
    required=True,
    type=str,
    nargs='+',
    help='''Control (input) file(s) in bam/bed/bed.gz/bed.bz2 format.''')

parser.add_argument(
    '--number-cores',
    '-cpu',
    required=False,
    default=1,
    type=int,
    help=
    '''Number of cpus to use. Can use at most one per chromosome. Default: 1.''')

parser.add_argument('--genome',
                    '-g',
                    required=False,
                    default="hg19",
                    type=str,
                    help='''Genome version to use.''')

parser.add_argument(
    '--min-z',
    '-mz',
    required=False,
    default=0.1,
    type=float,
    help=
    '''Minimum upper-tail z-value (default corresponds to standard normal p = 0.1)''')

parser.add_argument(
    '--min-shift',
    '-ms',
    required=False,
    default=10,
    type=int,
    help=
    '''Minimum inter-strand shift (lag) between peak coverage distributions (default 10 bp).''')

parser.add_argument(
    '--min-width',
    '-mw',
    required=False,
    default=10,
    type=int,
    help=
    '''Minimum number of bp (peak width) in peak-like region (default 10 bp).''')

parser.add_argument(
    '--read-width',
    '-rw',
    required=False,
    default=100,
    type=int,
    help=
    '''Read width w, symmetrically extended to a fixed value. Must be larger than the flank size. Default: 100 bp.''')

parser.add_argument(
    '--flank-distance',
    '-fd',
    required=False,
    default=150,
    type=int,
    help=
    '''Fixed spacing between central and flanking locations (must be > w). Default: 150 bp.''')

parser.add_argument(
    '--min-enrichment',
    '-mr',
    required=False,
    default=0.375,
    type=float,
    help=
    '''Minimum local enrichment ratio (default 3/8 quantile of the enrichment ratio)''')

parser.add_argument('--tmpdir',
                    '-td',
                    required=False,
                    default="/tmp/{process_id}/",
                    type=str,
                    help='''Directory to store intermediate results in.''')

parser.add_argument('--version',
                    '-v',
                    action='version',
                    version='%(prog)s {}'.format(__version__))


def preprocess(args):
    """Return {chrom: list of chromcovers} map."""
    treatment = _preprocess(args.treatment, args)
    control = _preprocess(args.control, args)
    return treatment, control


def _preprocess(files, args):

    chromosomes = set()

    chromosome_covers = defaultdict(list)
    for file in files:
        control_chromosome_dfs = bed_to_chromosome_dfs(file, args)
        ranged_data = make_ranged_data(control_chromosome_dfs, args)
        chrcovers = make_chromosome_cover_files(ranged_data, args)
        chromosome_covers[file].append(chrcovers)

        chromosomes.update(set(chrcovers.keys()))

    covers_per_chromosome = defaultdict(list)
    for file, chrcovers_list in chromosome_covers.items():
        for chrcovers in chrcovers_list:
            for chromosome, cover in chrcovers.items():
                covers_per_chromosome[chromosome].append(cover)

    return covers_per_chromosome


if __name__ == '__main__':
    args = parser.parse_args()
    print("# triform2 " + " ".join(argv[1:]))

    treatment, control = preprocess(args)
    rpy2.robjects.r["save"](treatment[0], file="chrcovers.RData")
