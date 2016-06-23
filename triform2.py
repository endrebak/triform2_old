# from collections import defaultdict
# from itertools import product
# from sys import argv
# import argparse
# import os
# from os.path import basename, join
import pkg_resources
from subprocess import check_output
from io import BytesIO

import pandas as pd
from joblib import Parallel, delayed
import rpy2
from rpy2.robjects.packages import importr

from triform.version import __version__
from triform.preprocess import (make_ranged_data, make_chromosome_cover_file)
from triform.chromosome import chromosome

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


def bed_to_chromosomes(bed_file, chromosome, args):

    # chromosomes = pkg_utils.

    text_dfs = Parallel(n_jobs=args.number_cores)(
        delayed(_bed_to_chromosomes)(bed_file, chromosome)
        for chromosome in chromosomes)


def _bed_to_chromosomes(bed_file, chromosome):
    command = "{grep} -E '^{chromosome}\\b' {bed_file} | cut -f 1-3,6".format(
        chromosome=chromosome,
        bed_file=bed_file)
    output = check_output(command, shell=True)

    return BytesIO(output)


if __name__ == '__main__':
    args = parser.parse_args()
    print("# triform2 " + " ".join(argv[1:]))
