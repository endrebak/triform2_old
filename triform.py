from sys import argv
import argparse
import os
import pkg_resources
from subprocess import call

from triform.version import __version__
from triform.preprocess import (make_ranged_data, make_chromosome_cover_file)

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
    '--min-lag',
    '-ml',
    required=False,
    default=10,
    type=int,
    help=
    '''Minimum inter-strand lag (shift) between peak coverage distributions (default 10 bp).''')

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

parser.add_argument('--tmpdir',
                    '-td',
                    required=False,
                    default="/tmp/",
                    type=str,
                    help='''Directory to store intermediate results in.''')

parser.add_argument('--version',
                    '-v',
                    action='version',
                    version='%(prog)s {}'.format(__version__))

if __name__ == '__main__':
    args = parser.parse_args()
    print("# triform " + " ".join(argv[1:]))
    infile = args.treatment[0]
    make_ranged_data(infile, "testing.rds")
    make_chromosome_cover_file("testing.rds", "deleteme.rds", "chr22", 100)

    # chromosome_script = pkg_resources.resource_filename("triform",
    #                                                     "R/chromosome.R")
    # call(
    #     ("Rscript {chromosome_script} {infile} {tempdir}"
    #      "{chromosome} {min_z}"
    #      " {min_shift} {min_width}").format(
    #          chromosome_script=chromosome_script,
    #          infile=infile,
    #          tempdir=args.tmpdir,
    #          chromosome="chrY",
    #          min_z=args.min_z,
    #          min_shift=args.min_lag,
    #          min_width=args.min_width),
    #     shell=True)
