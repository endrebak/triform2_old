import pkg_resources

from collections import defaultdict, OrderedDict
import os
import argparse
from sys import argv
import pandas as pd
from joblib import Parallel, delayed
import rpy2
from rpy2.robjects.packages import importr
from rpy2.robjects import r

from triform.version import __version__
from triform.preprocess.preprocess import preprocess
from triform.init import init
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

parser.add_argument('--genome',
                    '-g',
                    required=False,
                    default="hg19",
                    type=str,
                    help='''Genome version to use.''')

parser.add_argument(
    '--max-p',
    '-mp',
    required=False,
    default=0.1,
    type=float,
    help=
    '''Used to calculate minimum upper-tail z-value (default corresponds to standard normal p = 0.1)''')

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

if __name__ == '__main__':
    args = parser.parse_args()
    print("# triform2 " + " ".join(argv[1:]))

    treatment, control, treatment_sizes, control_sizes = preprocess(args)

    print(treatment_sizes, control_sizes)
    for k, v in treatment.items():
        print(k)
        print(v)
        print(type(v))
    init_treatment = init(treatment, False, args)
    init_control = init(control, True, args)
    print(init_treatment.keys())
    chromosome(init_treatment, init_control, treatment_sizes, control_sizes,
               args)

    print(init_treatment.keys(), "init_treatment.keys()")
    print(init_control.keys(), "init_control.keys()")

    # for k, v in init_treatment.items():
    #     print(k)
    #     print(v)
    # print(init_treatment)

    # rpy2.robjects.r["save"](treatment[0], file="chrcovers.RData")
