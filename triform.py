from sys import argv
import argparse
import os

from triform.version import __version__

parser = argparse.ArgumentParser(
    description="""Diffuse domain ChIP-Seq caller based on SICER.

(Visit github.com/endrebak/epic for examples and help.)

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
