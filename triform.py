from collections import defaultdict
from itertools import product
from sys import argv
import argparse
import os
from os.path import basename, join
import pkg_resources
from subprocess import call

from joblib import Parallel, delayed

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


def parallel_make_ranged_data(in_out_map, args):

    return Parallel(n_jobs=args.number_cores)(delayed(make_ranged_data)(
        infile, outfile) for (infile, outfile) in in_out_map.items())


def parallel_make_chromosome_cover_file(in_out_map2, args):

    for infile in in_out_map2:
        Parallel(n_jobs=args.number_cores)(delayed(make_chromosome_cover_file)(
            infile, outfiles, args.read_width)
                                           for (infile, outfiles) in product(
                                               [infile], in_out_map2[infile]))


def add_chromosome_to_filename(files, in_out_map):

    outfiles2 = defaultdict(list)

    for infile in files:
        infile2 = in_out_map[infile]
        path_start, path_end = infile2.rsplit("/", 1)
        for chromosome in ["chr22"]:
            chr_p_file = join(path_start,
                              "_".join([chromosome, "FORWARD", path_end]))
            chr_m_file = join(path_start,
                              "_".join([chromosome, "REVERSE", path_end]))
            outfiles2[infile2].append(",".join([chr_p_file, chr_m_file]))

    return outfiles2


if __name__ == '__main__':
    args = parser.parse_args()
    print("# triform " + " ".join(argv[1:]))
    chip = args.treatment
    input = args.control

    if args.tmpdir == "/tmp/{process_id}":
        pid = str(os.getpid())
        outpath = join("/tmp/", pid)
    else:
        outpath = args.tmpdir

    call("mkdir -p {}".format(outpath), shell=True)

    in_out_map = dict((infile, join(outpath, basename(infile) + ".RData"))
                      for infile in chip + input)
    # parallel_make_ranged_data(in_out_map, args)

    in_out_map2 = add_chromosome_to_filename(chip + input, in_out_map)

    print(in_out_map.keys())
    chip_files = in_out_map2[in_out_map[chip[0]]][0]
    input_files = in_out_map2[in_out_map[input[0]]][0]
    # chromosome(chip_files, input_files, args)
    # print(in_out_map2)
    # parallel_make_chromosome_cover_file(in_out_map2, args)
    # print(in_out_map2)

    # chromosome_files = {}

    # chip_files = defaultdict(list)
    # for infile in chip:
    #     parallel_make_chromosome_cover_file(infile, args)

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
