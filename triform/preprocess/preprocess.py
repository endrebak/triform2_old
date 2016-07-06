from collections import defaultdict

from triform.preprocess.make_ranged_data import make_ranged_data
from triform.preprocess.bed_to_chromosome_dfs import bed_to_chromosome_dfs
from triform.preprocess.make_chromosome_cover_files import (
    make_chromosome_cover_files)


def preprocess(args):
    """Return {chrom: list of chromcovers} map."""
    treatment, treatment_sizes = _preprocess(args.treatment, args)
    control, control_sizes = _preprocess(args.control, args)
    return treatment, control, treatment_sizes, control_sizes


def _preprocess(files, args):

    ranged_data_per_file = {}
    for file in files:
        chromosome_dfs = bed_to_chromosome_dfs(file, args)
        ranged_data_per_file[file] = make_ranged_data(chromosome_dfs, args)

    # print(ranged_data_per_file.keys(), 'ranged_data_per_file.keys()')
    chrcovers, sizes = make_chromosome_cover_files(ranged_data_per_file, args)

    py_sizes = defaultdict(dict)
    for c, d in sizes.items():
        for k, d2 in d.items():
            py_sizes[c][k] = int(d2[0])

    return chrcovers, py_sizes
