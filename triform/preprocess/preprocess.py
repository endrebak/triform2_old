from collections import defaultdict

from triform.preprocess.make_ranged_data import make_ranged_data
from triform.preprocess.bed_to_chromosome_dfs import bed_to_chromosome_dfs
from triform.preprocess.make_chromosome_cover_files import (
    make_chromosome_cover_files)


def preprocess(args):
    """Return {chrom: list of chromcovers} map."""
    treatment = _preprocess(args.treatment, args)
    control = _preprocess(args.control, args)
    return treatment, control


def _preprocess(files, args):

    chromosomes = set()

    chromosome_covers = defaultdict(list)
    for file in files:
        chromosome_dfs = bed_to_chromosome_dfs(file, args)
        ranged_data = make_ranged_data(chromosome_dfs, args)
        chrcovers = make_chromosome_cover_files(ranged_data, args)
        chromosome_covers[file].append(chrcovers)

        chromosomes.update(set(chrcovers.keys()))

    covers_per_chromosome = defaultdict(list)
    for file, chrcovers_list in chromosome_covers.items():
        for chrcovers in chrcovers_list:
            for chromosome, cover in chrcovers.items():
                covers_per_chromosome[chromosome].append(cover)

    return covers_per_chromosome
