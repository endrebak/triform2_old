from subprocess import check_output
from collections import OrderedDict

from joblib import Parallel, delayed

from triform.config.genomes import create_genome_size_dict


def bed_to_chromosome_dfs(bed_file, args):

    chromosomes = create_genome_size_dict(args.genome)

    text_dfs = Parallel(n_jobs=args.number_cores)(
        delayed(_bed_to_chromosomes)(bed_file, chromosome)
        for chromosome in chromosomes)

    text_dfs = OrderedDict([(
        c, df) for (c, df) in zip(chromosomes, text_dfs) if df != ""])
    return text_dfs


def _bed_to_chromosomes(bed_file, chromosome):
    command = "grep -E '^{chromosome}\\b' {bed_file} | cut -f 1-3,6".format(
        chromosome=chromosome,
        bed_file=bed_file)
    output = check_output(command, shell=True)

    return output
