from collections import OrderedDict

import pkg_resources

from natsort import natsorted


def get_genome_size_file(genome):
    return pkg_resources.resource_filename(
        "triform", "scripts/chromsizes/{}.chromsizes".format(genome))


def create_genome_size_dict(genome):
    """Creates genome size dict from string containing data."""

    size_file = get_genome_size_file(genome)
    size_lines = (l.split() for l in open(size_file).readlines())

    size_dict = OrderedDict()
    for chromosome, length in natsorted(size_lines):
        size_dict[chromosome] = int(length)

    return size_dict
