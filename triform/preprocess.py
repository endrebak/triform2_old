from sys import argv
import argparse
import os
import pkg_resources
from subprocess import call


def make_ranged_data(infile, outfile):

    make_ranged_data_script = pkg_resources.resource_filename(
        "triform", "R/makeRangedData.R")
    command = "Rscript {make_ranged_data_script} {infile} {outfile}".format(
        make_ranged_data_script=make_ranged_data_script,
        infile=infile,
        outfile=outfile)
    print(command)
    call(command, shell=True)


def make_chromosome_cover_file(infile, outfile, gapped_width):

    make_chromosome_cover_files_script = pkg_resources.resource_filename(
        "triform", "R/makeChromosomeCoverFiles.R")
    command = "Rscript {make_chromosome_cover_files_script} {infile} {outfile} {gapped_width}".format(
        make_chromosome_cover_files_script=make_chromosome_cover_files_script,
        infile=infile,
        outfile=outfile,
        gapped_width=gapped_width)
    print(command)
    call(command, shell=True)
