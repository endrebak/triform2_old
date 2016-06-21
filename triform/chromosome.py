from sys import argv
import argparse
import os
import pkg_resources
from subprocess import call


def chromosome(infile, is_input, args):

    chromosome_script = pkg_resources.resource_filename("triform",
                                                        "R/chromosome.R")
    all_vars = vars(args)
    all_vars["chromosome_script"] = chromosome_script
    all_vars["infile"] = infile
    all_vars["is_input"] = str(is_input).upper()
    command = "Rscript {chromosome_script} {infile} {is_input} {min_shift} {min_width} {min_enrichment} {min_z} {flank_distance}".format(
        **all_vars)
    print(command)
    call(command, shell=True)
