from sys import argv
import argparse
import os
import pkg_resources
from subprocess import call


def chromosome(chipfiles, inputfiles, args):

    chromosome_script = pkg_resources.resource_filename("triform",
                                                        "R/chromosome.R")
    all_vars = vars(args)
    all_vars["chromosome_script"] = chromosome_script
    all_vars["chipfiles"] = chipfiles
    all_vars["inputfiles"] = inputfiles
    command = "Rscript {chromosome_script} {chipfiles} {inputfiles} {min_shift} {min_width} {min_enrichment} {min_z} {flank_distance}".format(
        **all_vars)
    print(command)
    call(command, shell=True)
