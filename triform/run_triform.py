import logging
from sys import stdout

from triform.helper_functions import subset_RS4
from triform.chromosome import chromosome
from triform.find_peaks import find_peaks
from triform.preprocess.preprocess import preprocess
from triform.compute_fdr import compute_fdr
from triform.init import init_background, init_treatment
from triform.make_treatment_control_same_length import (
    make_treatment_control_same_length)
from triform.exclude_redundant_peaks import exclude_redundant_peaks
from triform.create_bigwig import create_bigwig


def run_triform_no_control(args):

    logging.info("Preprocessing bed files.")
    treatment, treatment_sizes = preprocess(args)

    logging.info("Initializing treatment data.")
    init_chip = init_treatment(treatment, args)

    logging.info("Computing region statistics.")
    results = chromosome(init_chip, None, treatment_sizes, None, args)

    logging.info("Finding enriched peaks.")
    peaks = find_peaks(results, args)

    logging.info("Excluding redundant peaks.")
    peak_info = exclude_redundant_peaks(peaks, args)

    logging.info("Computing FDR.")
    fdr_table = compute_fdr(peak_info)

    fdr_table.to_csv(stdout, sep=" ")

    if args.bedgraph:
        logging.info("Writing bedgraph to file " + args.bedgraph + ".")
        create_bigwig(init_chip, fdr_table, args)

    logging.info("Done.")


def run_triform(args):

    logging.info("Preprocessing bed files.")
    treatment, control, treatment_sizes, control_sizes = preprocess(args)

    logging.info("Initializing treatment data.")
    init_chip = init_treatment(treatment, args)

    logging.info("Initializing background data.")
    init_control = init_background(control, args)

    init_chip, init_control = make_treatment_control_same_length(init_chip,
                                                                 init_control)

    # for k, v in init_chip.items():
    #     print(k)
    #     for k2, v2 in init_chip.items():
    #         print(k2)
    logging.info("Computing region statistics.")
    results = chromosome(init_chip, init_control, treatment_sizes,
                         control_sizes, args)
    # print("results")
    # for k, v in results.items():
    #     print(k)
    #     for k2, v2 in v.items():
    #         print("  " + str(k2))
    #         try:
    #             print(v2.keys())
    #         except:
    #             pass

    logging.info("Finding enriched peaks.")
    peaks = find_peaks(results, args)

    # print("find_peaks")
    # for k, v in peaks.items():
    #     print(k)
    #     try:
    #         # print(v.keys())
    #         for k2, v2 in v.items():
    #             print("  " + str(k2))
    #             for k3, v3 in v2.items():
    #                 print("    " + str(k3))
    #     except:
    #         pass

    logging.info("Excluding redundant peaks.")
    peak_info = exclude_redundant_peaks(peaks, args)

    # print(peak_info)
    logging.info("Computing FDR.")
    fdr_table = compute_fdr(peak_info)
    # print(fdr_table)

    fdr_table.to_csv(stdout, sep=" ")

    if args.bedgraph:
        logging.info("Writing bedgraph to file " + args.bedgraph + ".")
        create_bigwig(init_chip, fdr_table, args)

    logging.info("Done.")
