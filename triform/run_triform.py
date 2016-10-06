import logging
from sys import stdout
from os.path import basename, dirname
from subprocess import call

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

from triform.matrix.create_matrix import create_matrix
from triform.create_bed import create_bed


def run_triform(args):

    logging.info("Preprocessing bed files.")
    (treatment, control, treatment_sizes, control_sizes, treatment_iranges,
     control_iranges) = preprocess(args)

    # for bed_file, d in treatment_iranges.items():
    #     print(bed_file, "bed_file")
    #     print(d)
    #     print(type(d))
    #     for chromosome in d:
    #         print(chromosome, "chromosome")
    #         print(d[chromosome])

    #         from rpy2.robjects import r, pandas2ri
    #         r["write.table"](d[chromosome], "iranges.bed", sep=" ")
    # raise

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
    #         print(v2)

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

    logging.info("Computing FDR.")
    fdr_table = compute_fdr(peak_info)

    if fdr_table.empty:
        print("No peaks found.")
    else:
        fdr_table.to_csv(args.outfile, sep=" ")

    logging.info("Done.")

    if args.bed:
        logging.info("Creating bed file: " + args.bed)
        outdir = dirname(args.bed)
        if outdir:
            call("mkdir -p {}".format(outdir), shell=True)
        bed_df = create_bed(treatment_iranges, control_iranges, fdr_table, args)
        bed_df.reset_index().to_csv(args.bed, sep="\t", header=False, index=False)


    if args.bigwig:
        call("mkdir -p {}".format(args.bigwig), shell=True)
        create_bigwig(treatment_iranges, args.bigwig, args)
        create_bigwig(control_iranges, args.bigwig, args)

    if args.matrix:
        matrix = create_matrix(treatment_iranges, control_iranges, fdr_table,
                               args)


def run_triform_no_control(args):

    logging.info("Preprocessing bed files.")
    treatment, treatment_sizes, treatment_iranges = preprocess(args)

    logging.info("Initializing treatment data.")
    init_chip = init_treatment(treatment, args)

    logging.info("Computing region statistics.")
    results = bed_file(init_chip, None, treatment_sizes, None, args)

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
