import logging
from sys import stdout
from os.path import basename

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

from rpy2.robjects import r

from triform.matrix.find_read_midpoints import find_read_midpoints


def run_triform(args):

    logging.info("Preprocessing bed files.")
    (treatment, control, treatment_sizes, control_sizes, treatment_iranges,
     control_iranges) = preprocess(args)

    if args.matrix:
        logging.info("Creating treatment matrix.")
        treatment_matrixes = find_read_midpoints(treatment_iranges, args)
        # logging.info("Creating control matrix.")
        # control_matrixes = find_read_midpoints(control_iranges, args)
        logging.info("Merging matrixes")

        # assert set(treatment_matrixes) == set(control_matrixes), "Chromosomes in ChIP and input differ"

        # for chromosome in set(treatment_matrixes).intersection(control_matrixes):
        chromosome = "chrY"
        all_granges = iter(
            treatment_matrixes[chromosome])  # + control_matrixes[chromosome])
        name, u = next(all_granges)

        for (fname, gr) in all_granges:
            _merge = r("""function(u, gr){{

             elementMetadata(u[overlapsAny(u, gr1)])${}= = elementMetadata(gr1)[,1]

             }}""".format(basename(fname)))

            u = _merge(u, gr)

        print(u)

        # need for loop to insert names of files (instead of gr1, gr2)
        # elementMetadata(u[overlapsAny(u, gr1)])$gr1 = elementMetadata(gr1)[,1]

        # elementMetadata(u[overlapsAny(u, gr2)])$gr2 = elementMetadata(gr2)[,1]

        # GRanges object with 4 ranges and 2 metadata columns:
        #       seqnames    ranges strand |       gr1       gr2
        #          <Rle> <IRanges>  <Rle> | <numeric> <numeric>
        #   [1]     chr1  [ 1,  9]      * |         6         0
        #   [2]     chr1  [11, 19]      * |         7         1
        #   [3]     chr1  [21, 29]      * |         8         0
        #   [4]     chr1  [31, 39]      * |         0         5
        # m = _create_matrix(list(treatment_iranges.values())[0]["chrY"], "chrY")

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

    # print(peak_info)
    logging.info("Computing FDR.")
    fdr_table = compute_fdr(peak_info)
    # print(fdr_table)

    fdr_table.to_csv(stdout, sep=" ")

    if args.bedgraph:
        logging.info("Writing bedgraph to file " + args.bedgraph + ".")
        create_bigwig(init_chip, fdr_table, args)

    logging.info("Done.")

    # if args.matrix:


def run_triform_no_control(args):

    logging.info("Preprocessing bed files.")
    treatment, treatment_sizes, treatment_iranges = preprocess(args)

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
