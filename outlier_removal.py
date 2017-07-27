#!/usr/bin/env python3

import argparse
import sys


def argument_parser(args):
    """
    Parses the list of arguments as implemented in argparse.
    """
    parser = argparse.ArgumentParser(
        description="A script to filter .VCF and .loci files based on loci list"
                    " files.",
        formatter_class=argparse.RawTextHelpFormatter)

    # Create subparsers for each main operation
    subparsers = parser.add_subparsers(
        help="Select which command you wish to execute.",
        dest="main_op")

    neutral_vcf = subparsers.add_parser("neutral", help="Outputs a VCF with "
                                        "only neutral loci.")
    selection_vcf = subparsers.add_parser("selection", help="Outputs a VCF "
                                          "with only loci under selection.")
    selection_fasta = subparsers.add_parser("fasta", help="Outputs a FASTA "
                                            "file with the sequences of "
                                            "loci under selection.")

    # ###################### Neutral VCF ARGUMENTS #############################
    # Group definition
    io_opts = neutral_vcf.add_argument_group("Input/Output options")

    io_opts.add_argument("-v", dest="vcf_filename", type=str, required=True,
                         metavar="vcf",
                         help="Location of the VCF file to filter.")
    io_opts.add_argument("-o", dest="outlier_loci", type=str, nargs="+",
                         metavar="list of files", default=None,
                         help="Path to outlier loci lists. Can provide more "
                         "than one file path.")
    io_opts.add_argument("-a", dest="assoc_loci", type=str, default=None,
                         metavar="path to associations summary",
                         help="Path to assocaitions summary file.")

    # ###################### Selection VCF ARGUMENTS ###########################
    # Group definition
    io_opts = selection_vcf.add_argument_group("Input/Output options")

    io_opts.add_argument("-v", dest="vcf_filename", type=str, required=True,
                         metavar="vcf",
                         help="Location of the VCF file to filter.")
    io_opts.add_argument("-o", dest="outlier_loci", type=str, nargs="+",
                         metavar="list of files", default=None,
                         help="Path to outlier loci lists. Can provide more "
                         " than one file path.")
    io_opts.add_argument("-a", dest="assoc_loci", type=str, default=None,
                         metavar="path to associations summary",
                         help="Path to assocaitions summary file.")

    # ################### Selection FASTA ARGUMENTS ############################
    # Group definition
    io_opts = selection_fasta.add_argument_group("Input/Output options")

    io_opts.add_argument("-v", dest="vcf_filename", type=str, required=True,
                         metavar="vcf",
                         help="Location of the VCF file to filter.")
    io_opts.add_argument("-o", dest="outlier_loci", type=str, nargs="+",
                         metavar="list of files", default=None,
                         help="Path to outlier loci lists. Can provide more "
                         " than one file path.")
    io_opts.add_argument("-a", dest="assoc_loci", type=str, default=None,
                         metavar="path to associations summary",
                         help="Path to assocaitions summary file.")

    # ################### END OF SPECIFIC CODE ###############################
    arguments = parser.parse_args(args)

    return arguments


def parse_outliers(outlier_filename):
    """
    Parses an outlier summary file and returns a list with the outlier SNPs.
    """
    infile = open(outlier_filename, 'r')
    loci = set()
    for lines in infile:
        lines = lines.strip()
        if lines == "":
            break
        else:
            lines = lines.strip().split()[1:]
            [loci.add(int(x)) for x in lines]

    infile.close()
    return loci


def parse_associations(association_summary_filename):
    """
    Parses an association summary file and return a list with the reported
    loci.
    """
    infile = open(association_summary_filename, 'r')
    loci = set()
    for lines in infile:
        lines = lines.split()
        loci.add(int(lines[1]))

    infile.close()
    return loci


def parse_vcf(vcf_filename, outliers, goal):
    """
    Parses a vcf file and prints only the lines not in the outliers set.
    """
    infile = open(vcf_filename, 'r')
    counter = 1
    for lines in infile:
        lines = lines.strip()
        if lines.startswith("#") and goal != "fasta":
            print(lines)
        else:
            if counter not in outliers and goal == "neutral":
                print(lines)
            elif counter in outliers and goal == "selection":
                print(lines)
            elif counter in outliers and goal == "fasta":
                print("{0}\t{1}".format(lines.split()[0], counter))
            if not lines.startswith("#"):
                counter += 1

    infile.close()


def runner(arg):
    """
    Manages what gets run and runs it.
    """
    loci_set = None
    if arg.outlier_loci is not None:
        for loci_list in arg.outlier_loci:
            outliers = parse_outliers(loci_list)
            if loci_set is None:
                loci_set = outliers
            else:
                loci_set = loci_set.intersection(outliers)

    if arg.assoc_loci is not None:
        assocs = parse_associations(arg.assoc_loci)
        if loci_set is None:
            loci_set = assocs
        else:
            loci_set = loci_set.union(assocs)

    parse_vcf(arg.vcf_filename, loci_set, arg.main_op)


def main():
    """
    Main function. Calls the argparser and the runner function.
    """
    # Make sure we provide an help message instead of an error
    if len(sys.argv) == 1:
        sys.argv += ["-h"]
    arg = argument_parser(sys.argv[1:])

    runner(arg)


if __name__ == "__main__":
    main()
