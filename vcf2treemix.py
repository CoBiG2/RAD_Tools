#!/usr/bin/python

# vcf2treemix.py
# Converts a vcf file into TreeMix input

import argparse
from collections import OrderedDict

parser = argparse.ArgumentParser(description="Parsing statistical output of"
                                 " VCFtools")

parser.add_argument("-vcf", dest="vcf_file", help="Path to VCF file",
                    required=True)
parser.add_argument("-pop", dest="pop_file", help="Path to pop file",
                    required=True)

arg = parser.parse_args()


def get_pops(pop_file):
    """
    Returns a dictionary with pop identifier as key and taxa as a list of
    strings. In the pop file, each populations should be in one line, starting
    withe pop name, a colon and the corresponding taxa separated by whitespace.
    E.g.:
    pop1: taxon1 taxon2 taxon3
    """

    pops = OrderedDict()

    with open(pop_file) as fh:

        for line in fh:
            fields = line.strip().split(":")
            pops[fields[0]] = fields[1].split()

    return pops


def vcf2treemix(vcf_file, pop_obj):
    """
    Converts a vcf file into treemix format.
    """

    vcf_fh = open(vcf_file)
    output_name = vcf_file.strip(".vcf") + ".tmix"
    output_fh = open(output_name, "w")

    # Write header for tmix file
    output_fh.write("{}\n".format(" ".join([x for x in pop_obj.keys()])))

    for line in vcf_fh:

        # Skip header
        if line.startswith("##"):
            pass

        # Get taxon positions
        elif line.startswith("#CHROM"):
            taxa_pos = line.strip().split()

        # Ignore empty lines
        elif line.strip() != "":

            fields = line.strip().split()

            # Ignore loci with more than two alleles
            if len(fields[4]) > 1:
                continue

            # Get allele counts for each populations
            temp_pop = OrderedDict((x, [0,0]) for x in pop_obj.keys())
            for pop, taxa in pop_obj.items():
                for taxon in taxa:
                    # Get taxon genotype
                    gen = fields[taxa_pos.index(taxon)]
                    # Skip if gen is missing data
                    if gen == "./.":
                        continue

                    temp_pop[pop][0] += gen.count("0")
                    temp_pop[pop][1] += gen.count("1")

            # Write current locus to file
            output_fh.write("{}\n".format(" ".join([str(x[0]) +  "," + str(x[1]) for x in temp_pop.values()])))

    vcf_fh.close()
    output_fh.close()


def main():
    # Args
    vcf_file = arg.vcf_file
    pop_file = arg.pop_file

    pop_obj = get_pops(pop_file)
    vcf2treemix(vcf_file, pop_obj)


main()
