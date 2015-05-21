#!/usr/bin/python

# Copyright 2015 Francisco Pina Martins <f.pinamartins@gmail.com>
# This file is part of Loci_counter.
# Loci_counter is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# Loci_counter is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with Loci_counter. If not, see <http://www.gnu.org/licenses/>.

# This script will compare technical replicates on a vcf file and output a
# number of statistics and error rates. The vcf and pair mapping files must
# be provided

import itertools
import numpy as np
import matplotlib.pyplot as plt
from collections import OrderedDict

def parse_pairs(pair_file):
    """
    Parses a txt file in which each line contains three columns with the name
    of the pair and the replicate names for each pair, which must match with the
    names in the VCF file.

    :param pair_file: string, file name containing the sample pairs

    e.g:
    pair_name1  TaxonA  TaxonB
    pair_name2  TaxonC  TaxonD
    (...)
    """

    pairs = {}
    file_handle = open(pair_file)

    for line in file_handle:
        fields = line.split()
        pairs[fields[0]] = fields[1:]

    return pairs


def compare_pairs(vcf_file, pairs):
    """
    Main function comparing the information contained in a VCF file between
    pairs of technical replicates. The relevant statistics will be stored as
    the VCF file is parsed. The statistics being stored for each pair are:

    .: For all VCF loci :.
    ..: total_allele: Total number of alleles (RAD alleles), here defined as
        entire stacks that may contain multiple SNPs, that contain data from
        both replicates
    ..: partial_loci: Total number of loci with date from at least one of the
        replicate pairs
    ..: shared_loci: Number of loci shared among both pair elements

    .: Only for lines with shared loci :.
    ..: snp_number: Total number of SNPs shared by both replicates
    ..: allele_mismatch: Number of mismatches between RAD alleles (consist of
        the entire stack, which may contain more than one SNP) for each locus
    ..: snp_mismatch: Number of mismatches at single sites between replicates

    This will allow the estimation of four error rates:

    - Total locus error rate: The proportion of loci from the full sample that
    are missing for both replicates of a pair
    - Partial locus error rate: The number of loci with data absent from one
    replicate divided by the number of loci with data for at least one replicate
    - Allele error rate: Proportion of allele mismatches divided by the total
    number of alleles for each pair
    - SNP error rate: Proportion of SNP mismatches at individual sites divided
    by the total number of SNPs for that pair

    :param vcf_file: string, path to VCF file
    :param pairs: dictionary, containing pair name as key and a list with the
    names of pair replicates as value

    """

    vcf_handle = open(vcf_file)

    # This will keep track of the previous locus number
    previous_locus = 0

    # This will store the total number of loci
    total_loci = 0

    # This variable will store the column number of each sample
    pos = {}

    # This attribute will store the statistics for each pair. Each pair will
    # contain a dictionary as a value where all relevant statistics will be
    # stored.
    # NOTE: The temp_genotype key will store the genotypes for the current
    # locus being analysed. When the locus changes, this list is analysed and
    # several other keys are updated. Then the temp_genotype list is reseted
    # for the next loci
    pair_statistics = dict((x, {"total_allele": 0,
                                "partial_loci": 0,
                                "shared_loci": 0,
                                "snp_number": 0,
                                "allele_mismatch": 0,
                                "snp_mismatch": 0,
                                "temp_genotype": []}) for x in pairs)

    for line in vcf_handle:
        # Skip header
        while line.startswith("##"):
            line = next(vcf_handle)

        # Get the fields of the samples and associate them to the pairs
        if line.startswith("#CHROM"):
            fields = [x.strip() for x in line.split("\t")]
            for sample in [x.strip() for y in pairs.values() for x in y]:
                pos[sample] = fields.index(sample)

        # Start parsing contents of VCF file for each pair. Most statistics,
        # with the exception of the snp_mismatch and snp_number will only be
        # updated when the locus number changes.
        else:
            fields = [x.strip() for x in line.split("\t")]
            # Get locus number
            current_locus = fields[0]
            # Update total loci number
            total_loci += 1

            # pname: pair name
            # reps: list containing replicate names
            for pname, reps in pairs.items():

                # Get 2 element list with the genotypes for each replicate
                genotype = [fields[pos[reps[0]]], fields[pos[reps[1]]]]

                # Process the first line of the VCF file with content or a line
                # from the same locus as the previous line. Here only
                # the snp_number and snp_mismatch will be updated if both
                # replicates have data
                if previous_locus == 0 or previous_locus == current_locus:

                    if set(genotype) != {"./."}:

                        pair_statistics[pname]["partial_loci"] += 1

                        # If both replicates contain data, update snps stats
                        if "./." not in genotype:
                            pair_statistics[pname]["snp_number"] += 1
                            pair_statistics[pname]["shared_loci"] += 1

                            # If genotypes differ, add to snp_mismatch
                            if len(set(genotype)) != 1:
                                pair_statistics[pname]["snp_mismatch"] += 1

                    # Store genotype in temporary list
                    pair_statistics[pname]["temp_genotype"].append(genotype)

                # This code will only run when the locus number changes. Here,
                # The genotypes contained in the "temp_genotype" will be
                # analyse to update allele and locus statistics. Finally, the
                # list will be reset for the next locus
                else:

                    # Analyses of alleles.
                    glist = pair_statistics[pname]["temp_genotype"]

                    # Allele analyses will only be performed when there is no
                    # missing data in any of the replicates.
                    if "./." not in [x for y in glist for x in y]:
                        pair_statistics[pname]["total_allele"] += 1
                        for gen1, gen2 in glist:
                            if gen1 != gen2:
                                pair_statistics[pname]["allele_mismatch"] += 1
                                break

                    if set(genotype) != {"./."}:

                        pair_statistics[pname]["partial_loci"] += 1

                        # If both replicates contain data, update snps stats
                        if "./." not in genotype:
                            pair_statistics[pname]["snp_number"] += 1
                            pair_statistics[pname]["shared_loci"] += 1

                            # If genotypes differ, add to snp_mismatch
                            if len(set(genotype)) != 1:
                                pair_statistics[pname]["snp_mismatch"] += 1

                    # Store genotype in temporary list
                    pair_statistics[pname]["temp_genotype"] = [genotype]

            previous_locus = current_locus

    return total_loci, pair_statistics


def plot_single_assembly(total_loci, stats, output_name):
    """
    This function will calculate and plot four error rates from a single
    assembly (i.e., vcf file). The error rates are:

    - Total locus error rate: The proportion of loci from the full sample that
    are missing for both replicates of a pair
    - Partial locus error rate: The number of loci with data absent from one
    replicate divided by the number of loci with data for at least one replicate
    - Allele error rate: Proportion of allele mismatches divided by the total
    number of alleles for each pair
    - SNP error rate: Proportion of SNP mismatches at individual sites divided
    by the total number of SNPs for that pair

    :param total_loci: int/float, the total number of loci from a given VCF file
    :param stats: dictionary, containing the statistics for each replicate pair
    obtained from compare_pairs function
    :param output_name: string, name of the output plot

    """

    # Setting plot style
    plt.style.use("ggplot")

    # Set figure with four axes as a 2-d array
    f, ax = plt.subplots(2, 2)

    total_loci = float(total_loci)

    # Data storage for quick ploting
    data = OrderedDict()

    # Total locus error rate
    locus_error_rate = tuple((x, float(total_loci - y["shared_loci"]) /
                             total_loci) for x, y in stats.items())

    # Partial locus error rate
    plocus_error_rate =  tuple((x, float(total_loci - y["partial_loci"]) /
                             total_loci) for x, y in stats.items())

    # Allele error rate
    allele_error_rate = tuple((x, float(y["allele_mismatch"]) /
                        float(y["total_allele"])) for x, y in stats.items())

    # SNP error rate:
    snp_error_rate = tuple((x, float(y["snp_mismatch"]) /
                            float(y["snp_number"])) for x, y in stats.items())

    titles = ["Total locus error", "Partial locus error", "Allele error",
              "SNP error"]

    # Gather plot data
    for i, l in enumerate([locus_error_rate, plocus_error_rate,
                          allele_error_rate, snp_error_rate]):
        data[titles[i]] = [[x[1] for x in l], [x[0] for x in l]]

    # Plot data
    for (k, val), a in zip(data.items(), [x for y in ax for x in y]):
        a.set_title(k)
        a.bar(np.arange(len(val[0])), val[0], align="center")
        a.set_xticks(np.arange(len(val[1])))
        a.set_xticklabels(val[1])

    f.tight_layout()
    plt.savefig(output_name + ".png")


def main():

    import sys

    # Get arguments
    args = sys.argv

    pairs_file = args[1]
    vcf_file = args[2]

    pairs = parse_pairs(pairs_file)
    total_loci, stats = compare_pairs(vcf_file, pairs)

    plot_single_assembly(total_loci, stats, vcf_file.split(".")[0])

main()
