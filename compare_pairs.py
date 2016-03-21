#!/usr/bin/python

# Copyright 2015 Diogo N. Silva <o.diogosilva@gmail.com>
# compare_pairs.py is free software: you can redistribute it and/or modify
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

import argparse
import itertools
import numpy as np
import matplotlib.pyplot as plt
from collections import OrderedDict

parser = argparse.ArgumentParser(description="Comparison of RAD assemblies "
                                             "using technical replicates")

parser.add_argument("-vcf", dest="vcf_infiles", nargs="+", help="Provide VCF "
                    " file(s).", required=True)
parser.add_argument("-pairs", dest="pairs_file", help="Provide the tab "
                    " separated file containing the techical pairs.")
parser.add_argument("-s", dest="single_assembly", action="store_const",
                    const=True, help="Use this flag to analyse each vcf file "
                    "independently.")
parser.add_argument("-o", dest="output_plot", help="Name of the output plot "
                    " file. Only has an effect when the -s option is not "
                    " specified")
parser.add_argument("-c", dest="convert_vcf", action="store_const",
                    const=True, help="Use this option to create"
                    " a filtered .vcf file without erroneous SNPs")
parser.add_argument("-f", dest="filter_vcf", help="Use this option to filter"
                    " a vcf file according to a text file containing the "
                    "chromossome positions to be excluded.")

arg = parser.parse_args()


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


def compare_pairs(vcf_file, pairs, filt_vcf=False):
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
    :param filt_vcf: boolean, whether a new filtered vcf file will be created
    without erroneous SNPs (True) or not (False)

    """

    vcf_handle = open(vcf_file)

    # Stores the locus number of the error containing loci
    bad_loci = []

    # Stores the position of bad snps for plotting
    bad_snps = []

    # Stores the number of bad SNPs for each bad_loci for plotting
    bad_loci_distribution = {}

    # Create filtered vcf file handle if specified
    if filt_vcf:
        vcf_name = vcf_file.split(".")[0] + "_filtered.vcf"
        filtered_vcf = open(vcf_name, "w")

    # This will keep track of the previous locus number
    previous_locus = 0

    # This will store the total number of loci
    total_loci = 0

    # Total number of chromossome
    total_chromossomes = 1
    record_chromossome = True

    # This variable will store the column number of each sample
    pos = {}

    # This variable is used when filtering a vcf file. When a vcf line should be
    # filtered, this variable is set to False
    filter_locus = True

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
            if filt_vcf:
                filtered_vcf.write(line)
                line = next(vcf_handle)
            else:
                line = next(vcf_handle)

        # Get the fields of the samples and associate them to the pairs
        if line.startswith("#CHROM"):
            if filt_vcf:
                filtered_vcf.write(line)
            fields = [x.strip() for x in line.split("\t")]
            for sample in [x.strip() for y in pairs.values() for x in y]:
                pos[sample] = fields.index(sample)

        # Start parsing contents of VCF file for each pair. Most statistics,
        # with the exception of the snp_mismatch and snp_number will only be
        # updated when the locus number changes.
        elif line.strip() != "":
            fields = [x.strip() for x in line.split("\t")]
            # Get locus number
            if fields[0] != "un":
                current_locus = fields[0]
                current_pos = fields[1]
            else:
                current_locus = fields[2]
            # Update total loci number
            total_loci += 1

            # pname: pair name
            # reps: list containing replicate names
            for pname, reps in pairs.items():

                # Get 2 element list with the genotypes for each replicate
                genotype = [fields[pos[reps[0]]].split(":")[0],
                            fields[pos[reps[1]]].split(":")[0]]

                if set(genotype) != {"./."}:

                    pair_statistics[pname]["partial_loci"] += 1

                    # If both replicates contain data, update snps stats
                    if "./." not in genotype:
                        pair_statistics[pname]["snp_number"] += 1
                        pair_statistics[pname]["shared_loci"] += 1

                        # If genotypes differ, add to snp_mismatch
                        if len(set(genotype)) != 1:
                            pair_statistics[pname]["snp_mismatch"] += 1
                            bad_loci.append((current_locus, current_pos))
                            bad_snps.append(int(fields[1]))
                            if current_locus in bad_loci_distribution:
                                bad_loci_distribution[current_locus] += 1
                            else:
                                bad_loci_distribution[current_locus] = 1
                            filter_locus = False

                # Process the first line of the VCF file with content or a line
                # from the same locus as the previous line. Here only
                # the snp_number and snp_mismatch will be updated if both
                # replicates have data
                if previous_locus == 0 or previous_locus == current_locus:

                    # Store genotype in temporary list
                    pair_statistics[pname]["temp_genotype"].append(genotype)

                # This code will only run when the locus number changes. Here,
                # The genotypes contained in the "temp_genotype" will be
                # analyse to update allele and locus statistics. Finally, the
                # list will be reset for the next locus
                else:

                    if record_chromossome:
                        total_chromossomes += 1
                        record_chromossome = False

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

                    # Store genotype in temporary list
                    pair_statistics[pname]["temp_genotype"] = [genotype]

            if filt_vcf and filter_locus:
                filtered_vcf.write(line)

            filter_locus = True
            record_chromossome = True

            previous_locus = current_locus

    # Write bad loci list to a file
    bad_loci = sorted(list(set(bad_loci)))
    with open("bad_loci.text", "w") as bad_fh:
        for l in bad_loci:
            bad_fh.write("{}\t{}\n".format(l[0], l[1]))

    # Plot error position and quantity distribution
    try:
        plt.hist(bad_snps)
        plt.savefig("bad_snp_positions.png")
        plt.close()
    except ValueError:
        print("No SNP mismatch found")

    try:
        plt.hist(list(bad_loci_distribution.values()))
        plt.savefig("bad_loci_distribution.png")
        plt.close()
    except ValueError:
        print("No locus with SNP mismatch found")

    print(total_chromossomes)

    return total_loci, pair_statistics, total_loci


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


def plot_multiple_assemblies(multi_total_loci, multi_stats, plot_name,
                             create_table=True, total_snps=None):
    """
    This function plots the four error rates for each assembly stats provided
    by the multi_stats argument. It calculates the mean error for each error
    rate and each assembly and displays a whisker plot.

    :param multi_total_loci: list, each element is the total number of loci
    for each assembly in the same position in multi_stats
    :param multi_stats: OrderedDict object, each key is the name of the assembly
    and the corresponding values are OrderedDict objects generated with the
    compare_pairs function.
    :param plot_name: string, name of the output
    :param create_table: Boolean, Whether (True) the results are saved in a
    .csv file in table format or (False) not.
    """

    # Setting plot style
    plt.style.use("ggplot")

    # Set figure with four axes as a 2-d array
    f, ax = plt.subplots(2, 2)

    # Data storage for quick ploting
    data = OrderedDict([("Total locus error", []),
                        ("Partial locus error", []),
                        ("Allele error", []),
                        ("SNP error", [])])
    assemblies = []

    if create_table:
        assembly_means = []
        table_handle = open("{}.csv".format(plot_name), "w")

    for total_loci, (a_name, stats) in zip(multi_total_loci,
                                           multi_stats.items()):
        total_loci = float(total_loci)
        assemblies.append(a_name)

        # Total locus error rate
        locus_error_rate = [float(total_loci - x["shared_loci"]) / total_loci
                            for x in stats.values()]
        data["Total locus error"].append(locus_error_rate)

        # Partial locus error rate
        plocus_error_rate = [float(total_loci - x["partial_loci"]) /
                             total_loci for x in stats.values()]
        data["Partial locus error"].append(plocus_error_rate)

        # Allele error rate
        allele_error_rate = [float(x["allele_mismatch"]) /
                             float(x["total_allele"]) for x in stats.values()]
        data["Allele error"].append(allele_error_rate)

        # SNP error rate
        snp_error_rate = [float(x["snp_mismatch"]) / float(x["snp_number"])
                          for x in stats.values()]
        data["SNP error"].append(snp_error_rate)

    for (k, val), a in zip(data.items(), [x for y in ax for x in y]):
        # Creating plot
        if total_snps:
            a2 = a.twinx()
            a2.plot([x + 1 for x in range(len(total_snps))], total_snps,
                    lw=1.5, color="green", alpha=.3)
            a2.set_ylabel("SNPs")
            for t1 in a2.get_yticklabels():
                t1.set_color("green")

        a.set_title(k)
        a.boxplot(val)
        a.set_xticklabels(assemblies, rotation=45, ha="right")
        a.set_ylabel("Proportion")
        a.set_xlim([0.5, len(assemblies) + .5])

        # Create table
        if create_table:
            # Get mean value for each error rate
            er_mean = np.mean([x for y in val for x in y])
            table_handle.write("{}:; {}\n".format(k, er_mean))

            assembly_means.append([np.mean(x) for x in val])

    if create_table:
        table_handle.write("\nAssembly; Total locus error; Partial locus error"
                           "; Allele error; SNP error\n")
        for aname, vals in zip(assemblies, [x for x in zip(*assembly_means)]):
            table_handle.write("{}; {}\n".format(aname,
                                            ";".join([str(x) for x in vals])))

        table_handle.close()

    f.tight_layout()
    plt.savefig(plot_name + ".png")


def filter_vcf_chromossmes(vcf_file, bad_loci_file):
    """
    Filters the chromossome positions of a vcf_file that are contained in the
    bad_loci_file.

    :param vcf_file: string, path to the vcf file
    :param bad_loci_file: string, path to the text file containing the
    chromossomes to be excluded
    """

    with open(vcf_file) as vcf_fh, open(bad_loci_file) as bad_fh, \
        open("{}_filtered.vcf".format(vcf_file.split(".")[0]), "w") as filt_fh:

        bad_loci = bad_fh.read().split("\n")

        bad_loci = [x.split() for x in bad_loci]

        print(bad_loci)

        for line in vcf_fh:
            if line.startswith("#"):
                filt_fh.write(line)
            elif line.strip() != "":
                loci_number = line.split()[0]
                pos_number = line.split()[1]
                if [loci_number, pos_number] not in bad_loci:
                    filt_fh.write(line)

def main():

    import sys

    # Parse arguments
    pairs_file = arg.pairs_file
    vcf_files = arg.vcf_infiles
    plot_file = arg.output_plot
    convert_vcf = arg.convert_vcf
    filter_vcf = arg.filter_vcf

    # Parse pairs file
    pairs = parse_pairs(pairs_file)

    if arg.filter_vcf:
        filter_vcf_chromossmes(vcf_files[0], filter_vcf)

    elif arg.single_assembly:
        for vfile in vcf_files:
            total_loci, stats, snps = compare_pairs(vfile, pairs, convert_vcf)
            plot_single_assembly(total_loci, stats, vfile.split(".")[0])

    else:
        multi_stats = OrderedDict()
        multi_loci = []
        total_snps = []
        for vfile in vcf_files:
            total_loci, stats, snps = compare_pairs(vfile, pairs)
            multi_stats[vfile] = stats
            multi_loci.append(total_loci)
            total_snps.append(snps)
        plot_multiple_assemblies(multi_loci, multi_stats, plot_file,
                                 total_snps=total_snps)

main()
