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

# vcftools_stats.py contains several functions that parse statistics tests
# from vcftools. The idea is to provide any type of vcftools' output with the
# -in option and then specify the type of statistical output it is with a
# second argument

import argparse
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict

# Setting plot style
plt.style.use("ggplot")

parser = argparse.ArgumentParser(description="Parsing statistical output of"
                                 " VCFtools")

parser.add_argument("-in", dest="infile", help="The vcftools output file",
                    required=True)
parser.add_argument("--weir-fst", dest="weir_fst", help="Parse a Weir FST"
                    " output file. Provide a minimum fst threshold to save"
                    " loci.")
parser.add_argument("--singleton", dest="singleton", const=True,
                    action="store_const", help="Parse a singleton/doubleton"
                    " output file")
parser.add_argument("--alternative-snps", dest="altsnp", const=True,
                    action="store_const", help="Parse a VCF file and plot"
                    " the distribution of the number of alternative allele"
                    " SNPs per taxa")

arg = parser.parse_args()


def weir_fst(infile, fst_threshold=None):
    """
    Parses the Weir and Cockerham FST output file
    """

    fst_vals = []

    if fst_threshold:
        fst_fh = open("fst_chosen_loci.txt", "w")

    with open(infile) as fh:

        # Skip header
        next(fh)

        for line in fh:
            if line.strip().split()[2] != "-nan":
                fst = float(line.strip().split()[2])
                fst_vals.append(fst)

                if fst_threshold:
                    if fst >= fst_threshold:
                        fst_fh.write("{}\t{}\n".format(line.strip().split()[0],
                                                       line.strip().split()[1]))

    if fst_threshold:
        fst_fh.close()

    # Creating plots
    # Histogram
    f, ax = plt.subplots()

    plt.hist(fst_vals)

    f.tight_layout()
    plt.savefig("fst_distribution.png")
    plt.close()

    # Plot
    f, ax = plt.subplots()

    plt.plot(fst_vals, "bo")

    f.tight_layout()
    plt.savefig("fst_vals.png")


def singletons(infile):
    """
    Parses the singleton/doubleton output file
    """

    data = {}

    with open(infile) as fh:

        # Skip header
        next(fh)

        for line in fh:
            taxon = line.strip().split()[-1]
            if taxon in data:
                data[taxon] += 1
            else:
                data[taxon] = 1

    # Plot
    f, ax = plt.subplots()

    data = [(x, y) for x, y in data.items()]

    data.sort(key=lambda tup: tup[1])

    plt.bar([x for x in range(len(data))], [x[1] for x in data])
    plt.xticks([x for x in range(len(data))], [x[0] for x in data],
               rotation=45)
    plt.xlim([0, len(data)])

    plt.tight_layout()
    plt.savefig("singleton_distribution.svg")


def alternative_snp_distribution(vcf_file):
    """
    Parses a VCF file and plots the frequency of the number of alternative SNPs
    for each taxon
    """

    vcf_fh = open(vcf_file)

    # Maximum number of total alternative SNPs for plotting
    max_alt = 0

    # Skip header
    for line in vcf_fh:
        if line.startswith("##"):
            pass
        elif line.startswith("#CHROM"):
            # Get taxa names and positions
            taxa_names = line.strip().split()[9:]
            storage = dict((x, defaultdict(int)) for x in taxa_names)
        else:
            fields = line.strip().split()
            # Get genotypes
            genotypes = [x.split(":")[0] for x in fields[9:]]

            # Get total number of alternative SNPs
            alt_num = 0
            for i in genotypes:
                if "1" in i:
                    alt_num += 1

            # Update max alternative alleles
            if alt_num > max_alt:
                max_alt = alt_num

            # Get number of alternative SNPs for each taxa
            for p, i in enumerate(genotypes):
                if "1" in i:
                    # Add 1 to the taxa counter for the number of alt_num
                    storage[taxa_names[p]][alt_num] += 1

    # Prepare data for plotting
    plot_data = []
    for i in range(1, max_alt + 1):
        temp_list = []
        for taxon in taxa_names:
            temp_list.append(storage[taxon][i])
        plot_data.append(temp_list)

    # Creating plot
    fig, ax = plt.subplots()

    colors = plt.cm.jet(np.linspace(0, 2, len(taxa_names)))

    ind = np.arange(len(taxa_names))
    y_offset = np.array([0.0] * len(taxa_names))

    plots = []
    w = 0.80
    for i in range(len(plot_data)):
        p = plt.bar(ind, plot_data[i], w, bottom=y_offset, color=colors[i])
        y_offset = y_offset + plot_data[i]
        plots.append(p)

    plt.xticks(ind + w/2., taxa_names, rotation=45, ha="right")

    plt.legend(plots, range(1, max_alt + 1))

    fig.tight_layout()
    plt.savefig("Alternative_SNPs_distribution.png")

def main():

    # Arguments
    infile = arg.infile

    if arg.weir_fst:
        weir_fst(infile, float(arg.weir_fst))

    if arg.singleton:
        singletons(infile)

    if arg.altsnp:
        alternative_snp_distribution(infile)

main()
