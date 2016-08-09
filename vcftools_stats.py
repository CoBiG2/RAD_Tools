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
from collections import defaultdict, Counter, OrderedDict
import itertools
import statsmodels.sandbox.stats.multicomp as multi_correction


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
parser.add_argument("--shared-snps", dest="shared", help="Plots putatively"
                    " introgressed loci for each taxa provided in the first"
                    " file against the taxa in the second file", nargs=2)
parser.add_argument("--fst-vals", dest="fst_vals", help="Provide FST output"
                    " output file from VCFtools. Required or --shared-snps"
                    " option.")
parser.add_argument("--remove-int", dest="remove_int", const=True,
                    action="store_const", help="Use this option to create a new"
                    " vcf file without the putatively introgressed loci")
parser.add_argument("--filter-fst", dest="filter_fst", nargs="+", help="Filter"
                    " a vcf file provided in '-in' according to the fst values"
                    " provided in '--fst-vals' tha fall between the values"
                    " provided in this option (e.g. '--filter-fst 0.8 1')."
                    " Alternatively, only one value can be provided and Only"
                    " the SNPs with that exact FST value are saved.")
parser.add_argument("--hwe", dest="hwe", help="Parse a hardy-weinberg .hardy"
                    " output file. Provide a threshold for the q-value")
parser.add_argument("--filter-hwe", dest="filter_hwe", help="Filters a vcf file"
                    " according to the hardy-weinberg .hardy file. The p-values"
                    " of the hwe test are corrected with a FDR approach.")
parser.add_argument("--count-loci", dest="count_loci", const=True,
                    action="store_const", help="Prints the number of loci")


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
    plt.savefig("fst_distribution.pdf")
    plt.close()

    # Plot
    f, ax = plt.subplots()

    plt.plot(fst_vals, "bo")

    f.tight_layout()
    plt.savefig("fst_vals.pdf")


def parse_hwe(f, alpha, vcf_file):
    """
    Parses a hardy-weinberg output file, corrects p-values according to a FDR
    and generates several plots to visualize the hwe results
    """

    vcf_outfile = vcf_file.split(".")[0] + "_filtered.vcf"

    snp_pos = []
    pvals = []
    het_deficit = []
    het_excess = []

    with open(f) as fh:
        #Skip header
        next(fh)

        for line in fh:
            fields = line.strip().split()
            snp_pos.append((fields[0], fields[1]))
            pvals.append(float(fields[5]))
            het_deficit.append(float(fields[6]))
            het_excess.append(float(fields[7]))

    fdr_bool_list, fdr_pvalue_list, alpha_S, alpha_B = \
        multi_correction.multipletests(pvals, alpha=float(alpha),
                                       method="fdr_bh")

    snp_pvals = OrderedDict()
    for pos, pval in zip(snp_pos, fdr_pvalue_list):
        snp_pvals["-".join(pos)] = pval

    with open(vcf_file) as vcf_fh, open(vcf_outfile, "w") as ofh:
        for line in vcf_file:
            if line.startswith("#"):
                ofh.write(line)
            elif line.strip() != "":
                fields = line.split()
                # Check pval for locus
                pos = "-".join(fields[0], fields[1])
                if snp_pvals[pos] <= 0.05:
                    ofh.write(line)



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


def parse_taxa_file(f):
    """
    Returns a list of taxa from a taxa file
    """

    fh = open(f)

    return [x.strip() for x in fh.readlines() if x.strip() != ""]


def parse_fst(fst_file, fst_range=None):
    """
    Parses an Fst file from vcftools. Returns a dictionary with
    the chromosome and postition as key and fst value as value. If he fst_range
    argument is provided, it will only store SNPs within the fst_range
    """

    fh = open(fst_file)
    fst_storage = {}

    # Skip header
    next(fh)

    for line in fh:

        if line.strip() != "":
            fields = line.strip().split()
            if fields[2] != "-nan":
                fst = float(fields[2])
                if fst_range:
                    if len(fst_range) == 1 and fst == float(fst_range[0]):
                        fst_storage["{}_{}".format(fields[0], fields[1])] = fst
                    else:
                        if float(fst_range[0]) <= fst <= float(fst_range[1]):
                            fst_storage["{}_{}".format(fields[0], fields[1])] = fst
                else:
                    fst_storage["{}_{}".format(fields[0], fields[1])] = fst
            else:
                # Only save SNPs with -nan when the fst_range is not specified
                if not fst_range:
                    fst_storage["{}_{}".format(fields[0], fields[1])] = 0

    return fst_storage


def filter_fst(vcf_file, fst_storage):
    """
    Filters a vcf_file so that it includes only the SNPs from the fst_storage
    """

    vcf_fh = open(vcf_file)
    out_vcf = open(vcf_file.split(".")[0] + "_filtered.vcf", "w")

    for line in vcf_fh:
        if line.startswith("#"):
            out_vcf.write(line)
        elif line.strip() != "":
            fields = line.split()
            # Get chrom and position
            coord = "{}_{}".format(fields[0], fields[1])

            if coord in fst_storage:
                out_vcf.write(line)

    vcf_fh.close()
    out_vcf.close()


def introgressed(vcf_file, p1, p2, fst_storage):
    """
    :param vcf_file: path to vcf file
    :param p1: list, taxa to count shared polymorphisms
    :param p2: list, reference taxa
    """

    vcf_fh = open(vcf_file)

    if arg.remove_int:
        filtered_vcf = open(vcf_file.split(".")[0] + "_filtered.vcf", "w")

    het_taxa_storage = OrderedDict((x, 0) for x in p1)
    hom_taxa_storage = OrderedDict((x, 0) for x in p1)
    prop_taxa_storage = OrderedDict((x, 0) for x in p1)

    chrom_vals = OrderedDict((x, {}) for x in p1)

    introgressed_chroms = OrderedDict((x, []) for x in p1)

    # Counter of full diagnostic SNPs
    diagnostic = 0

    # Variable that will determine whether the current SNP is to be filtered
    # (False) or not (True)
    flag = True

    for line in vcf_fh:
        if line.startswith("##"):
            if arg.remove_int:
                filtered_vcf.write(line)
            else:
                pass
        elif line.startswith("#CHROM"):
            if arg.remove_int:
                filtered_vcf.write(line)

            taxa_list = line.strip().split()

        elif line.strip() != "":
            fields = line.strip().split()

            # Get locus position
            loc_pos = "{}_{}".format(fields[0], fields[1])

            # Evaluate fst value
            if fst_storage[loc_pos] > 0.8:
                diagnostic += 1
                # Get genotypes for p2
                p2_geno = [fields[taxa_list.index(x)].split(":")[0] for x in p2]
                # Get most common allele from p2
                try:
                    p2_al = Counter("".join(p2_geno).replace("|","").replace(".","").replace("/","")).most_common(1)[0][0]
                except IndexError:
                    continue
                # Get shared alleles for each taxa in p1
                for taxon in p1:
                    if fields[0] not in chrom_vals[taxon]:
                        chrom_vals[taxon][fields[0]] = []
                    # Get genotype for taxon
                    gen = fields[taxa_list.index(taxon)].split(":")[0]
                    al_count = gen.count(p2_al)
                    if al_count == 1:
                        # For shared Heterozygous SNPs set flag so that they
                        # are filtered from the VCF
                        flag = False
                        het_taxa_storage[taxon] += 1
                        chrom_vals[taxon][fields[0]].append(1)
                        introgressed_chroms[taxon].append(fields[0])
                    elif al_count == 2:
                        hom_taxa_storage[taxon] += 1
                        chrom_vals[taxon][fields[0]].append(1)
                        introgressed_chroms[taxon].append(fields[0])
                    else:
                        chrom_vals[taxon][fields[0]].append(0)

            if flag and arg.remove_int:
                filtered_vcf.write(line)

            # Reset flag value for next iteration
            flag = True

    # This piece of code is highly ad-hoc to the three Hv taxa with putatively
    # introgressed signals.
    try:
        from matplotlib_venn import venn3

        venn3([set(introgressed_chroms["999pu1_3"]),
               set(introgressed_chroms["3624_1"]),
               set(introgressed_chroms["2377pu1_2"])],
              ["999pu1_3", "3624_1", "2377pu1_2"])
        plt.savefig("introgressed_venn.pdf")

        shared_tb = open("shared_int.csv", "w")
        for comb in itertools.combinations(["999pu1_3", "3624_1", "2377pu1_2"], 2):
            s1 = set(introgressed_chroms[comb[0]])
            s2 = set(introgressed_chroms[comb[1]])
            shared = s1.intersection(s2)
            unique1 = s1.difference(s2)
            unique2 = s2.difference(s1)

            shared_tb.write("Combination {}-{}\n".format(comb[0], comb[1]))
            shared_tb.write("Shared: {}\nUnique {}: {}\nUnique {}: {}\n".format(
                len(shared), comb[0], len(unique1), comb[1], len(unique2)
            ))

        shared_tb.close()

    except:
        pass

    # Writes to a tabular file the number and proportions of SNPs in the same
    # chromosome that match the signal to each other
    chrom_tb = open("chromosome_match.csv", "w")
    chrom_tb.write("Taxon; Total; Match;%; Mismatch;%;\n")

    for taxon in chrom_vals:
        total = 0
        match = 0
        mismatch = 0
        for val in chrom_vals[taxon].values():
            if len(val) > 1:
                total += 1
                if len(set(val)) == 1:
                    match += 1
                else:
                    mismatch += 1
        chrom_tb.write("{};{};{};{};{};{}\n".format(
            taxon, total, match, (float(match)/float(total))*100,
            mismatch, (float(mismatch)/float(total))*100))

    chrom_tb.close()

    for t, het, hom in zip(p1, het_taxa_storage.values(), hom_taxa_storage.values()):
        prop_taxa_storage[t] = ((het + hom) / diagnostic) * 100

    # Generate table
    output = open("Shared_alleles.csv", "w")
    output.write("Taxon; Shared (Heterozygous); Shared (Homozygous); %\n")
    for t in p1:
        output.write("{}; {}; {}; {}\n".format(t, het_taxa_storage[t],
                                               hom_taxa_storage[t],
                                               prop_taxa_storage[t]))
    output.close()

    # Plot bar plot with shared allele count for each taxon
    fig, ax1 = plt.subplots()

    ind = np.arange(len(p1))
    w = 0.8

    # Heterozygous data
    het_data = [x for x in het_taxa_storage.values()]
    # Homozygous data
    hom_data = [x for x in hom_taxa_storage.values()]

    bar1 = ax1.bar(ind, het_data, w, color="blue")
    bar2 = ax1.bar(ind, hom_data, w, color="green", bottom=het_data)

    plt.xticks(ind + w / 2., p1, rotation=45, ha="right")
    ax1.set_xlabel("Taxa")
    ax1.set_ylabel("Frequency")

    ax2 = ax1.twinx()
    # Percentage data
    perc_data = [x for x in prop_taxa_storage.values()]
    ax2.plot([x + w/2. for x in ind], perc_data, marker="+", ls="-")
    ax2.set_ylabel("Percentage")

    fig.tight_layout()
    plt.savefig("Shared_alleles.pdf")

def print_loci(vcf_file):
    """
    Prints the number of loci in the VCF file
    """

    vcf_fh = open(vcf_file)
    loci = []

    for line in vcf_fh:
        if line.startswith("#"):
            pass
        elif line.strip() != "":
            # Get locus number
            l = line.split()[0]

            if l not in loci:
                loci.append(l)

    print(len(loci))


def main():

    # Arguments
    infile = arg.infile

    if arg.weir_fst:
        weir_fst(infile, float(arg.weir_fst))

    if arg.singleton:
        singletons(infile)

    if arg.altsnp:
        alternative_snp_distribution(infile)

    if arg.shared:
        p1_file, p2_file = arg.shared[0], arg.shared[1]
        fst_file = arg.fst_vals

        # Parse taxa files
        p1 = parse_taxa_file(p1_file)
        p2 = parse_taxa_file(p2_file)
        fst = parse_fst(fst_file)

        # Get introgressed loci
        introgressed(infile, p1, p2, fst)

    if arg.filter_fst:
        fst_storage = parse_fst(arg.fst_vals, arg.filter_fst)
        filter_fst(arg.infile, fst_storage)

    if arg.hwe:
        parse_hwe(infile, arg.hwe)

    if arg.count_loci:
        print_loci(infile)

main()
