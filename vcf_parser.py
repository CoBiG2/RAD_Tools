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

# vcf_parser.py is that performs filtering and transformations steps on
# a VCF file that are not possible with vcftools

# Usage: python3 vcf_parser.py -h (will show all available options)

import argparse
import random
from collections import Counter

PARSER = argparse.ArgumentParser(description="Filtering of VCF files")

PARSER.add_argument("-vcf", dest="vcf_infile", help="Provide VCF "
                    " file(s).", required=True)
PARSER.add_argument("--remove-inv", dest="remove_inv", const=True,
                    action="store_const", help="Filters invariable SNP sites"
                    " from the VCF file (These sites may occurr when "
                    "individuals are removed from a VCF file).")
PARSER.add_argument("--remove-singletons", dest="remove_singletons", const=True,
                    action="store_const", help="Filters singletons SNP sites"
                    " from the VCF file.")
PARSER.add_argument("--one-snp", dest="one_snp", const=True,
                    action="store_const", help="Filters the VCF file so that"
                    " only one SNP per locus is retained - The first one")
PARSER.add_argument("--random-snp", dest="rnd_snp", const=True,
                    action="store_const", help="Filters the VCF file so that"
                    " only one random SNP per locus is retained.")
PARSER.add_argument("--center-snp", dest="center_snp", const=True,
                    action="store_const", help="Filters the VCF file so that"
                    " only the SNP closest to the locus center is retained.")

ARG = PARSER.parse_args()

def remove_sites(vcf_file, mode="inv", suffix="_NoInv.vcf"):
    """
    Removes invariable sites from a VCF file. This assumes that the genotype
    columns start at the 10th column until the last column
    :param vcf_file: string, path to vcf file
    :param mode: string, specifies the removal operation. May be one of the
    following:
        .:"inv": Removes invariable sites
        .:"singletons": Removes singleton sites
    :param suffix: string, the suffix of the output file.
    """

    vcf_output = vcf_file.split(".")[0] + suffix

    with open(vcf_file) as vcf_handle, open(vcf_output, "w") as vcf_out:

        for line in vcf_handle:

            if line.startswith("#"):
                vcf_out.write(line)

            elif line.strip() != "":

                # Get genotypes. Remove genotypes with no data.
                genotypes = [x.split(":")[0] for x in line.split()[9:] if x.split(":")[0] != "./."]

                if mode == "inv":
                    # If number of unique genotypes higher than 1, save SNP
                    if len(set(genotypes)) > 1:
                        vcf_out.write(line)

                elif mode == "singletons":
                    gens = Counter(genotypes)
                    # Remove most common genotype
                    del gens[gens.most_common()[0][0]]
                    # Check if remaining genotype count is higher than 1
                    if sum(gens.values()) > 1:
                        vcf_out.write(line)


def filter_one_snp(vcf_file):
    """
    Filters a VCF file so that only one SNP per locus (the first) is retained
    """

    vcf_output = vcf_file.split(".")[0] + "OneSNP.vcf"

    chrom_list = []

    with open(vcf_file) as vcf_handle, open(vcf_output, "w") as vcf_out:

        for line in vcf_handle:

            if line.startswith("#"):
                vcf_out.write(line)

            elif line.strip() != "":

                # Get chrom number
                chrom = line.split()[0]

                if chrom not in chrom_list:
                    vcf_out.write(line)
                    chrom_list.append(chrom)


def filter_random_snp(vcf_file):
    """
    Filters a VCF file so that only one random SNP per locus is retained.
    """

    vcf_output = vcf_file.split(".")[0] + "RandSNP.vcf"

    current_chrom = 0
    loci_snps = []

    with open(vcf_file) as vcf_handle, open(vcf_output, "w") as vcf_out:

        for line in vcf_handle:

            if line.startswith("#"):
                vcf_out.write(line)

            elif line.strip() != "":

                # Get chrom number
                chrom = line.split()[0]

                if chrom != current_chrom and loci_snps != []:
                    choosen = random.choice(loci_snps)
                    vcf_out.write(choosen)
                    loci_snps = [line]
                    current_chrom = chrom

                else:
                    loci_snps += [line]
        vcf_out.write(random.choice(loci_snps))


def filter_center_snp(vcf_file):
    """
    Filters a VCF file so that only one SNP per locus (the first) is retained
    """

    vcf_output = vcf_file.split(".")[0] + "CenterSNP.vcf"

    current_chrom = ""
    line_list = []

    with open(vcf_file) as vcf_handle, open(vcf_output, "w") as vcf_out:

        for line in vcf_handle:

            if line.startswith("#"):
                vcf_out.write(line)

            elif line.strip() != "":

                # Get chrom number
                chrom = line.split()[0]

                # Get SNP position
                pos = int(line.split()[1])

                if chrom != current_chrom and current_chrom != "":
                    closest = min(pos_list, key=lambda x: abs(x - 45))
                    vcf_out.write(line_list[pos_list.index(closest)])
                    pos_list = [pos]
                    line_list = [line]
                    current_chrom = chrom
                elif chrom != current_chrom:
                    pos_list = [pos]
                    line_list = [line]
                    current_chrom = chrom
                else:
                    pos_list += [pos]
                    line_list += [line]

        closest = min(pos_list, key=lambda x: abs(x - 45))
        vcf_out.write(line_list[pos_list.index(closest)])


def main():
    """
    Main function that controls what to do.
    """
    # Args
    vcf_file = ARG.vcf_infile

    if ARG.remove_inv:
        remove_sites(vcf_file)

    if ARG.remove_singletons:
        remove_sites(vcf_file, mode="singletons", suffix="_NoSing.vcf")

    if ARG.one_snp:
        filter_one_snp(vcf_file)

    if ARG.rnd_snp:
        filter_random_snp(vcf_file)

    if ARG.center_snp:
        filter_center_snp(vcf_file)

main()
