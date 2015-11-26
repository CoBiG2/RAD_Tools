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

import argparse

parser = argparse.ArgumentParser(description="Filtering of VCF files")

parser.add_argument("-vcf", dest="vcf_infile", help="Provide VCF "
                    " file(s).", required=True)
parser.add_argument("--remove-inv", dest="remove_inv", const=True,
                    action="store_const", help="Filters invariable SNP sites"
                    " from the VCF file (These sites may occurr when "
                    "individuals are removed from a VCF file).")

arg = parser.parse_args()

def remove_invariable(vcf_file):
    """
    Removes invariable sites from a VCF file. This assumes that the genotype
    columns start at the 10th column until the last column
    """

    vcf_output = vcf_file.split(".")[0] + "_NoInv.vcf"

    with open(vcf_file) as vcf_handle, open(vcf_output, "w") as vcf_out:

        for line in vcf_handle:

            if line.startswith("#"):
                vcf_out.write(line)

            elif line.strip() != "":

                # Get genotypes. Remove genotypes with no data.
                genotypes = [x.split(":")[0] for x in line.split()[9:] if x.split(":")[0] != "./."]

                # If number of unique genotypes higher than 1, save SNP
                if len(set(genotypes)) > 1:
                    vcf_out.write(line)


def main():

    # Args
    vcf_file = arg.vcf_infile

    if arg.remove_inv:
        remove_invariable(vcf_file)

main()
