#!/usr/bin/python3
# Copyright 2015 Francisco Pina Martins <f.pinamartins@gmail.com>
# This file is part of filter_loci_by_vcf.
# filter_loci_by_vcf is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# filter_loci_by_vcf is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with filter_loci_by_vcf.  If not, see <http://www.gnu.org/licenses/>.

# Discards everything in a .loci file that is not present in a .vcf file.

import Loci_filter_from_vcf


def main(loci_filename, vcf_filename):
    """
    Handle everything.
    """
    loci_list = Loci_filter_from_vcf.VCF_parser(vcf_filename)
    loci_parser(loci_filename, loci_list)


def loci_parser(loci_filename, loci_list):
    """
    Parses the .loci file and prints the lines that belong to a loci from the
    vcf file.
    """

    infile = open(loci_filename, 'r')

    if 1 in loci_list:
        interested = True
    else:
        interested = False

    for lines in infile:
        if interested == True:
            print(lines, end="")
        if lines.startswith("//"):
            loci_num = int(lines.split("|")[1])
            if (loci_num + 1) in loci_list:
                interested = True
            else:
                interested = False



if __name__ == "__main__":
    # Usage python3 filter_loci_by_vcf.py file.loci file.vcf
    from sys import argv
    main(argv[1], argv[2])
