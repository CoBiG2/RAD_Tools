#!/usr/bin/python3

# Copyright 2015 Francisco Pina Martins <f.pinamartins@gmail.com>
# This file is part of Loci_filter_from_vcf.
# Loci_filter_from_vcf is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# Loci_filter_from_vcf is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with Loci_filter_from_vcf. If not, see <http://www.gnu.org/licenses/>.

#Usage:python3 Loci_filter_from_vcf.py file.vcf file.loci filtered_file.loci

def VCF_parser(vcf_filename):
    """Parse the VCF and return a list with the loci that contain SNPs."""
    infile = open(vcf_filename, 'r')
    loci = []
    for lines in infile:
        if lines.startswith("#") == False:
            loci.append(int(lines.split()[0]))
    loci = sorted(list(set(loci)))
    
    return loci
    
if __name__ == "__main__":
    from sys import argv
    print(VCF_parser(argv[1]))
