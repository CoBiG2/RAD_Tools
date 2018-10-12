#!/usr/bin/env python3

# Copyright 2018 Francisco Pina Martins <f.pinamartins@gmail.com>
# segregating_loci_finder.py is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# Loci_counter is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with Loci_counter. If not, see <http://www.gnu.org/licenses/>.

# This script will compare two groups of individuals and highlight any
# loci that segregate both groups

import re

def vcf_parser(vcf_filename, group_split_point):
    """
    Parses a vcf file and returns TODO
    """
    infile = open(vcf_filename, 'r')
    loci = {}
    for lines in infile:
        if lines.startswith("#"):  # Skip headers
            pass
        else:
            lines = lines.split()
            locus = lines[0]
            data = lines[8:]
            groups = []
            freqs = []
            groups.append(data[:group_split_point])
            groups.append(data[group_split_point:])
            abs_freqs = [re.match(".*:", x) for x in groups[0]]
