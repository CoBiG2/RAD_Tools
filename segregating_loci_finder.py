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
from collections import Counter, OrderedDict


def vcf_parser(vcf_filename, group_split_point):
    """
    Parses a vcf file and returns TODO
    """
    infile = open(vcf_filename, 'r')
    loci = {}
    group_split_point = int(group_split_point)
    for lines in infile:
        if lines.startswith("#"):  # Skip headers
            if lines.startswith("#CHROM"):  # Group checker lines
                lines = lines.split()
                data = lines[9:]
                groups = [data[:group_split_point],
                          data[group_split_point:]]
                print(groups)
        else:
            lines = lines.split()
            locus = lines[0]
            data = lines[9:]
            groups = [data[:group_split_point], data[group_split_point:]]
            gr_freqs = [get_freqs(x) for x in groups]
            loci[locus] = gr_freqs

    return loci


def get_freqs(vcf_data):
    """
    Gets relative frequencies from VCF data
    """
    abs_freqs = [re.match(".*?:", x).group(0)[:-1] for x in vcf_data]
    dummy_freqs = {"0/0": 0, "0/1": 0, "1/0": 0, "1/1": 0, "./.": 0}
    rel_freqs = Counter(abs_freqs)
    try:
        mvs = rel_freqs.pop("./.")
    except KeyError:
        mvs = 0
    dummy_freqs.update(Counter(abs_freqs))
    rel_freqs = dummy_freqs
    rel_freqs["0/1"] += rel_freqs.pop("1/0")
    try:
        non_missing = len(abs_freqs) - mvs
        rel_freqs = {k: v/non_missing for k, v in rel_freqs.items()}
    except ZeroDivisionError:
        rel_freqs = None
    # print(rel_freqs)

    return rel_freqs


def segregating_freqs(loci):
    """
    Defines wether a locus segregates the two groups
    For now only works with full segregation
    """
    segregators = []
    for locus, data in loci.items():
        try:
            segregators += [locus for k, v in data[0].items()
                            if (data[1][k] == 0 and v == 1)
                            or (data[1][k] == 1 and v == 0)]
        except AttributeError:
            pass

    return segregators


if __name__ == "__main__":
    from sys import argv
    seg_loci = vcf_parser(argv[1], argv[2])
    for i in segregating_freqs(seg_loci):
        print(i)
