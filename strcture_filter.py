#!/usr/bin/python3

# Copyright 2015 Francisco Pina Martins <f.pinamartins@gmail.com>
# This file is part of structure_filter.
# structure_filter is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# structure_filter is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with structure_filter. If not, see <http://www.gnu.org/licenses/>.

#Usage:python3 structure_filter.py subset_file.txt file.structure > file.filtered.structure

def get_list_of_subset(subset_filename):
    """Reads the subset loci from a file and returns a list with it."""
    subset_file = open(subset_filename, 'r')
    subset = []
    for lines in subset_file:
        subset.append(lines.strip())
    
    subset_file.close()
    return tuple(subset)


def structure_breaker(structure_filename, subset):
    """Filters the columns from a structure file based on the subset tuple."""
    structure_file = open(structure_filename, 'r')
    indeces = []
    header = [" ", " "] + structure_file.readline().strip().split()
    for snps in subset:
        indeces.append(header.index(snps))
    indeces.sort()
    print(" ".join(header[0:2] + [header[i] for i in indeces]))
    for lines in structure_file:
        lines = lines.strip().split()
        print(" ".join(lines[0:2] + ["\t"] + [lines[i] for i in indeces]))
    
    structure_file.close()


if __name__ == "__main__":
    from sys import argv
    subset = get_list_of_subset(argv[1])
    structure_breaker(argv[2], subset)
    
