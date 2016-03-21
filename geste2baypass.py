#!/usr/bin/python3

# Copyright 2016 Francisco Pina Martins <f.pinamartins@gmail.com>
# This file is part of geste2baypass.
# geste2baypass is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# geste2baypass is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with geste2baypass. If not, see <http://www.gnu.org/licenses/>.

#Usage: python3 geste2baypass.py

# Shamelessly stole code from  https://github.com/Telpidus/omics_tools/blob/master/GESTE2BayEnv.py

from collections import OrderedDict

def parse_geste(infile_name): # Parses the  file (Bayscan input)
    """
    Parses a GESTE file and retuns a dict with:
    {"SNP_num":["FreqAp1\tFreqBp1","FreqAp2\tFreqBp2",...]}.
    """
    infile = open(infile_name, "r")
    snp = OrderedDict()
    for line in infile:
        line = line.split()
        # Neat trick to ignore data that is not SNP info
        try:
            int(line[0])
        except (ValueError, IndexError):
            continue
        try:
            snp[line[0]].append("\t".join(line[3:]))
        except KeyError:
            snp[line[0]] = ["\t".join(line[3:])]

    infile.close()

    return snp


def write_baypass(snp, baypass_filename):
    """
    Write a Baypass inpt file based on the GESTE dict.
    """
    outfile = open(baypass_filename, "w")
    for snps in snp.values():
        line = "\t".join(snps) + "\n"
        outfile.write(line)

    outfile.close()



if __name__ == "__main__":
    from sys import argv
    snp_dict = parse_geste(argv[1])
    write_baypass(snp_dict, argv[2])
