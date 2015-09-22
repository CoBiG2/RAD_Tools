#!/usr/bin/python3
# Copyright 2015 Francisco Pina Martins <f.pinamartins@gmail.com>
# This file is part of bayenv2_results_miner.
# bayenv2_results_miner is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# bayenv2_results_miner is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with bayenv2_results_miner.  If not, see <http://www.gnu.org/licenses/>.

# Usage: bayenv2_results_miner.py file.bf env_vars.txt [BF_value] [p-value]


import re


def main(bf_filename, envfile_name, bayesfactor=10, p_value=0.05):
    """
    Parses the ".bf" file outputted by Bayenv2 and returns the SNPs with
    putative associations.
    """

    env_vars = parse_env_vars(envfile_name)

    infile = open(bf_filename, 'r')

    for lines in infile:
        lines = lines.strip().split()
        snp_name = name_sanitizer(lines[0])
        bfs = [float(x) for x in lines[1::3]]
        spearman = [float(x) for x in lines[2::3]]
        # pearson = [float(x) for x in lines[3::3]] # Unused so far

        for i in range(len(bfs)):
            if bfs[i] >= bayesfactor and spearman[i] < p_value:
                print(snp_name + ": " + env_vars[i])

    infile.close()

def name_sanitizer(snp_name):
    """
    Replaces the full text on the SNP name for something sensible & returns it.
    """
    replacement = re.search("/\\d+.", snp_name).group(0)[1:-1]
    replacement = "SNP_" + replacement

    return replacement

def parse_env_vars(envfile_name):
    """
    Pareses a txt file with one environment variable name per line and returns
    a list with these values.
    """
    env_vars = []
    infile = open(envfile_name, 'r')
    for lines in infile:
        env_vars.append(lines.strip())

    return env_vars

if __name__ == "__main__":
    from sys import argv
    ARGUMENTS = argv[1:]
    main(*ARGUMENTS)
