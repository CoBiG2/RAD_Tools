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

# Setting plot style
plt.style.use("ggplot")

parser = argparse.ArgumentParser(description="Parsing statistical output of"
                                 " VCFtools")

parser.add_argument("-in", dest="infile", help="The vcftools output file",
                    required=True)
parser.add_argument("--weir-fst", dest="weir_fst", help="Parse a Weir FST"
                    " output file", const=True, action="store_const")
parser.add_argument("--singleton", dest="singleton", const=True,
                    action="store_const", help="Parse a singleton/doubleton"
                    " output file")

arg = parser.parse_args()


def weir_fst(infile):
    """
    Parses the Weir and Cockerham FST output file
    """

    fst_vals = []

    with open(infile) as fh:

        # Skip header
        next(fh)

        for line in fh:
            if line.strip().split()[2] != "-nan":
                fst_vals.append(float(line.strip().split()[2]))

    # Creating plots
    # Histogram
    f, ax = plt.subplots()

    plt.hist(fst_vals)

    f.tight_layout()
    plt.savefig("fst_distribution.png")
    plt.close()

    # Plot
    f, ax = plt.subplots()

    plt.plot(fst_vals, "bo")

    f.tight_layout()
    plt.savefig("fst_vals.png")


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


def main():

    # Arguments
    infile = arg.infile

    if arg.weir_fst:
        weir_fst(infile)

    if arg.singleton:
        singletons(infile)

main()
