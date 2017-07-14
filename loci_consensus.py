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

# loci_consensus will collapse sequences from each locus in a PyRAD .loci file
# into a consensus, witch SNPs being replaced by ambiguity codes


import argparse

parser = argparse.ArgumentParser(description="Creates a consensus sequence for"
                                             " each locus in a .loci file")

parser.add_argument("-in", dest="loci_file", help="Provide the path to .loci "
                    "file", required=True)
parser.add_argument("-o", dest="output_file", help="Name of the output"
                    " file containing the consensus sequences")

arg = parser.parse_args()

iupac = {"AG": "R", "CT": "Y", "CG": "S", "AT": "W", "GT": "K", "AC": "M",
         "CGT": "B", "AGT": "D", "ACT": "H", "ACG": "V", "ACGT": "N"}


def consensus(sequence_list):
    """
    Takes a list of sequences and returns a single consensus sequence. There is
    an assumption that sequences have the same sequence lenght.

    :param sequence_list: list, with sequence strings as elements

    :returns con: string, the consensus sequence of sequence_list
    """

    con = ""

    for column in zip(*sequence_list):

        # Remove gaps, missing characters and duplicate entries in column
        column = set([x for x in column if x != "-" and x != "N"])

        if len(column) == 1:
            con += "".join(column)

        elif len(column) > 1:
            try:
                amb = "".join(sorted(column))
                con += iupac[amb]
            except KeyError:
                con += "N"

    return con


def create_consensus(loci_file, output_file):
    """
    The parsing of the .loci file and creation and writting of consensus
    sequences to output_file are performed together to improve performance.

    :param loci_file: string, path to .loci file
    :param output_file: string, path to output file
    """

    loci_handle = open(loci_file)
    output_handle = open(output_file, "w")

    locus_storage = []

    for line in loci_handle:

        if line.strip().startswith("//") and locus_storage:
            # Get locus number
            locus_number = line.split("|")[-1].strip()
            # Get sequence consensus
            con_seq = consensus(locus_storage)
            # Write to output
            output_handle.write(">locus{}\n{}\n".format(locus_number, con_seq))

            # Reset locus storage. Ready for next locus
            locus_storage = []

        else:
            locus_storage.append(line.split()[-1].strip())

def main():
    # Args
    loci_file = arg.loci_file
    output_file = arg.output_file

    create_consensus(loci_file, output_file)


main()
