#!/usr/bin/python3

# Copyright 2015 Diogo N. Silva <o.diogosilva@gmail.com>
# Copyright 2017-2018 Francisco Pina-Martins <f.pinamartins@gmail.com>
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
import re

PARSER = argparse.ArgumentParser(description="Creates a consensus sequence for"
                                             " each locus in a .loci file")

PARSER.add_argument("-in", dest="loci_file", help="Provide the path to .loci "
                    "file", required=True)
PARSER.add_argument("-o", dest="output_file", help="Name of the output"
                    " file containing the consensus sequences")
PARSER.add_argument("-l", dest="loci_list_file",
                    help="Provide the path to file with loci numbers.",
                    default=None)

ARG = PARSER.parse_args()

IUPAC = {"AG": "R", "CT": "Y", "CG": "S", "AT": "W", "GT": "K", "AC": "M",
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
        column = {x for x in column if x not in ("-", "N")}

        if len(column) == 1:
            con += "".join(column)

        elif len(column) > 1:
            try:
                amb = "".join(sorted(column))
                con += IUPAC[amb]
            except KeyError:
                con += "N"

    return con


def list_parser(list_filename):
    """
    Parses a file with a loci list (one entry per line).
    Returns a set with the loci numbers on the list.
    """
    loci_file = open(list_filename, "r")
    loci_set = {}
    for lines in loci_file:
        lines = lines.strip().split()
        locus_num_id = int(re.search(r"\d+$", lines[0]).group()) - 1
        loci_set[str(locus_num_id)] = lines[-1]

    loci_file.close()

    return loci_set


def standard_write(output_handle, locus_number, con_seq, loci_set):
    """
    Standard fasta writter function. To be used when no loci list is provided.
    """
    locus_name = str(int(locus_number) + 1)
    output_handle.write(">vcf_locus{}\n{}\n".format(locus_name, con_seq))


def conditional_write(output_handle, locus_number, con_seq, loci_set):
    """
    Writter function to be used when a loci list is provided.
    """
    if locus_number in loci_set:
        locus_name = str(int(locus_number) + 1)
        output_handle.write(">vcf_locus{0}#{1}\n{2}"
                            "\n".format(locus_name,
                                        loci_set[locus_number],
                                        con_seq))


def create_consensus(loci_file, output_file, loci_list_file):
    """
    The parsing of the .loci file and creation and writting of consensus
    sequences to output_file are performed together to improve performance.

    :param loci_file: string, path to .loci file
    :param output_file: string, path to output file
    """

    loci_handle = open(loci_file)
    output_handle = open(output_file, "w")

    if loci_list_file is not None:
        loci_set = list_parser(loci_list_file)
        writer_func = conditional_write
    else:
        loci_set = None
        writer_func = standard_write

    locus_storage = []

    for line in loci_handle:

        if line.strip().startswith("//") and locus_storage:
            # Get locus number
            locus_number = line.split("|")[-2].strip()
            # Get sequence consensus
            con_seq = consensus(locus_storage)
            # Write to output
            writer_func(output_handle, locus_number, con_seq, loci_set)

            # Reset locus storage. Ready for next locus
            locus_storage = []

        else:
            locus_storage.append(line.split()[-1].strip())


def main():
    # Args
    loci_file = ARG.loci_file
    output_file = ARG.output_file
    loci_list_file = ARG.loci_list_file

    create_consensus(loci_file, output_file, loci_list_file)


main()
