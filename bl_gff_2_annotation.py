#!/usr/bin/python

# Copyright 2017 Francisco Pina-Martins <f.pinamartins@gmail.com>
# bl_gff_2_annotation.py is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# Loci_counter is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with Loci_counter. If not, see <http://www.gnu.org/licenses/>.

from collections import OrderedDict


def blast_tab_parser(blast_tab_filename):
    """
    Parses a BLAST tabular format with the following fields:
    "qseqid sseqid qlen evalue nident length sstart send"
    Returns a dict identified by "qseqid" whose value is a list with the other
    fields.
    """
    infile = open(blast_tab_filename, "r")
    blast_dict = {}
    for lines in infile:
        lines = lines.strip().split()
        blast_dict[lines[1]] = [lines[0]] + lines[2:]

    infile.close()
    return blast_dict


def gff_parser(gff_filename, blast_data):
    """
    Parses a GFF file to look for the annotation corresponding to the zone
    provided in the BLAST data dictionary.
    Prints the output in the process.
    """
    # Print the header
    print("SNP_name\tScaffold name\tSequence length\tevalue\tidentity\tmatch "
          "length\tscaffold start\tscaffold end\tdatabase reference\tNote\t"
          "GO term")
    infile = open(gff_filename, "r")
    for lines in infile:
        lines = lines.strip().split("\t")
        if lines[0] in blast_data:
            annot_range = [int(lines[3]), int(lines[4])]
            seq_range = [int(blast_data[lines[0]][-2]),
                         int(blast_data[lines[0]][-1])]
            if min(seq_range) > min(annot_range) and \
                    max(seq_range) < max(annot_range):
                annot = lines[8].split(";")
                tags = OrderedDict([("dbxref", ""),
                                    ("note", ""),
                                    ("ontology_term", "")])
                for param in annot:
                    if param.lower().startswith(tuple(tags.keys())):
                        tags[param[:param.index("=")].lower()] = param
                annot = [x for x in tags.values()]
                if tags["note"] != "":
                    print("{0}\t{1}\t{2}\t{3}"
                          "".format(blast_data[lines[0]][0],
                                    lines[0],
                                    "\t".join(blast_data[lines[0]][1:]),
                                    "\t".join(annot)))


if __name__ == "__main__":
    from sys import argv
    BLAST_DATA = blast_tab_parser(argv[1])
    gff_parser(argv[2], BLAST_DATA)
