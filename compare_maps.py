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

# This script will compare technical replicates on a vcf file and output a
# number of statistics and error rates. The vcf and pair mapping files must
# be provided

# compare_maps.py takes two SAM files produced by sequencing mapping tools,
# such as  bowtie, and compares the mapping location/depth (among other stats)
# between them. This will compare all non-header entries in the SAM files,
# so any previous filtering should be done with samtools.

# This was made with technical replicates in mind but can compare
# any pair of SAM files. Eventually multiple SAM files may be supported and
# compared through multiple regression analyses

import os
import argparse
from collections import defaultdict

parser = argparse.ArgumentParser(description="Comparison of RAD assemblies "
                                             "using technical replicates")

parser.add_argument("-in", dest="infiles", required=True, nargs="+",
                    help="Provide the SAM files separated by whitespace")
parser.add_argument("-mode", dest="mode", choices=["replicates", "assemblies"],
                    help="replicates: Provide a SAM file for each replicate;"
                          " assemblies: Provide a VCF file for each assembly")

arg = parser.parse_args()


def parse_vcf(vcf_file):
    """
    Parses a VCF file and retains the location of all SNPs in the genome as a list of tupples
    """

    vcf_fh = open(vcf_file)
    snp_locations = []

    for line in vcf_fh:
        #Skip header
        if line.startswith("#"):
            pass
        else:
            fields = line.strip().split()
            chrom = fields[0]
            pos = fields[1]
            snp_locations.append([chrom, pos])

    return snp_locations


def vcf_overlap(stor1, stor2):
    """
    Determines the SNPs and chromosome overlap between two vcf storage objects
    """

    # Get SNP position overlap and data set exclusives
    stor1_str = ["".join(x) for x in stor1]
    stor2_str = ["".join(x) for x in stor2]
    snp_overlap = len(set(stor1_str).intersection(set(stor2_str)))
    stor1_exclusive = len(set(stor1_str) - set(stor2_str))
    stor2_exclusive = len(set(stor2_str) - set(stor1_str))

    # Get chromosome overlap
    chrom1 = [x[0] for x in stor1]
    chrom2 = [x[0] for x in stor2]
    chrom_overlap = len(set(chrom1).intersection(set(chrom2)))

    log_handle = open("vcf_overlap.txt", "w")

    log_handle.write("SNP intersection stats:\nSNP overlap: {}\nExclusive to VCF1: {}\nExclusive to VCF2: "
                     "{}\n\nChromosome overlap: {}".format(snp_overlap, stor1_exclusive, stor2_exclusive, chrom_overlap))


def parse_sam(sam_file):
    """
    Parses a SAM file and retains each mapping location in the genome, along
    with the coverage.

    :param sam_file: string, path to SAM file

    :return map_locations: dictionary, with the genome contig/chromossome name
    plus position in the contig as value, and number of reads as value
    Example: map_locations = {"contig1_511": 14, "contig1_1902": 10}
    """

    sam_handle = open(sam_file)
    map_location = defaultdict(int)

    for line in sam_handle:

        # Skip headers
        if line.startswith("@"):
            continue

        else:
            fields = line.strip().split()
            genome_loc = fields[2]

            # Skip unmapped reads
            if genome_loc != "*":
                pos = fields[3]
                read = fields[0]

                map_location["{}_{}".format(genome_loc, pos)] += 1

    return map_location


def get_overlap(map_storage, output_file="map_overlap.csv"):
    """
    Takes a list of tuples with SAM file names as first element and their
    correspoding  dictionaries as second element. Returns basic statistics on
    the total  of mapping locations for each file, their overlap and the
    exclusive locations

    :param map_storage: list, SAM file name as 1st element; map_location
    dictionary as 2nd element
    """

    output_handle = open(output_file, "w")
    lstorage = []

    # Get individual stats
    output_handle.write("File; Mapping locations\n")
    for f, mstorage in map_storage:
        output_handle.write("{}; {}\n".format(f, len(mstorage)))

    # Get overlap and exclusive locations
    output_handle.write("\n")
    sam1 = set(map_storage[0][1])
    sam1_f = map_storage[0][0]
    sam2 = set(map_storage[1][1])
    sam2_f = map_storage[1][0]
    overlap = sam1.intersection(sam2)
    output_handle.write("Overlapping locations: {}({})\n".format(len(overlap),
                    float(len(overlap)) / float(min([len(sam1), len(sam2)]))))

    # Exclusive locations
    exclusive_sam1 = sam1.difference(sam2)
    exclusive_sam2 = sam2.difference(sam1)
    output_handle.write("Exclusive from {}: {}({})\n"
                        "Exclusive from {}: {}({})\n".format(
                    sam1_f, len(exclusive_sam1),
                    round(float(len(exclusive_sam1)) / float(len(sam1)), 2),
                    sam2_f, len(exclusive_sam2),
                    round(float(len(exclusive_sam2)) / float(len(sam2)), 2)))

    output_handle.close()


def main():
    # Arguments
    infiles = arg.infiles

    if arg.mode == "replicates":
        storage = []

        for f in infiles:
            map_loc = parse_sam(f)
            storage.append((os.path.basename(f), map_loc))

        get_overlap(storage)

    elif arg.mode == "assemblies":
        vcf1_storage = parse_vcf(infiles[0])
        vcf2_storage = parse_vcf(infiles[1])
        vcf_overlap(vcf1_storage, vcf2_storage)

main()
