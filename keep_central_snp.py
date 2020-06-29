#!/usr/bin/env python3

# Copyright 2019 Duarte Teomoteo Balata <duarte.balata@gmail.com>
# keep_central_snps.py is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# keep_central_snps is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with keep_central_snps. If not, see <http://www.gnu.org/licenses/>.

import argparse
import sys


def get_args(args):
    parser = argparse.ArgumentParser(prog='python3')

    parser.add_argument("input", metavar="input.vfc", type=str,
                        help="VCF file with all the SNPs.")

    parser.add_argument("-o", "--output", metavar="output.vcf", dest="output", type=str,
                        help="output file with the SNPs closest to the centre of each locus (default = ./center_snps_filtered.vcf)")

    parser.add_argument("-l", "--len", metavar="int", dest="length", type=str,
                        help="length of each locus (required)", required=True)

    arguments = parser.parse_args(args)
    return arguments

# writes vcf header lines to output


def write_vcf_headers(vcf_path, output_path):
    vcf_file = open(vcf_path, "r")
    out_file = open(output_path, "w")

    for line in vcf_file:
        if line.startswith("#"):
            out_file.write(line)
        else:
            break

    vcf_file.close()
    out_file.close()

# Chooses the most central SNP based on loci ID and SNP read position


def write_vcf_body(vcf_path, output_path, locus_length):

    vcf_file = open(vcf_path, "r")
    out_file = open(output_path, "a")

    buff = []
    dist = 0
    last_id = ""
    last_name = ""
    loci_id = ""

    for line in vcf_file:
        try:
            name = line.split("\t")[0]
            dist = int(line.split("\t")[1])
            loci_id = line.split("\t")[2].split(":")[0]

            # if locus/ chromossome is different from last line, writes the
            # previous best snps and creates a new list
            if last_id != loci_id or last_name != name:

                if len(buff) > 0:
                    prev_best = int(locus_length)/2
                    for loci in buff:
                        loci_position = int(loci.split("\t")[2].split(":")[1])
                        # the snp ocurring closer to the read center is considered optimal.
                        if abs(int(locus_length)/2-loci_position) <= prev_best:
                            center = loci
                            prev_best = abs(int(locus_length)/2-loci_position)

                    out_file.write(center)

                buff = []
                buff.append(line)
            else:
                buff.append(line)

            last_id = loci_id
            last_name = name

        except (ValueError, IndexError) as error:
            pass

    # writes last iteration
    prev_best = int(locus_length)/2
    for loci in buff:
        loci_position = int(loci.split("\t")[2].split(":")[1])
        if abs(int(locus_length)/2-loci_position) < prev_best:
            center = loci
            prev_best = abs(int(locus_length)/2-loci_position)

    out_file.write(center)

    vcf_file.close()
    out_file.close()


if __name__ == "__main__":
    args = get_args(sys.argv[1:])
    try:
        write_vcf_headers(args.input, args.output)
        write_vcf_body(args.input, args.output, args.length)
    except IndexError:
        write_vcf_headers(args.input, "central_snps_filtered.vcf")
        write_vcf_body(args.input, "central_snps_filtered.vcf", args.length)
