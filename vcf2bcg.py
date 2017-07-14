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

# Conversion tool between the VCF format and bcg format

import argparse

parser = argparse.ArgumentParser(description="Comparison of RAD assemblies "
                                             "using technical replicates")

parser.add_argument("-vcf", dest="vcf_file", help="Provide vcf file",
                    required=True)
parser.add_argument("-p1", dest="parent_1", help="Provide file with taxa from"
                    " parent 1")
parser.add_argument("-p2", dest="parent_2", help="Provide file with taxa from"
                    " parent 2")
parser.add_argument("-hib", dest="hibrid", help="Provide file with taxa from"
                    " hybrids")

arg = parser.parse_args()

def parse_taxa_file(f):
    """
    Returns a list of taxa from a taxa file
    """

    fh = open(f)

    return [x.strip() for x in fh.readlines() if x.strip() != ""]


def convert_vcf(vcf_file, p1, p2, h):
    """
    Converts a vcf file into bcg format according to the taxa provided for
    parental populations (p1 and p2) and hybrid population
    """

    vcf_fh = open(vcf_file)

    parent1_fh = open("parent1.txt", "w")
    parent2_fh = open("parent2.txt", "w")
    hybrid_fh = open("hybrid.txt", "w")

    pos = {}

    c = 0

    for line in vcf_fh:
        if line.startswith("##"):
            pass
        elif line.startswith("#CHROM"):
            taxa_list = line.strip().split()

            # Check if all taxa in p1, p2 and h are in taxa_list
            for i in p1 + p2 + h:
                if i not in taxa_list:
                    print(taxa_list)
                    print("taxa {} not in vcf file".format(i))
                    raise SystemExit

        elif line.strip() != "":
            fields = line.strip().split()
            # Get genotype list for p1
            p1_geno = [fields[taxa_list.index(x)].split(":")[0] for x in p1]
            # Get genotype list for p2
            p2_geno = [fields[taxa_list.index(x)].split(":")[0] for x in p2]
            # Get genotype list for h
            h_geno = [fields[taxa_list.index(x)].split(":")[0] for x in h]

            # Write data for p1
            gen1_count = "".join(p1_geno).count("0")
            gen2_count = "".join(p1_geno).count("1")
            parent1_fh.write("locus_{}\n{} {}\n".format(c, gen1_count,
                                                        gen2_count))

            # Write data for p2
            gen1_count = "".join(p2_geno).count("0")
            gen2_count = "".join(p2_geno).count("1")
            parent2_fh.write("locus_{}\n{} {}\n".format(c, gen1_count,
                                                        gen2_count))

            # Write data for h
            hybrid_fh.write("locus_{}\npop_0\n".format(c))
            for i in h_geno:
                print(i)
                if [i.count("0"), i.count("1")] == [0, 0]:
                    gen = ["-9", "-9"]
                else:
                    gen = [i.count("0"), i.count("1")]
                hybrid_fh.write("{} {}\n".format(gen[0], gen[1]))

            c += 1

    # Close file handles
    parent1_fh.close()
    parent2_fh.close()
    hybrid_fh.close()


def main():

    # Args
    vcf_file = arg.vcf_file
    p1_file = arg.parent_1
    p2_file = arg.parent_2
    h_file = arg.hibrid

    p1 = parse_taxa_file(p1_file)
    p2 = parse_taxa_file(p2_file)
    h = parse_taxa_file(h_file)

    convert_vcf(vcf_file, p1, p2, h)

main()
