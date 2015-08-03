#!/usr/bin/python3
# Copyright 2015 Francisco Pina Martins <f.pinamartins@gmail.com>
# This file is part of loci2phy.
# loci2phy is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# loci2phy is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with loci2phy.  If not, see <http://www.gnu.org/licenses/>.

# Usage: python3 loci2phy.py file.vcf file.loci file.phy

def vcf_parser(vcf_filename):
    """Parses a VCF file and returns a dict with loci names and a sortd list
    with taxa names."""
    vcf = open(vcf_filename, 'r')
    loci = []
    seqnames = {}

    for line in vcf:
        if line.startswith("##"):
            pass
        elif line.startswith("#CHROM"):
            for names in line.split()[9:]:
                seqnames[names] = ""
        else:
            loci.append(line.split()[0])

    vcf.close()
    loci = sorted(list(set(loci)))

    return loci, seqnames


def loci_parser(loci_filename, loci, seqnames):
    """Gets a loci list, a .loci file and sequence names and filters the .loci
    file according to the loci list. Returns a dict {seqname: sequence}"""
    loci_file = open(loci_filename, 'r')

    if loci[0] == "1":
        gather_stuff = 1
    else:
        gather_stuff = 0

    seqlen = 0
    totlen = 0
    locus_number = 1
    vcfseqs = set(seqnames.keys())
    taxaset = set()
    seqlines= "  "
    seqlens = [0]
    for lines in loci_file:
        if gather_stuff == 1 and lines.startswith("//") == False:

            seqlines = lines.strip(">\n").split()

            if seqlines[0] in seqnames:
                seqnames[seqlines[0]] += seqlines[1]
                taxaset.add(seqlines[0])

        elif lines.startswith("//") and gather_stuff == 1:
            seqlen = len(seqlines[1])
            totlen += seqlen
            seqlens.append(totlen)
            
            difset = vcfseqs.difference(taxaset)

            for t in difset:
                seqnames[t] += "N" * seqlen

            print(locus_number)

            taxaset = set()
            gather_stuff = 0

            try:
                 if str(int(lines[lines.find("|"):]) + 1) in loci:
                     gather_stuff = 1
            except:
                 locus_number += 1
                 if str(locus_number) in loci:
                    gather_stuff = 1

        elif gather_stuff == 0 and lines.startswith("//"):

            try:
                if str(int(lines[lines.find("|"):]) + 1) in loci:
                    gather_stuff = 1
            except:
                locus_number += 1
                if str(locus_number) in loci:
                    gather_stuff = 1

    loci_file.close()

    return seqnames, seqlens


def phy_writer(phy_filename, seqnames, seqlens):
    """Writes the output ready to submit to RAxML or other phylogeny
    programs, along with a partitions file. Based on seqnames dict 
    {seqname: sequence} and on [seqlens]"""
    phy = open(phy_filename, 'w')

    seqnum = len(seqnames)
    bpnum = len(list(seqnames.values())[0])
    phy.write(str(seqnum) + " " + str(bpnum) + "\n")
    for k, v in seqnames.items():
        phy.write(k + "\t" + v + "\n")

    phy.close()
    
    part = open(phy_filename[:-4] + ".part", 'w')
    for counts in list(range(len(seqlens)))[:-1]:
        part.write("DNA, p" + str(counts + 1) + "=" + str(seqlens[counts] + 1) + "-" + str(seqlens[counts + 1]) + "\n")
    
    part.close()

if __name__ == "__main__":
    from sys import argv
    loci, seqnames = vcf_parser(argv[1])
    seqnames, seqlens = loci_parser(argv[2], loci, seqnames)
    phy_writer(argv[3], seqnames, seqlens)
