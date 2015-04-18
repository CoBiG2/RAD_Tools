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
    """Parses a VCF file and returns a sorted list with loci names and a list
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
    """Gets a loci list, and a loci file and sequence names and filters the loci
    file according to the loci list. Returns a dict {seqname: sequence}"""
    loci_file = open(loci_filename, 'r')

    if loci[0] == "0":
        estouinteressado = 1  
    else:
        estouinteressado = 0  
    
    c = 0
    seqlen = 0
    totlen = 0
    vcfseqs = set(seqnames.keys())
    #taxaset = set(vcfseqs)
    taxaset = set()
    seqlines= "  "
    for lines in loci_file:
        if estouinteressado == 1 and lines.startswith("//") == False:

            seqlines = lines.strip(">\n").split()

            seqnames[seqlines[0]] += seqlines[1]
            taxaset.add(seqlines[0])

        elif lines.startswith("//") and estouinteressado == 1:
            seqlen = len(seqlines[1])
            totlen += seqlen

            difset = vcfseqs.difference(taxaset)

            for t in difset:
                seqnames[t] += "N" * seqlen

            taxaset = set()
            c += 1
            estouinteressado = 0

            if str(c) in loci:
                estouinteressado = 1
        
        elif estouinteressado == 0 and lines.startswith("//"):
            c += 1
            if str(c) in loci:
                estouinteressado = 1
                
    loci_file.close()

    return seqnames


def phy_writer(phy_filename, seqnames):
    """Writes the output ready to submit to RAxML or other phylogeny
    program. Based on seqnames dict {seqname: sequence}"""
    phy = open(phy_filename, 'w')
    seqnum = len(seqnames)
    bpnum = len(list(seqnames.values())[0])
    phy.write(str(seqnum) + " " + str(bpnum) + "\n")
    for k, v in seqnames.items():
        phy.write(k + "\t" + v + "\n")
        
    phy.close()

if __name__ == "__main__":
    from sys import argv
    loci, seqnames = vcf_parser(argv[1])
    seqnames = loci_parser(argv[2], loci, seqnames)
    phy_writer(argv[3], seqnames)
