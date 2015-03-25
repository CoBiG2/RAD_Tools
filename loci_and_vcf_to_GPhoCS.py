#!/usr/bin/python3
# Copyright 2015 Francisco Pina Martins <f.pinamartins@gmail.com>
# This file is part of loci_and_vcf_toGPhoCS.
# loci_and_vcf_toGPhoCS is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# loci_and_vcf_toGPhoCS is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with loci_and_vcf_toGPhoCS.  If not, see <http://www.gnu.org/licenses/>.

# Usage: python3 loci_and_vcf_toGPhoCS.py file.vcf file.loci file.GPhoCS

def vcf_parser(vcf_filename):
    """Parses a VCF file and returns a sorted list with loci names"""
    vcf = open(vcf_filename, 'r')
    loci = []

    for line in vcf:
        if line.startswith("#") == False:
            loci.append(line.split()[0])
    
    vcf.close()
    loci = sorted(list(set(loci))) 
    
    return loci


def GPhoCS_writer(loci_filename, GPhoCS_filename, loci):
    """Gets a loci list, and a loci file and filters it. It then saves the data
    in GPhoCS format."""
    loci_file = open(loci_filename, 'r')
    GPhoCS = open(GPhoCS_filename, 'w')
    
    if loci[0] == "0": 
        estouinteressado=1  
    else:
        estouinteressado=0  
    
    c=0
    sequences=[]
    GPhoCS.write(str(len(loci)) + "\n\n")
    for lines in loci_filename:
        if estouinteressado==1 and lines.startswith("//")==False:
            lines=lines.replace("-","N")
            sequences.append(lines)
        elif lines.startswith("//"):
            if sequences!=[]:
                name=str(c)
                name=str(c)
                numbseq=str(len(sequences))
                seqlen=str(len(sequences[0].split()[1]))
                GPhoCS.write(name + " " + numbseq + " " + seqlen + "\n")
                for seqs in sequences:
                    GPhoCS.write(seqs)
                GPhoCS.write("\n")
            c+=1
            estouinteressado=0
            sequences=[]
            if str(c) in loci:
                estouinteressado=1
                
    loci_file.close()
    GPhoCS.close()


if __name__ == "__main__":
    from sys import argv
    loci = vcf_parser(argv[1])
    GPhoCS_writer(argv[2], argv[3], loci)
