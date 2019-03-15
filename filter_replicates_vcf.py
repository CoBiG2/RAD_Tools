#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 14 14:04:41 2019

@author: duartb
"""
from statistics import mode
import argparse

parser = argparse.ArgumentParser(prog='python3')

parser.add_argument("input_vcf", metavar="input.vcf", type=str,
                        help="VCF input file")
parser.add_argument("input_replicates", metavar="replicates.txt", type=str,
                        help="Input file with all the replicates' names")
parser.add_argument("output", metavar="output.vcf", type=str,
                        help="Output filtered VCF")
arguments = parser.parse_args()

def filter_replicate_vcf(vcf_input,replicates_input,vcf_output):
    
    vcf_file = open(vcf_input,"r")
    replicates_file = open(replicates_input,"r")
    output=open(vcf_output,"w")
    
    replicates = {}
    individuals=[]
    
    # cria um dicionário que associa o nome da amostra aos seus replicados
    for line in replicates_file:
        fields = line.split()
        individuals.append(fields[0])
        replicates[fields[0]] = fields[1:]
        
    replicate_number= len(fields)-1
    
    # creates header variable with the names of all replicates from the vcf    
    for line in vcf_file:
        if  line.startswith("#CHROM") == False:
            output.write(line)
        if line.startswith("#CHROM"):
            header = line.split()
            # stores the part of the header that doesn't contain the names of the replicates
            header_keep= header[0:9]
            break
    
    # creates a dictionary with the indexed position of all replicates in the header
    pos={}
    for rep_group in replicates.values():
        for rep in rep_group:
            pos[rep]= header.index(rep)
    
    # writes the new_header (joined replicates) to the output file
    new_header="\t".join(header_keep + individuals)
    output.write(new_header + "\n")
    
    genes=[]
    
    for line in vcf_file:
        
        write= True
        line = line.split()
        ind_dictio= {}    
        ind=0
        
        # iteares over the columns of a line in a vcf file
        for i in range(9,len(line)):
            genes.append(line[i].split(":")[0])
            # creates list of the genes in N columns, where N is the number of replicates used for each sample
            if len(genes) == replicate_number:
                # if an individual has contraditory information between its replicates the present snp is skipped
                if len(set(genes)) > 1 and "./." not in set(genes):
                    write= False
                    pass
                # generates dictiorary with genotypes of all individuals in present line
                else:
                    try:
                        ind_dictio[individuals[ind]]= mode(genes)
                    except:
                        ind_dictio[individuals[ind]]= "./."
                # resets list of replicate genes       
                genes=[]
                # keeps track of the index of the current individual
                ind += 1
        
        # line is writen only if none of the replicates had contraditory information
        if write == True:        
            new_line= line[0:9] + list(ind_dictio.values())
            line_write= "\t".join(new_line)
            output.write(line_write + "\n")
    
    output.close()
    
filter_replicate_vcf(arguments.input_vcf,arguments.input_replicates,arguments.output)