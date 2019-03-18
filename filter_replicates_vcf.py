#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 14 14:04:41 2019

@author: duartb
"""
import argparse

parser = argparse.ArgumentParser(prog='python3')

parser.add_argument("input_vcf", metavar="input.vcf", type=str,
                        help="VCF input file")
parser.add_argument("input_replicates", metavar="replicates.txt", type=str,
                        help="Input file with all the replicates' names")
parser.add_argument("output", metavar="output.vcf", type=str,
                        help="Output filtered VCF")
parser.add_argument("-m","--missing_data_percentage", metavar="25",type=int, default=25,
                        help= "Maximum percentage of missing values allowed in a set of replicates (default=25%")
arguments = parser.parse_args()

missing_ratio= arguments.missing_data_percentage/100

def filter_replicate_vcf(vcf_input,replicates_input,vcf_output):
    
    non_missing=0
    missing=0 
    copies=0
    
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
    deleted_loci=0
    loci=0
    
    for line in vcf_file:
         
        write= True
        line = line.split()
        ind_dictio= {}    
        ind=0
        loci +=1
        
        # iteares over the columns of a line in a vcf file
        for i in range(9,len(line)):
            genes.append(line[i].split(":")[0])
            # creates list of the genes in N columns, where N is the number of replicates used for each sample
            if len(genes) == replicate_number:
                # if an individual has contraditory information between its replicates the present snp is skipped
                if (len(set(genes)) > 1 and "./." not in set(genes)) or len(set(genes)) > 2:
                    write = False
                # generates dictiorary with genotypes of all individuals in present line
                else:
                    if genes.count("./.") <= len(genes) * missing_ratio:
                        while "./." in genes: genes.remove("./.")
                        ind_dictio[individuals[ind]]= genes[0]
                    else:
                        ind_dictio[individuals[ind]]= "./."
                # resets list of replicate genes       
                genes=[]
                # keeps track of the index of the current individual
                ind += 1
        
        # line is writen only if none of the replicates had contraditory information
        if write == False:
            deleted_loci +=1
        
        if write == True:
            new_line= line[0:9] + list(ind_dictio.values())
            line_write= "\t".join(new_line)
            output.write(line_write + "\n")
            
            non_missing += (list(ind_dictio.values()).count("0/1") + list(ind_dictio.values()).count("1/1") + list(ind_dictio.values()).count("0/0"))
            missing += list(ind_dictio.values()).count("./.")
            copies += len(ind_dictio.values())
    
    output.close()
    
    print("\n"+ str(deleted_loci) + " loci out of " + str(loci) + " were removed due to contraditory information between replicates.")
    print(str(missing) + " Snp's out of " + str(non_missing + missing) + " were considered missing data." + "\n")
    
filter_replicate_vcf(arguments.input_vcf,arguments.input_replicates,arguments.output)