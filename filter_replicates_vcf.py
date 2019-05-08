#!/usr/bin/env python3

# Copyright 2019 Duarte Teomoteo Balata <duarte.balata@gmail.com>
# filter_replicates_vcf.py is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# filter_replicates_vcf is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with keep_central_snps. If not, see <http://www.gnu.org/licenses/>.

import argparse
from statistics import mode

parser = argparse.ArgumentParser(prog='python3')

parser.add_argument("input_vcf", metavar="input.vcf", type=str,
                        help="VCF input file")
parser.add_argument("input_replicates", metavar="replicates.txt", type=str,
                        help="Input file with all the replicates' names")
parser.add_argument("output", metavar="output.vcf", type=str,
                        help="Output filtered VCF")
parser.add_argument("-m","--missing_data_percentage", metavar="25",type=int, default=25,
                        help= "Maximum percentage of missing values allowed in a set of replicates (default=25%)")
parser.add_argument("-e","--error_percentage", metavar="10",type=int, default=10,
                        help= "Maximum percentage of individuals per SNP with erroneous replicates.")
arguments = parser.parse_args()

#missing_ratio= arguments.missing_data_percentage/100

def filter_replicate_vcf(vcf_input,replicates_input,vcf_output,missing_percentage,error_percentage):

    # keeps track of the total number of loci
    total=0
    # keeps track of the number of loci that are considered missing data
    missing=0
    # keeps track of the number of loci that are deleted
    deleted_loci=0
    # keeps track of the total number of loci
    loci=0
    #keeps track of ambiguities
    amb=0

    genes=[]

    # opens input vcf file
    vcf_file = open(vcf_input,"r")

    # open file with the replicate names association
    replicates_file = open(replicates_input,"r")

    #opens output file
    output=open(vcf_output,"w")

    replicates = {}
    individuals=[]

    # creates a dictionary that associates every sample name with the names of its replicates
    for line in replicates_file:
        fields = line.split()
        individuals.append(fields[0])
        replicates[fields[0]] = fields[1:]

    replicate_number= len(fields)-1

    # writes the unnaltered header lines of the vcf to output file
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

    # writes the new_header (per individual) to the output file
    new_header="\t".join(header_keep + individuals)
    output.write(new_header + "\n")

    # iterates over vcf file
    for line in vcf_file:

        #writing to output is enabled
        write= True
        error_count= 0
        line = line.split()
        ind_dictio= {}
        # keeps track of the index of the individual that the current gene group belongs to
        ind=0
        loci +=1

        # iteares over the columns of a line in a vcf file
        for i in range(9,len(line)):
            genes.append(line[i].split(":")[0])
            # creates list of the genes in N columns, where N is the number of replicates used for each sample
            if len(genes) == replicate_number:
                # if an individual has contraditory information between its replicates
                if (len(set(genes)) > 1 and "./." not in set(genes)) or len(set(genes)) > 2:
                    # writing is disabled for the current line
                    error_count += 1
                    if error_count > len(individuals) * error_percentage/ 100:
                        write = False
                else:
                    # if the percentage of missing data in the present group of replicates is
                    # below provided threshold, the allele information is added to the dictionary
                    if genes.count("./.") <= len(genes) * missing_percentage/ 100:
                        #while "./." in genes:
                            #genes.remove("./.")
                        try:
                            ind_dictio[individuals[ind]]= mode(genes)
                        except:
                            try:
                                genes.remove("./.")
                                ind_dictio[individuals[ind]]= mode(genes)
                            except:
                                amb += 1
                                ind_dictio[individuals[ind]]= "./."
                    # if the percentage of missing data in the present group of replicates is too
                    # high, the whole sample is considered missing data
                    else:
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

            # adds the number of samples in the present locus that were not considered missing data to the counter
            total += len(list(ind_dictio.values()))
            # adds the number of samples in the present locus that were considered missing data to the counter
            missing += list(ind_dictio.values()).count("./.")
        else:
            deleted_loci +=1

    output.close()

    print("\n"+ str(deleted_loci) + " loci out of " + str(loci) + " were removed due to contraditory information between replicates.")
    print(str(missing) + " Snp's out of " + str(total) + " were considered missing data." + "\n")
    print(str(amb) + " ambiguities found between replicates.")

filter_replicate_vcf(arguments.input_vcf,arguments.input_replicates,arguments.output, arguments.missing_data_percentage, arguments.error_percentage)
