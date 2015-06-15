#!/usr/bin/python3

#Usage:python3 VCF_corresponder.py infile.vcf > outfile.modified_vcf

def VCF_modifier(vcf_filename):
    """Grabs a VCF file and adds "SNP_##" to the end of the line, to match
    those outputed by PDGspider when converting to to other formats."""
    vcf = open(vcf_filename)
    counter = 1
    for lines in vcf:
        if lines.startswith("#") == False:
            lines = lines.strip() + "\tSNP_" + str(counter) + "\n"
            counter += 1
        print(lines, end="")


if __name__ == "__main__":
    from sys import argv
    VCF_modifier(argv[1])

#What is VCF_corresponder.py good for?
#To produce subsets of loci under positive or balancing selection, or evolving under neutrality, downstream of LOSITAN (Beaumont and Nichol, 1996; Antao et al., 2008) utilization you use the following line of Shell code: 
#1. Edit your loci_list file adding "Positive" "Neutral" and "Balancing" terms to the corresponding locus
#(You can use our script FDR.R to correct your loci_list from LOSITAN) 
#2. grep "Balancing" loci_list.csv | sed 's/,.*//g' > Balancing.txt

#3.for i in $(cat Balancing.txt)
#do
#grep -e "${i}$"  outfile.modified.vcf |sed 's/\tSNP.*//g' >> Balancing.vcf
#done

#4.If you are gonna use these on STRUCTURE or other software, you should do:
#sort -n -k 1,2 Balancing.vcf > Balacing_sorted.vcf

#5. add the typical .vcf header to your files with: 
#cat header.vcf Balancing_sorted.vcf >> Balancing_ready.vcf


