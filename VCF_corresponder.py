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
