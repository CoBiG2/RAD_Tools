#!/usr/bin/python3

# Usage:python3 VCF_corresponder.py loci_list.txt infile.vcf > loci_list_matching_vcf.txt


def vcf_parser(vcf_filename, loci_set):
    """
    Grabs a VCF file and a loci list and matches the SNP order number to that
    of the VCF.
    If they match, the "real" SNP name is retrieved.
    Prints all the loci "true" names.
    """
    vcf = open(vcf_filename)
    counter = 1
    for lines in vcf:
        if lines.startswith("#") is False:
            name = lines.split()[0]
            if str(counter) in loci_set:
                print(name)
            counter += 1

    vcf.close()


def list_parser(loci_list_filename):
    """
    Parses a list with loci order numbers and returns a set with those values.
    """
    loci_file = open(loci_list_filename, "r")
    loci_set = set()
    for line in loci_file:
        loci_set.add(line.strip())

    loci_file.close()

    return loci_set

if __name__ == "__main__":
    from sys import argv
    LOCI_SET = list_parser(argv[1])
    vcf_parser(argv[2], LOCI_SET)
