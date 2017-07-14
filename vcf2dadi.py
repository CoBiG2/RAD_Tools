#!/usr/bin/python

# Simple conversion tool from VCF to dadi's SNP input format

import vcf
import argparse
from collections import OrderedDict

parser = argparse.ArgumentParser(description="Conversion tool from VCF to SNP"
					      " input format from dadi")

parser.add_argument("-vcf", dest="vcf_file", help="Path to vcf file")
parser.add_argument("-p", dest="population_file", help="Populations file. There"
		    " should be one line per population, with the name of the"
		    " population separated from the samples by a colon (':') "
		    " and each sample separated by whitespace")
parser.add_argument("-o", dest="output_file", help="Name of the output file")

arg = parser.parse_args()


def parse_populations(pop_file):
    """
    Parses a population file and returns an orderedDict object with the
    population name as keys and the corresponding list of samples as values
    """

    pops = OrderedDict()

    fh = open(pop_file)

    for line in fh:
        fields = line.split(":")
        pops[fields[0]] = [x for x in fields[1].split() if x]

    return pops


def get_data_from_vcf(vcf_file, pops_dic):
    """
    Uses the vcf module to parse and retrieve data from a VCF file. The data
    is returned as a list of lists
    """

    storage = [["Reference", "Alternative", "Allele1"] + list(pops_dic.keys()) + ["\tAllele2"] + list(pops_dic.keys()) + ["\tChr", "Pos"]]

    # Get vcf handle
    vcf_handle = vcf.Reader(open(vcf_file))

    for record in vcf_handle:

        # Ignore if more than two alleles
        if len(record.alleles) > 2:
            continue

    	# Get major and minor alleles for first two columns
        ref_allele = record.alleles[0]
        alt_allele = str(record.alleles[1])

		# Get allele counts for populations
        pop_ref_counts = dict((x, 0) for x in pops_dic)
        pop_alt_counts = dict((x, 0) for x in pops_dic)

        for pop, samples in pops_dic.items():
            for ind in samples:
                # Count ref alleles
                try:
                    pop_ref_counts[pop] += \
                        record.genotype(ind).gt_bases.count(ref_allele)
                except AttributeError:
                    pass
                try:
                    pop_alt_counts[pop] += \
                        record.genotype(ind).gt_bases.count(alt_allele)
                except AttributeError:
                    pass

        # Create line record
        ref_snippet = "-{}-".format(ref_allele)
        alt_snippet = "-{}-".format(alt_allele)

        line_record = [ref_snippet, alt_snippet, ref_allele] + list([str(x) for x in pop_ref_counts.values()]) + [alt_allele] + list([str(x) for x in pop_alt_counts.values()]) + [str(record.CHROM), str(record.POS)]

        storage.append(line_record)

    return storage


def write_snp_file(records, output_file):
    """
    Writes the records into a SNP file format
    """

    fh = open(output_file, "w")

    for rec in records:
        fh.write("\t".join(rec) + "\n")

    fh.close()

def main():
    # Get arguments
    vcf_file = arg.vcf_file
    output_file = arg.output_file
    pop_file = arg.population_file

    pop_dic = parse_populations(pop_file)

    records = get_data_from_vcf(vcf_file, pop_dic)

    write_snp_file(records, output_file)


main()
