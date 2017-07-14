#!/usr/bin/python

# vcf2fas.py converts a vcf file with multiple samples into a fasta sequence for ech locus based on SAM files for each sample in the vcf file

import argparse


parser = argparse.ArgumentParser(description="Converts VCF files into Fasta files from SAM")

parser.add_argument("-vcf", dest="vcf_infile", help="Provide the VCF file.", required=True)

arg = parser.parse_args()

def get_sam_sequences(fname, chrom, pos):
    """
    Gets a list of sequences from a SAM file according to a chromosome name and position. It also assumes that reads are 90bp long
    """

    ref_seqs = []

    sam_fh = open(fname)

    saved = False

    for line in sam_fh:
        if line.startswith("@"):
            pass
        else:
            fields = line.split()
            ref_chrom = fields[2]
            ref_pos = fields[3]
            seq = fields[0]

        if ref_chrom == chrom:
            save = True
            
            #Check if sequence overlaps position
            if pos in range(ref_pos, ref_pos + len(seq)):
                ref_seqs.append(seq)
            
        elif ref_chrom != chrom and not saved:
            pass

        elif ref_chrom != chrom and saved:
            return ref_seqs

def vcf2fasta(vcf_file, mode):
    """
    This checks the coverage of each SNP in the VCF file according to the SAM files
    """

    vcf_fh = open(vcf_file)

    line = next(vcf_fh)

    if mode == "check_coverage":
        bad_snps = open("bad_coverage_loci.csv", "w")
        bad_counts = 0
        total_counts = 0

    taxa_pos = {}

    # Skip VCF header
    for line in vcf_file:
        if line.startswith("#"):
            pass
        
        elif line.startswith("CHOM"):

            fields = line.split()

            taxa_names = [x.strip() for x in fields[9:]]

            for sample in taxa_names:
                taxa_pos[sample] = fields.index(sample)

        else:
            # Iterate over every SNP
            fields = line.split()

            chrom = fields[0]
            pos = fields[1]

            genotypes = [x.strip().split(":")[0] for x in fields[9:]]
            coverage = [x.strip().split(":")[-1] for x in fields[9:]]

            ref_coverage = []

            # Check SNP coverage for each sample
            for taxon, gen in zip(taxa_names, genotypes):
                # Only do this for present genotypes
                if gen != "./.":
                    sequences = get_sam_sequences("{}.sam".format(taxon), chrom, pos)
                    ref_coverage.append(len(sequences))
                    if mode == "check_coverage":
                        total_counts += 1
                else:
                    if mode == "check_coverage":
                        ref_coverage.append(None)

                if mode == "convert":


            if mode == "check_coverage":
                # Check coverage
                for ref_cov, cov, taxon in zip(ref_coverage, coverage, taxa_names):
                    if ref_cov:
                        if ref_cov != cov:
                            bad_snps.write("{}; {}; {}\n".format(ref_cov, cov, taxon))
                            bad_counts += 1

    if mode == "check_coverage":
        print("There were %s bad SNP coverage counts (%s%%)" % (bad_counts, round(float(bad_counts)/float(total_counts), 2)))

def main():
    # Arguments
    vcf_file = arg.vcf_infile


main()
