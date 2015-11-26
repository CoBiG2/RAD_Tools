#!/usr/bin/python


# pairwise_ld.py. Calculates the D' and r2 values for each SNP pair in a VCF file


import argparse
import itertools
from scipy.stats import chisquare
import numpy
import matplotlib.pyplot as plt


parser = argparse.ArgumentParser(description="Comparison of RAD assemblies "
                                             "using technical replicates")

parser.add_argument("-vcf", dest="vcf_infile", help="Provide VCF "
                    " file.", required=True)

arg = parser.parse_args()


def parse_vcf(vcf_file):
    """
    Parses a VCF file and returns a list of genotypes
    """

    genotypes = []
    x= 0
    vcf_handle = open(vcf_file)

    for line in vcf_handle:

        if line.startswith("#"):
            pass
        else:
            fields = line.strip().split()

            # Check if SNP is biallelic
            ref_snp = fields[3]
            alt_snp = fields[4]

            if len(ref_snp) > 1 or len(alt_snp) > 1:
                continue
            
            # This will add something like ["00", "01", "11", "01", "00"]
            gen = ["".join(x.split(":")[0].split("/")) if x != "./." else ".." for x in fields[9:]]
            genotypes.append((x, gen))
            x += 1

    return genotypes


def compute_pairwise_ld(genotype_list):
    """
    Computes a triangular matrix for D' and r2 for each combinaiton in genotype_list
    """

    D_matrix = numpy.empty((len(genotype_list), len(genotype_list)), dtype=float)
    r2_matrix = numpy.empty((len(genotype_list), len(genotype_list)), dtype=float)

    n = 0
    Dvals = []
    r2vals = []

    for gen_l1, gen_l2 in itertools.combinations(genotype_list, 2):

        pos1 = gen_l1[0]
        pos2 = gen_l2[0]

        gen1 = gen_l1[1]
        gen2 = gen_l2[1]

        # Get total number of alleles
        total_alleles = 0
        p_count = 0
        q_count = 0
        for i,j in zip(gen1, gen2):
            if i == ".." or j == "..":
                pass
            else:
                total_alleles += 1

                # Get allele frequencies for both genotypes
                p_count += i.count("0")
                q_count += j.count("0")

        if total_alleles == 0:
            Dr = None
            r2 = None
            D_matrix[pos1][pos2] = Dr
            r2_matrix[pos1][pos2] = r2
            continue

        total_alleles = total_alleles * 2

        p1 = float(p_count) / float(total_alleles)
        p2 = 1. - p1

        q1 = float(q_count) / float(total_alleles)
        q2 = 1. - q1

        # Get observed genotypes
        gen_counts = {"00": 0., "01": 0., "10": 0., "11": 0.}
        for i,j in zip("".join(gen1), "".join(gen2)):
            if i != "." and j != ".":
                gen_counts["{}{}".format(i,j)] += 1.

        # Get expected genotypes
        exp_gen = {}
        exp_gen["00"] = p1 * q1 * total_alleles
        exp_gen["01"] = p1 * q2 * total_alleles
        exp_gen["10"] = p2 * q1 * total_alleles
        exp_gen["11"] = p2 * q2 * total_alleles

        # Compute chi-square
        var = ["00", "01", "10", "11"]
        obs_list = [gen_counts[x] for x in var]
        exp_list = [exp_gen[x] for x in var]

        chi2 = chisquare(obs_list, exp_list, 1)

        # Quantification of LD
        # Frequencies of genotypes
        freq_00 = float(gen_counts["00"] / total_alleles)
        freq_01 = float(gen_counts["01"] / total_alleles)
        freq_10 = float(gen_counts["10"] / total_alleles)
        freq_11 = float(gen_counts["11"] / total_alleles)

        # Compute D
        D = (freq_00*freq_11) - (freq_01*freq_10)

        # Compute r2
        try:
            r2 = D**2 / (q1 * q2 * p1 * p2)
            r2vals.append(r2)
        except ZeroDivisionError:
            r2 = None

        # Compute Dr (D')
        try:
            if D < 0:
                Dr = abs(D / min([p1*q1, p2*q2]))
            else:
                Dr = abs(D / min([p1*q2, p2*q1]))
            Dvals.append(Dr)
        except ZeroDivisionError:
            Dr = None
            n += 1

        D_matrix[pos1][pos2] = Dr
        r2_matrix[pos1][pos2] = r2

    print("Mean D' value of {} +- {}".format(numpy.mean(Dvals), numpy.std(Dvals)))
    print("Mean r2 value of {} +- {}".format(numpy.mean(r2vals), numpy.std(r2vals)))

    # Plot D' and r2
    fix, ax = plt.subplots()

    heat1 = ax.imshow(D_matrix, interpolation="nearest", cmap="jet")
    cbar = plt.colorbar(heat1)
    cbar.set_label("D'")

    plt.savefig("test.png")
    plt.close()

    print(n)

    fig, ax = plt.subplots()

    hist = plt.hist(Dvals, bins=20)

    plt.savefig("hist.png")



def main():
    # Args
    vcf_file = arg.vcf_infile

    storage = parse_vcf(vcf_file)
    compute_pairwise_ld(storage)

main()
