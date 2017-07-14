#!/usr/bin/python3
# Copyright 2015 Diogo N Silva <o.diogosilva@gmail.com>
# This file is part of VCP2phy.
# VCP2phy is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# VCP2phy is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with VCP2phy. If not, see <http://www.gnu.org/licenses/>.

import argparse
from collections import defaultdict, OrderedDict

parser = argparse.ArgumentParser(description="Converter of VCF into phylip "
											 "format")
parser.add_argument("-vcf", dest="vcf_file", help="VCF file")
parser.add_argument("-loci", dest="loci_file", help="If provided, the output "
					"file will contain the full sequences contained in the "
					" loci file.")
parser.add_argument("-o", dest="output_file", help="Name of output file")

arg = parser.parse_args()

iupac = {"AG": "R", "CT": "Y", "CG": "S", "AT": "W", "GT": "K", "AC": "M",
         "CGT": "B", "AGT": "D", "ACT": "H", "ACG": "V", "ACGT": "N", "AA": "A",
		 "CC": "C", "GG": "G", "TT":"T"}

def parse_vcf(vcf_file):
	"""
	Parses the VCF file and returns a dictionary with the loci (chromossomes)
	as keys and a list of positions as values.
	"""

	loci = defaultdict(list)

	vcf_handle = open(vcf_file)

	for line in vcf_handle:

		if line.startswith("##"):
			pass

		elif line.startswith("#CHROM"):
			taxa_list = line.strip().split()[9:]

		elif line.strip() != "":
			fields = line.strip().split()
			loci[int(fields[0].split("_")[-1])].append(int(fields[1]) - 1)

	vcf_handle.close()

	return loci, taxa_list


def convert_genotype(genotype_string, variants):
	"""
	Converts a genotype string from the vcf ("0|0") into an actual genotype
	("A"), provided the variants for that particular site
	"""

	# Remove separators
	genotype_string =  genotype_string.replace("|", "").replace("/", "").upper()

	if genotype_string == "..":
		return "N"

	# Convert to actual bases
	actual_genotype = "".join(sorted([variants[int(x)] for x in
									  genotype_string]))

	return iupac[actual_genotype]



def parse_vcf_variants(vcf_file):
	"""
	Parses the VCF but retains the variation information. This is done when
	converting a VCF directly to phylip without using a .loci file. The
	resulting phylip file will contain only variable sites, though.
	"""

	vcf_handle = open(vcf_file)
	# Will store the sequence data for each taxon in the VCF


	for line in vcf_handle:

		if line.startswith("##"):
			pass

		elif line.startswith("#CHROM"):
			# List of taxa that will retain their order in the VCF
			taxa_list = line.strip().split()[9:]
			# Will store the sequence data for each taxon in the VCF
			seq_data = dict((tx,[]) for tx in taxa_list)

		# Skip empty lines and read the VCF data
		elif line.strip() != "":
			fields = line.strip().split()
			# e.g. ("A","T")
			variants = [fields[3]] + fields[4].split(",")
			# e.g. ["0|0", ".|.", "1|1"]
			genotypes = [x.split(":")[0] for x in fields[9:]]

			# Add information to seq_data for each data point
			for p,gen in enumerate(genotypes):

				# get the genotype
				genotype = convert_genotype(gen, variants)

				# Get taxon
				tx = taxa_list[p]

				seq_data[tx].append(genotype)

	return seq_data


def mask_alignment(aln, var_positions):
	"""
	Parses an alignment in an OrderedDict format {taxon: seq} and returns a
	similar alignment OrderedDict with the variable positions not present in the
	var_positions list masked
	"""

	masked_dic = OrderedDict((tx, []) for tx in aln)

	for p, column in enumerate(zip(*aln.values())):

		gapless_column = len(set([x for x in column if x.lower() not in ["-", "n"]]))

		# If this represents a real variable position or a column with no
		# variation (excluding missing data), save sequences as they
		# are
		if p in var_positions or (column not in var_positions and
							      gapless_column == 1):

			# This appends the existing nucleotide of each taxon to its
			# corresponding taxon in masked_dic
			list(map(lambda char, tx: masked_dic[tx].append(char),
			                column, list(aln.keys())))

		elif p not in var_positions and gapless_column > 1:

			# This appends an "n" to each taxon in masked_dic
			list(map(lambda tx: masked_dic[tx].append("n"), list(aln.keys())))

	masked_dic = OrderedDict((tx, "".join(seq)) for tx, seq in masked_dic.items())

	return masked_dic


def parse_loci(loci_file, vcf_loci, taxa_list):
	"""
	Writes a phylip file from a loci file, according to the positions in the
	VCF. Variable positions that are not present in the VCF are soft masked
	"""

	fh = open(loci_file)

	current_aln = OrderedDict()

	# Keeps track of the sequence lenght for each loci. For the partitions file.
	seq_lens = []

	# Using list values instead of strings makes the append process much faster
	total_aln = dict((x, []) for x in taxa_list)

	for line in fh:

		if line.startswith("//"):

			# Get locus number
			locus = int(line.split("|")[1]) + 1

			print("\rParsing alignment {}".format(locus), end="")

			if locus in vcf_loci:

				if current_aln:
					masked_aln = mask_alignment(current_aln, vcf_loci[locus])
					current_len = len(list(masked_aln.values())[0])

					# Append masking result to total alignment
					for tx in taxa_list:
						try:
							total_aln[tx].append(masked_aln[tx])
						except KeyError:
							total_aln[tx].append("n" * current_len)

					seq_lens.append(current_len)

			current_aln = OrderedDict()

		# Loci present in VCF. Begin alignment parsing routine
		elif line.strip() != "":
		# Get alignment

			fields = line.strip().split()
			# Removes ">" from taxon name
			taxon = fields[0]
			current_aln[taxon] = fields[1]

	fh.close()

	# Finalize alignment object by converting list values into strings
	return dict((tx, "".join(seq)) for tx, seq in total_aln.items()), seq_lens


def write_to_phy(aln_dict, seq_lens, output_file):
	"""
	Writes a alignment dictionary {taxon: seq} to a phylip file
	"""

	# Write main phylip file
	fh = open(output_file + ".phy", "w")

	fh.write("{} {}\n".format(len(aln_dict),
	                          len(list(aln_dict.values())[0])))

	for tx, seq in aln_dict.items():
		fh.write("{}\t{}\n".format(tx, seq))

	fh.close()

	# Write partitions file
	fh = open(output_file + ".partitions", "w")
	locus_n = 1
	c = 1

	for l in seq_lens:
		fh.write("DNA, p{}={}-{}\n".format(locus_n, c, c + l - 1))
		c += l
		locus_n += 1

	fh.close()


def write_vcf_variants(seq_data, output_file):
	"""
	Writes the variable sites of the VCF file directly to phylip format
	"""

	fh = open(output_file + ".phy", "w")

	# Get number taxa
	ntaxa = len(seq_data)

	# Get alignment length
	nlen = len(list(seq_data.values())[0])

	# Write header
	fh.write("{}\t{}\n".format(ntaxa, nlen))

	# Write alignment
	for tx, seq_list in seq_data.items():
		fh.write("{}\t{}\n".format(tx, "".join(seq_list)))

	fh.close()


def main():

	# Arguments
	vcf_file = arg.vcf_file
	loci_file = arg.loci_file
	output_file = arg.output_file

	if not loci_file:
		print("Loci file not provided. Converting only variable sites in VCF.")
		seq_data = parse_vcf_variants(vcf_file)
		write_vcf_variants(seq_data, output_file)
		return

	# Parse VCF
	loci_storage, taxa_list = parse_vcf(vcf_file)

	# Parse loci
	aln_obj, seq_lens = parse_loci(loci_file, loci_storage, taxa_list)

	# Write phylip file
	write_to_phy(aln_obj, seq_lens, output_file)

if __name__ == "__main__":
	main()
