#!/usr/bin/python3

from sys import argv

def Ambiguifier(bases):
	"""Take a list or tuple of bases and returns the corresondig ambiguity."""
	bases = list(bases)
	ambigs = {"AA": "A", "CC": "C", "TT": "T", "GG": "G", "AC": "M", "AG": "R",
	"AT": "W", "CG": "S", "CT": "Y", "GT": "K", "ACG": "V",
	"ACT": "H", "AGT": "D", "CGT": "B", "ACGT": "N",
	"..": "N"}
	
	bases.sort()
	
	return ambigs["".join(bases)]
	
def VCF_parser(vcf_file_name):
	"""Parse a VCF file and return {name:seq}."""
	infile = open(vcf_file_name,'r')
	sequences = {}
	for lines in infile:
		if lines.startswith("##"):
			pass
		elif lines.startswith("#"):
			taxa = lines.split()[9:]
			for i in taxa:
				sequences[i] = ""
		else:
			sequences = VCF_line_parser(lines, sequences, taxa)
			
	return sequences
			

def VCF_line_parser(line, sequences, taxa):
	"""Parse a single data line of a VCF file and adds each individual
	base to the sequences dict, which is returned."""
	line = line.strip().split()
	aleles = ",".join(line[3:5])
	aleles = aleles.split(",")
	bases = line[9:]
	for b, t in zip(bases, taxa):
		b = b.replace("|", "")
		for a in range(len(aleles)):
			b = b.replace(str(a), aleles[a])
		b = b.replace("./.", "..")
		b = Ambiguifier(b)
		sequences[t] += b
		
	return sequences


def Phylip_writer(sequences):
	"""Take a dict and write a phylip file (no header so far)."""
	for names in sequences:
		print(names + "\t", end="")
		print(sequences[names], end="\n")
	

if __name__ == "__main__":
	# Usage: python3 VCF2phy.py file.vcf
	sequences = VCF_parser(argv[1])
	Phylip_writer(sequences)
