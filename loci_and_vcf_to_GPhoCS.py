#!/usr/bin/python3
from sys import argv
vcf= open(argv[1],"r") #input = o output .vcf do vcftools
raw= open (argv[2], "r") #input= ficheiro.loci do pyrad
poGPhosCS= open (argv[3], "w") #output pronto po GPhoCS

d=[]

for n in vcf:
	if n.startswith("#") ==False:
		s=n.split()
		d.append(s[0]) #saca o numero do locus do vcf
		

d= sorted(list(set(d))) #vamos ficar com a lista de loci de interesse sem as repeticoes devido aos SNPs em diferentes posicoes desse locus

if d[0] == "0": 
	estouinteressado=1	#variável temporária para selecionar 	
else:
	estouinteressado=0   	#isto e porque o vcf nesta versão do pyRad começa a contar do 0, e sem isto nao apanharia as seq. do primeiro locus.

c=0 # numero do loci
sequences=[]
poGPhosCS.write(str(len(d)) + "\n\n")
for lines in raw:
	if estouinteressado==1 and lines.startswith("//")==False:
		lines=lines.replace("-","N") #O GPhoCS lê missing data como N, so aceita IUPAC e N no other characters
		sequences.append(lines) #lines e ja o nome do individuo e a seq. correspondente
	elif lines.startswith("//"): #conta os locis e faz o increase do c
		if sequences!=[]:
			name=str(c) #nome do locus, que vem do vcf
			numbseq=str(len(sequences)) #conta o numero de inviduos(sequencias) daquele locus
			seqlen=str(len(sequences[0].split()[1])) #quantifica o comprimento das sequencias
			poGPhosCS.write(name + " " + numbseq + " " + seqlen + "\n")
			for seqs in sequences:
				poGPhosCS.write(seqs)
			poGPhosCS.write("\n")
		c+=1
		estouinteressado=0
		sequences=[]
		if str(c) in d: #ver se o loci está na lista do vcf 
			estouinteressado=1 #se o loci estiver no vcf então append no prontopoGPhoCS 


vcf.close()
raw.close()
poGPhosCS.close()

