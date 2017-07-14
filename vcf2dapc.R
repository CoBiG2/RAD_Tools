#!/usr/bin/Rscript

# vcf2dapc.R performs an interactive DAPC analyses from a VCF input file

# This is not a standalone script, due to interactive queries of some adegenet
# functions. Use the interactive mode to run these commands.

suppressMessages(library("vcfR"))
suppressMessages(library("adegenet"))

args <- commandArgs(trailingOnly = T)

# Get arguments
vcf_file = args[1]
max_cl = args[2]
output_file = args[3]

# Reading VCF file
print("Reading VCF file")
vcf <- read.vcfR(vcf_file)
print("Converting VCF to GenInd format")
vcf_gi <- vcfR2genlight(vcf)

# Determine optimum number of clusters
print("Determining optimum number of clusters")
grp <- find.clusters(vcf_gi, choose.n.clust=F, criterion="diffNgroup")

# Compute DAPC
print("Computing DAPC")
dapc1 <- dapc(vcf_gi, grp$grp)

# Printing plot
print("Generating plot")

pdf(paste(output_file, "_assignments.pdf", sep=""))
assignplot(dapc1)
dev.off()

pdf(paste(output_file, "_dapc.pdf", sep=""))
scatter(dapc1, scree.da=F, bg="white", pch=20, cell=0, cstar=0,
        solid=.4, cex=3, clab=0, leg=T,
        txt.leg=paste("Cluster", 1:length(grp$size)))

# BW version
#scatter(dapc1, scree.da=F, bg="white", pch=c(5,7, 6), cell=0, cstar=0,
#         solid=.4, cex=2, clab=0, leg=T, col="grey", lwd=2)
