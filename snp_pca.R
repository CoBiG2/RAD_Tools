#!/usr/bin/Rscript
# snp_pca.R performs a PCA using the SNPRelate R package using a VCF file
# and an option populations files

# Usage:
# snp_pca.R vcf_file output_file_name popupations_file[optional]

library("SNPRelate")

args <- commandArgs(trailingOnly = TRUE)

# Get arguments
vcf_file <- args[1]
output_name <- args[2]
pops_file <- args[3]

# Convert VCF to gds
snpgdsVCF2GDS(vcf_file, "temp.gds", method="biallelic.only")

# Open GDS file
genofile <- snpgdsOpen("temp.gds")

# Run PCA
pca <- snpgdsPCA(genofile, num.thread=1, autosome.only=F)
print(pca)

pc.percent<- pca$varprop * 100
print(round(pc.percent, 2))

# Open figure driver
pdf(paste(output_name, ".pdf", sep=""))

# Plots PCA
if (!is.na(pops_file)) {
  sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
  pop_code <- scan(pops_file, what=character())
  tab <- data.frame(sample.id = pca$sample.id,
    pop = factor(pop_code)[match(pca$sample.id, sample.id)],
    EV1 = pca$eigenvect[,1],
    EV2 = pca$eigenvect[,2],
    stringsAsFactors=F)
  save(tab, file=paste(output_name, ".Rdata", sep=""))
  plot(tab$EV1, tab$EV2, col=adjustcolor(as.integer(tab$pop), alpha=.5), xlab="eigenvector 1",
    ylab="eigenvector 2", pch=20, cell=0, cstar=0, solid=.4, cex=3,
     clab=0, leg=T, scree.da=F, bg="white")
  legend("topleft", legend=levels(tab$pop), pch="o", col=1:nlevels(tab$pop))
} else {
  tab <- data.frame(sample.id = pca$sample.id,
    EV1 = pca$eigenvect[, 1],
    EV2 = pca$eigenvect[, 2],
    stringsAsFactors=F)
  plot(tab$EV1, tab$EV2, xlab="eigenvector 1", ylab="eigenvector 2")
}

# remove temporary gds file
file.remove("temp.gds")

dev.off()
