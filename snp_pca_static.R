#!/usr/bin/Rscript
# snp_pca.R performs a PCA using the SNPRelate R package using a VCF file
# and an option populations files

# Altered from the original: 
#   displays the percentage of variance explained by each principal component in the axis labels
#   slightly altered the pch and cls variables

# Usage:
# snp_pca.R vcf_file output_file_name popupations_file[optional]

library("SNPRelate")
library("RColorBrewer")

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
  pop_code <- read.table(pops_file, sep="\t")
  sorted_pops <- pop_code$V2[order(match(pop_code$V1, sample.id))]
  tab <- data.frame(sample.id = pca$sample.id,
    pop = factor(sorted_pops),
    EV1 = pca$eigenvect[,1],
    EV2 = pca$eigenvect[,2],
    stringsAsFactors=F)
  cls <-rep(brewer.pal(n = 5, name = "Set1"), times=5)
  pch_v <- rep(c(16, 15, 17, 18), each=5)
  save(tab, file=paste(output_name, ".Rdata", sep=""))
  # par(mar =  c(5, 4, 4, 6) + 1.8)
  
  plot.new()
  leg = legend(0, 0, legend=levels(tab$pop), bty='n', pch=pch_v, plot=FALSE,
               col=cls[0:length(tab$pop)])
  # calculate right margin width in ndc
  leg_wid <- grconvertX(leg$rect$w, to='ndc') - grconvertX(0, to='ndc')
  
  par(omd=c(0, 1-leg_wid, 0, 1))
  plot1 <- plot(tab$EV1, tab$EV2, col=cls[as.integer(tab$pop)], xlab=paste("PC 1: ", print(round(pc.percent[1],2)),"%", sep = ""),
    ylab=paste("PC 2: ", print(round(pc.percent[2],2)),"%", sep = ""), pch=pch_v[as.numeric(tab$pop)], solid=.2, cex=1.2,
    clab=1, leg=T, bg="white")
  
  text(tab$EV1 - 1, tab$EV2, labels=tab$sample.id)

  legend(par('usr')[2], par('usr')[4], legend=levels(tab$pop), bty='n', xpd=NA,
         col=cls[0:length(tab$pop)], pch=pch_v)
  
  print(cls[as.numeric(sorted_pops)])
  
} else {
  tab <- data.frame(sample.id = pca$sample.id,
    EV1 = pca$eigenvect[, 1],
    EV2 = pca$eigenvect[, 2],
    stringsAsFactors=F)
  plot(tab$EV1, tab$EV2, xlab="PC 1", ylab="PC 2")
}

# remove temporary gds file
file.remove("temp.gds")

dev.off()