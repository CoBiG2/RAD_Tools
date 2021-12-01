#!/usr/bin/Rscript
# snp_pca.R performs a PCA using the SNPRelate R package using a VCF file
# and an option populations files

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
  cls <-rep(brewer.pal(n = 12, name = "Set3"), times=5)
  pch_v <- rep(c(16, 15, 17, 18), each=12)
  save(tab, file=paste(output_name, ".Rdata", sep=""))
  # par(mar =  c(5, 4, 4, 6) + 1.8)
  
  plot.new()
  leg = legend(0, 0, legend=levels(tab$pop), bty='n', pch=pch_v, plot=FALSE,
               col=cls[0:length(tab$pop)])
  # calculate right margin width in ndc
  leg_wid <- grconvertX(leg$rect$w, to='ndc') - grconvertX(0, to='ndc')
  
  par(omd=c(0, 1-leg_wid, 0, 1))
  plot(tab$EV1, tab$EV2, col=cls[as.integer(tab$pop)], xlab="PC 1",
    ylab="PC 2", pch=pch_v[as.numeric(tab$pop)], solid=.2, cex=1,
    clab=1, leg=T, bg="white")

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

##create matrix with pairwise euclidean distance
mat <- dist(tab, method="euclidean", diag = F, upper =F, p=2)
distance_matrix <- as.matrix(dist(mat))
write.table(distance_matrix, file = "pairwise_distance_matrix.txt", sep = "\t", row.names = TRUE, col.names = NA)

#list first 6 eigen values
head_eigen <- as.data.frame(head(round(pc.percent, 2))) 
head_eigen$index <- c("1","2","3","4","5","6")
colnames(head_eigen) <- c("varexp","pcindex")

#combined plot - paiwise pca + eigen values plot
library(pryr)
library(grid)
png(paste("pca_combined_plot", ".png", sep=""))
  #plot pairwise pca using the first four vectors
p1 %<a-% pairs(pca$eigenvect[,1:4], labels = lbls, col=tab$pop, 
               lower.panel = NULL,
               oma=c(4,4,6,10), pch=16)
  #plot the first 6 eigen values
p2 <-ggplot(head_eigen, aes(x=pcindex, y = varexp), group=1) +
  geom_bar(stat = "identity", width=0.6, fill="steelblue") +
  theme_bw() +
  geom_point(color="darkblue")+
  geom_text(aes(label=varexp), vjust=-0.3, size=4)+
  geom_line(color="darkblue" ,group=1)+
  ylim(0,5)+
  ylab("Explained variance (%)") + xlab("Principle Component Number")
  #two plots in the same figure
par(xpd=TRUE)
vp <- viewport(width = 0.4, height = 0.3, x = 0.23, y = 0.2)
print(p1)
print(p2, vp=vp)
legend(0.9,0.7, fill=unique(tab$pop), legend = c(levels(tab$pop)), cex = 0.6)
dev.off()
