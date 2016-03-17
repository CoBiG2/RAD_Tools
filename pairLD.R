#!/usr/bin/Rscript

# pairLD.R calculates pair-wise SNP linkage disequilibrium statistics using a
# a ML approach and outputs several statistics and heatmap plots

# Usage:
# pairLD.R vcf_input_file output_prefix

suppressMessages(library("pegas"))
suppressMessages(library("fields"))
suppressMessages(library("genetics"))
suppressMessages(library("reshape2"))
suppressMessages(library("ggplot2"))

args <- commandArgs(trailingOnly = TRUE)

# Get arguments
vcf_file = args[1]
output_prefix = args[2]

# Read VCF file. Read up to 1M SNPs
print("Reading VCF file")
vcf <- read.vcf(vcf_file, to=1000000, quiet=T)
vcf_info <- VCFloci(vcf_file)

# Replace missing data with NA
vcf[vcf=="./."] <- NA

# Build a data frame of genotypes
print("Building genotypes")
geno <- makeGenotypes(data.frame(vcf))

# Compute pairwise LD statistics
print("Calculating LD")
ldres <- LD(geno)

# Save LD results into file
save(ldres, file=paste(output_prefix, "_ldres.RData", sep=""))

print("Adjusting P-values")
# Create a new ldres slot with FDR corrected p-values
ldres$"Q-values" <- p.adjust(ldres$"P-value", method="fdr")
# Convert qvalues to matrix
ldres$"Q-values" <- t(matrix(ldres$"Q-values", nrow=dim(ldres$"P-value")[1],
                             byrow=T))
print("Plotting")
# Generate heatmap plots for D', R^2, p-value and q-values
# For D'
pdf(paste(output_prefix, "_Dc.pdf", sep=""))
image.plot(ldres$"D'", legend.lab="D'")
dev.off()
# For R^2
pdf(paste(output_prefix, "_R2.pdf", sep=""))
image.plot(ldres$"R^2", legend.lab="R^2")
dev.off()
# For p-value
pdf(paste(output_prefix, "_pval.pdf", sep=""))
image.plot(ldres$"P-value", legend.lab="P-value")
dev.off()
# For q-values
pdf(paste(output_prefix, "_qval.pdf", sep=""))
image.plot(ldres$"Q-values")
dev.off()

# Generate histograms with D' and R2 values
Dc_hist_data <- melt(ldres$"D'")
pdf(paste(output_prefix, "_Dc_hist.pdf", sep=""))
g <- ggplot(Dc_hist_data, aes(x=value))
g + geom_histogram() + xlab("D'")
dev.off()

R2_hist_data <- melt(ldres$"R^2")
pdf(paste(output_prefix, "_R2_hist.pdf", sep=""))
g <- ggplot(R2_hist_data, aes(x=value))
g + geom_histogram() + xlab("R^2")
dev.off()

# Generate plots for D' and R^2, but show only data for pairs with pvalues below
# 0.05
Dc_sig <- ldres$"D'"
Dc_sig[ldres$"Q-values" > 0.05 ] <- NA
pdf(paste(output_prefix, "_Dc_sig.pdf", sep=""))
image.plot(Dc_sig, zlim=c(0,1), legend.lab="D'")
dev.off()

R2_sig <- ldres$"R^2"
R2_sig[ldres$"Q-values" > 0.05 ] <- NA
pdf(paste(output_prefix, "_R2_sig.pdf", sep=""))
image.plot(R2_sig, zlim=c(0,1), legend.lab="R^2")
dev.off()

# Get proportions of statistically significant departures from LD
all_pairwise <- length(which(!is.na(ldres$"Q-values")))
sig_results <- length(which(ldres$"Q-values" <= 0.05))
sig_prop <- sig_results / all_pairwise

# Get mean and sd from D' and R2
Dc_mean = mean(ldres$"D'", na.rm=T)
Dc_stdev = sd(ldres$"D'", na.rm=T)
R2_mean = mean(ldres$"R^2", na.rm=T)
R2_stdev = sd(ldres$"R^2", na.rm=T)

print("Writing significant SNP pairs")
# Get SNP pairs coordinates with significant q-values
snp_coord <- which(ldres$"Q-values" <= 0.05, arr.ind=TRUE)
# Convert SNPs pairs coordinates into the actual chromosome and position
snp_pairs <- t(apply(snp_coord, 1, function(r) c(vcf_info$CHROM[r[1]],
  vcf_info$POS[r[1]], vcf_info$CHROM[r[2]], vcf_info$POS[r[2]], ldres$"R^2"[r[1], r[2]], ldres$"D'"[r[1], r[2]])))

# Write pairs of SNPs with significant LD into csv
table_data <- matrix(c("CHROM1", "POS1", "CHROM2", "POS2", "R^2", "D'"), ncol=6)
table_data <- rbind(table_data, snp_pairs)
write.matrix(table_data, file="Significant_snp_pairs.csv", sep=";")

# Write info to file
write(c(paste("Total number of pair wise comparisons: ", all_pairwise, sep=""),
        paste("Number of significant LD comparisons: ", sig_results, sep=""),
        paste("Porportion of significant LD comparisons: ", sig_prop, sep=""),
        paste("Mean D': ", Dc_mean, sep=""),
        paste("Stdev D': ", Dc_stdev, sep=""),
        paste("Mean R2: ", R2_mean, sep=""),
        paste("Stdev R2: ", R2_stdev, sep="")),
      file=paste(output_prefix, ".log", sep=""))
