# Copyright INRA
# author: Renaud VITALIS (2013)
#  
# renaud.vitalis@inra.fr
# 
# Mariana MORA (2023)
#
#fc52541@alunos.fc.ul.pt
# This file is part of SelEstim.
# 
# SelEstim is a computer program whose purpose is to is to detect
# and measure selection from gene frequency data.
# 
# SelEstim is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# SelEstim is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

library(gplots)

plot.delta <- function(file = "summary_delta.out",map = "") {
  
  # 'map' is the name of a file that contains two columns only: the SNP ID and its position (in bp). Be careful to keep the same order of the SNPs as in the data file)
  
  palette(rich.colors(64))
  output <- read.table(file,header = TRUE)
  mean.delta <- output$mean
  relative.delta <- (mean.delta / max(mean.delta))
  n.snps <- length(mean.delta)
  if (map != "") {
    if(!file.exists(map)) {
      stop(paste("The file ",map," does not exist...",sep = ""))
    }
    physical.map <- read.table(map,header = FALSE)
    position <- physical.map[,2] / 1e6
    sorted <- order(position)
    plot(position[sorted],mean.delta[sorted],type = "n",xlab = "Position (Mb)",ylab = expression(Locus-specific~selection~coefficient~delta[j]),ylim = c(0,max(mean.delta) + 10))
    segments(position[sorted],0,position[sorted],mean.delta[sorted],col = 1 + 63 * relative.delta[sorted],lwd = relative.delta[sorted])
  }
  else {
    plot(seq(1,n.snps),mean.delta,type = "n",xlab = "Markers",ylab = expression(Locus-specific~selection~coefficient~delta[j]),ylim = c(0,max(mean.delta) + 10))
    segments(seq(1,n.snps),0,seq(1,n.snps),mean.delta,col = 1 + 63 * relative.delta,lwd = relative.delta)
  }
}

plot.kld <- function(file = "summary_delta.out",map = "",calibration_file = "calibration/summary_delta.out",limit,window.size,n.markers) {
  
  # 'map' is the name of a file that contains two columns only: the SNP ID and its position (in bp). Be careful to keep the same order of the SNPs as in the data file)
  # 'limit' (optional) is used to compute the threshold value of the empirical distribution of the KLD used to calibrate the KLD;
  #  e.g., if you chose limit = 0.01, then the 99\%-quantile of the KLD distribution from the pod analysis will be used as a decision criterion to discriminate between selection and neutrality.
  # 'window.size (optional) is the size of the sliding window, which may be used to determine the regions that contain at least n.markers with KLD above the threshold determined by the 'limit' argument
  # 'n.markers' (optional) is uded to determine the regions that contain at least 'n.markers' with KLD above the threshold determined by the 'limit' argument
  
  output <- read.table(file,header = TRUE)
  kld <- output$KLD
  n.snps <- length(kld)
  if (map != "") {
    if (!file.exists(map)) {
      stop(paste("The file ",map," does not exist...",sep = ""))
    }
    physical.map <- read.table(map,header = F)
    position <- physical.map[,2]
    sorted <- order(position)
    sorted.position <- position[sorted]
    sorted.position.mb <- sorted.position / 1e6
    plot(sorted.position.mb,kld[sorted],cex = 0.25,xlab = "Position (Mb)",ylab = "Kullback-Leibler divergence (KLD)",pch = 16,col = "grey",type = "n")
    points(sorted.position.mb,kld[sorted],cex = 0.25,pch = 16,col = "grey")
    if (!missing(window.size) & !missing(n.markers) & !missing(limit)) {
      if (!file.exists(calibration_file)) {
        stop(paste("The file ",calibration_file," does not exist...",sep = ""))
      }
      calibration <- read.table(calibration_file,header = TRUE)
      kld.calibration <- calibration$KLD
      threshold <- quantile(kld.calibration,(1 - limit))
      outstanding.region <- vector("numeric",n.snps)
      for (i in 1:n.snps) {
        window <- abs(sorted.position[i] - sorted.position) <= (window.size / 2)
        outstanding.region[i] <- length(kld[window][kld[window] >= threshold])
      }
      within <- (outstanding.region >= n.markers & kld >= threshold)
      if (length(within[within]) > 0) {
        segments(sorted.position.mb[within],0, sorted.position.mb[within],max(kld),col = rich.colors(1,alpha = 0.2),lwd = 2)
      }
      points(sorted.position.mb[within][kld[within] >= threshold],kld[within][kld[within] >= threshold],cex = 0.25,col = "black",pch = 8)
    }
    else {
      if (!missing(limit)) {
        if (!file.exists(calibration_file)) {
          stop(paste("The file ",calibration_file," does not exist...",sep = ""))
        }
        calibration <- read.table(calibration_file,header = TRUE)
        kld.calibration <- calibration$KLD
        threshold <- quantile(kld.calibration,(1 - limit))
        #				points(position.mb[sorted][kld[sorted] >= threshold],kld[sorted][kld[sorted] >= threshold],cex = 0.25,col = "black",pch = 8)
        points(sorted.position.mb[kld[sorted] >= threshold],kld[sorted][kld[sorted] >= threshold],cex = 0.25,col = "black",pch = 8)
      }
    }
  }
  else {
    plot(seq(1,n.snps),kld,cex = 0.25,xlab = "Markers",ylab = "Kullback-Leibler divergence (KLD)",pch = 16,col = "grey",type = "n")
    points(seq(1,n.snps),kld,cex = 0.25,pch = 16,col = "grey")
    if (!missing(limit)) {
      if (!file.exists(calibration_file)) {
        stop(paste("The file ",calibration_file," does not exist...",sep = ""))
      }
      calibration <- read.table(calibration_file,header = TRUE)
      kld.calibration <- calibration$KLD
      threshold <- quantile(kld.calibration,(1 - limit))
      points(seq(1,n.snps)[kld >= threshold],kld[kld >= threshold],cex = 0.25,col = "black",pch = 8)
    }
  }
}

randomize.reference.allele <- function (infile = "",outfile = "",pool = FALSE) {
  
  # 'infile' is the original dataset
  # 'outfile' is the dataset where reference alleles are chosen randomly
  # 'pool' is an indicator variable that says whether the data consist in allele counts or reads (pooled data)
  # the 'reference.allele' output file contains, for each locus, the reference allele chosen from the original data
  
  if (infile != "") {
    if(!file.exists(infile)) {
      stop(paste("The file ",infile," does not exist...",sep = ""))
    }
  }
  else {
    stop(paste("\n\tThe argument \"infile\" is missing, with no default value",sep = ""))
  }
  if (outfile == "") {
    stop(paste("\n\tThe argument \"outfile\" is missing, with no default value",sep = ""))
  }
  skip.lines <- 2
  if (pool) {
    skip.lines <- skip.lines + 1
  }
  orig.data <- read.table(infile,skip = skip.lines)
  dummy <- scan(infile,nmax = 2)
  number.of.populations <- dummy[1]
  number.of.loci <- dummy[2]
  if (number.of.populations != (ncol(orig.data) / 2)) {
    stop(paste("\tProblem reading file \"infile\": perhaps the \"pool\" argument is misspecified",sep = ""))
  }
  if (number.of.loci != nrow(orig.data)) {
    stop(paste("\tProblem reading file \"infile\": perhaps the \"pool\" argument is misspecified",sep = ""))
  }
  new.data <- matrix(nrow = nrow(orig.data),ncol = ncol(orig.data))
  flip <- vector(mode = "numeric",length = ncol(orig.data))
  cpt <- 0
  for (i in 1: number.of.populations) {
    flip[c(1,2) + cpt] <- c(2,1) + cpt
    cpt <- cpt + 2
  }
  mask <- sample(c(TRUE,FALSE),size = number.of.loci,replace = TRUE)
  same <- seq(1,number.of.loci)[mask]
  anti <- seq(1,number.of.loci)[!mask]
  new.data[same,] <- as.matrix(orig.data[same,])
  new.data[anti,] <- as.matrix(orig.data[anti,flip])
  write.table(number.of.populations,file = outfile,row.names = FALSE,col.names = FALSE,sep = '\t')
  write.table(number.of.loci,file = outfile,row.names = FALSE,col.names = FALSE,sep = '\t',append = TRUE)
  if (pool) {
    sample.size <- read.table(infile,nrows = 1,skip = 2)
    write.table(sample.size,file = outfile,row.names = FALSE,col.names = FALSE,sep = '\t',append = TRUE)
  }
  write.table(new.data,file = outfile,row.names = FALSE,col.names = FALSE,sep = '\t',append = TRUE)
  list.alleles <- 2 - mask
  write.table(cbind(seq(1,number.of.loci),list.alleles),file = "reference.allele",quote = FALSE,row.names = FALSE,col.names = c("locus","allele"),sep = '\t')
}

compute.F_ST <- function(infile = "",pool = FALSE) {
  
  # 'infile' is the original dataset (in SelEstim format)
  # 'pool' is an indicator variable that says whether the data consist in allele counts or reads (pooled data)
  # The F_ST for pool-seq data is compuyed following Hivert et al. (in prep.)
  # The compute.F_ST function results in a list of two elements: F_ST (a vector of locus-specific F_ST estimates) and F_ST_multilocus (the multilocus estimate of F_ST)
  
  if (infile != "") {
    if(!file.exists(infile)) {
      stop(paste("The file ",infile," does not exist...",sep = ""))
    }
  }
  else {
    stop(paste("\n\tThe argument \"infile\" is missing, with no default value",sep = ""))
  }
  if (!pool) {
    skip.lines <- 2
    counts <- read.table(infile,skip = skip.lines)
    dummy <- scan(infile,nmax = 2)
    number.of.populations <- dummy[1]
    number.of.loci <- dummy[2]
    if (number.of.populations != (ncol(counts) / 2)) {
      stop(paste("\tProblem reading file \"infile\": perhaps the \"pool\" argument is misspecified",sep = ""))
    }
    if (number.of.loci != nrow(counts)) {
      stop(paste("\tProblem reading file \"infile\": perhaps the \"pool\" argument is misspecified",sep = ""))
    }
    r <- ncol(counts) / 2
    l <- seq(1,(2 * r),2)
    ss <- counts[,l] + counts[,(l + 1)]	
    ss2 <- rowSums((counts[,l] + counts[,(l + 1)])^2)
    n <- rowSums(ss)
    n_c <- (n - ss2 / n) / (r - 1.0)
    p <- counts[,l] / ss
    q <- counts[,(l + 1)] / ss
    pbar <- rowSums(counts[,l]) / rowSums(ss)
    qbar <- rowSums(counts[,(l + 1)]) / rowSums(ss) 
    SSI <- rowSums(ss * (p - p^2) + ss * (q - q^2))
    SSP <- rowSums(ss * (p - pbar)^2 + ss * (q - qbar)^2)		
    MSI <- SSI / (n - r)
    MSP <- SSP / (r - 1.0)
  }
  else {
    skip.lines <- 3	
    reads <- read.table(infile,skip = skip.lines)
    dummy <- scan(infile,nmax = 2)
    number.of.populations <- dummy[1]
    number.of.loci <- dummy[2]
    if (number.of.populations != (ncol(reads) / 2)) {
      stop(paste("\tProblem reading file \"infile\": perhaps the \"pool\" argument is misspecified",sep = ""))
    }
    if (number.of.loci != nrow(reads)) {
      stop(paste("\tProblem reading file \"infile\": perhaps the \"pool\" argument is misspecified",sep = ""))
    }
    n_i <- as.vector(read.table(infile,skip = (skip.lines - 1),nrows = 1),mode = "numeric")
    if (length(n_i) != (ncol(reads) / 2)) {
      stop(paste("\tProblem reading file \"infile\": the pool sizes are misspecified",sep = ""))
    }
    nbr_loci <- nrow(reads)
    n_d <- ncol(reads) / 2
    l <- seq(1,(2 * n_d),2)
    mtrx.n_i <- matrix(n_i,nrow = nbr_loci,ncol = n_d,byrow = TRUE)
    R_1_i <- reads[,l] + reads[,(l + 1)]
    R_1 <- rowSums(R_1_i)
    R_2 <- rowSums(R_1_i^2)
    C_1 <- rowSums(R_1_i / mtrx.n_i + (mtrx.n_i - 1) / mtrx.n_i)
    C_1.star <- rowSums(R_1_i * (R_1_i / mtrx.n_i + (mtrx.n_i - 1) / mtrx.n_i)) / R_1
    n_c <- (R_1 - R_2 / R_1) / (C_1 - C_1.star)
    SSI <- rowSums(reads[,l] - reads[,l]^2 / R_1_i + reads[,(l + 1)] - reads[,(l + 1)]^2 / R_1_i)	
    SSP <- rowSums(R_1_i * ((reads[,l] / R_1_i) - (rowSums(reads[,l]) / R_1))^2 + R_1_i * ((reads[,(l + 1)] / R_1_i) - (rowSums(reads[,(l + 1)]) / R_1))^2)
    MSI <- SSI / (R_1 - C_1)
    MSP <- SSP / (C_1 - C_1.star)
  }
  F_ST <- (MSP - MSI)  / (MSP + (n_c - 1) * MSI)
  F_ST_multilocus <- sum(MSP - MSI)  / sum(MSP + (n_c - 1) * MSI)
  rslt <- list(F_ST = F_ST,F_ST_multilocus = F_ST_multilocus)
  return(rslt)
}

#SelEstim_plots plots the results of a selection detection analysis performed with the program SellEstim using
#the summary_delta.out obtained from the main run, the summary_delta.out file obtained from the calibration run, the summary_sigma.out.
#from the main run
#optionally a data map
#Usage:
#summary_delta.out [main run] summary_delta.out [calibration run] [summary_sigma.out [main run] data map [optional]
#Note: uncomment the lines using the data map as arguments to use it.

#Variable to allow using only the files as inputs
args <- commandArgs(trailingOnly = TRUE)

# Get arguments
delta.out <- args[1]
delta.out_calibration <- args[2]
sigma.out <- args[3]
output_namedelta <- args[4]
output_namekld <- args[5]

#data.map <- args[6]
#Plot delta values
pdf(paste(output_namedelta, ".pdf", sep=""))
plot.delta(file = delta.out,)
# map = data.map)
dev.off()
#Plot kld values
pdf(paste(output_namekld, ".pdf", sep=""))
plot.kld(file = delta.out,)
dev.off()
# map = data.map)
#Update kld values plot to highlight loci under selection
pdf("~/Desktop/selestim3.pdf")
plot.kld(file = delta.out,
         #  map = data.map
         calibration_file =
          delta.out_calibration,
         limit = 0.001)
rslt <- read.table(delta.out,
                   header = TRUE)
top.snp <- which(rslt$KLD == max(rslt$KLD))
top.snp
dev.off()
#Add line to divide manhatan plot in chromossomes
pdf("~/Desktop/selestim4.pdf")
plot.kld(file = delta.out,
         #  map = data.map
         calibration_file =
           delta.out_calibration,
         limit = 0.01)
abline(v = 4.867859,lty = 2)
rslt$mean[top.snp]
sigma <- read.table(sigma.out, header = TRUE)
sigma$mean[which(sigma$locus == top.snp)]
dev.off()
