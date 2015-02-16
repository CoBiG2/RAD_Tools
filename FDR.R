#!/usr/bin/Rscript

#Usage: Rscript FDR.R infile.csv outfile.csv

library(stats)

args <- commandArgs(trailingOnly = TRUE)

#The csv file conatins an error (-100.0) which has to be replaced by (0.5).
cs <- read.csv(args[1], sep="\t", stringsAsFactors=FALSE)

pos <- subset(cs, P.Simul.Fst.sample.Fst.>=0.5)
bal <- subset(cs, P.Simul.Fst.sample.Fst.<0.5)

Qbal <- p.adjust(bal$P.Simul.Fst.sample.Fst., method="fdr")

invpos <- 1 - pos$P.Simul.Fst.sample.Fst.
Qpos <- 1 - p.adjust(invpos, method="fdr")

pos[[5]] <- Qpos
names(pos)[5] <- "Q-values"

bal[[5]] <- Qbal
names(bal)[5] <- "Q-values"

All <- rbind(bal, pos)

All <- All[with(All, order(as.numeric(row.names(All)))), ]

write.table(All, file=args[2], row.names=FALSE, sep="\t")
