#!/usr/bin/Rscript

library(Demerelate)

args <- commandArgs(trailingOnly = TRUE)

# Read input file
input_file <- read.csv(args[1], sep=" ", header=TRUE)

# Transform into data frame
input_data <- as.data.frame(input_file)

# Fis calculations
Fis.calc(input_data, 1000, args[2], "object", args[3], "Fis")