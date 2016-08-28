#!/usr/bin/Rscript
# Copyright 2015 Francisco Pina Martins <f.pinamartins@gmail.com>
# This file is part of FDR.R.
# FDR.R is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# FDR.R is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with FDR.R. If not, see <http://www.gnu.org/licenses/>.


#Usage: Rscript FDR.R infile.csv outfile.csv

library(stats)

args <- commandArgs(trailingOnly = TRUE)

#The locilist.csv file conatins an error at the last column "-100.0" which has to be replaced by "0.5"
cs <- read.csv(args[1], sep="\t", stringsAsFactors=FALSE)

pos <- subset(cs, P.Simul.Fst.sample.Fst.>=0.5)
bal <- subset(cs, P.Simul.Fst.sample.Fst.<0.5)
#bal <- cs[(cs$P.Simul.Fst.sample.Fst.<0.5) | (is.na(cs$P.Simul.Fst.sample.Fst.)), ] # Usefull for a different purpose.


Qbal <- p.adjust(bal$P.Simul.Fst.sample.Fst., method="fdr")

invpos <- 1 - pos$P.Simul.Fst.sample.Fst.
Qpos <- 1 - p.adjust(invpos, method="fdr")

pos["Q-values"] <- Qpos

bal["Q-values"] <- Qbal

All <- rbind(bal, pos)

All <- All[with(All, order(as.numeric(row.names(All)))), ]

write.table(All, file=args[2], row.names=FALSE, sep="\t")
