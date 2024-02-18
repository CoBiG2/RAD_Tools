#automatedbayescanplots: plots the FST of each loci and the density of loci with a given selection coefficient.
#based on the file.g_fst and the file.g.sel outpputed by the selection detection program BayeScan.
#Usage:
# BayeScan_plots path to [plot_R.r] [file.g_fst] [g.sel] output_namefst output_namedensity
#Variable to allow using only the files as inputs
plot_bayescan <- function(res, FDR = 0.05, size = 1, pos = 0.35, highlight = NULL, name_highlighted = F, add_text = T)
{
  if (is.character(res))
    res=read.table(res)
  
  colfstat=4
  colq=colfstat-1
  
  highlight_rows=which(is.element(as.numeric(row.names(res)),highlight))
  non_highlight_rows=setdiff(1:nrow(res),highlight_rows)
  
  outliers=as.integer(row.names(res[res[,colq]<=FDR,]))
  
  ok_outliers=TRUE
  if (sum(res[,colq]<=FDR)==0)
    ok_outliers=FALSE;
  
  res[res[,colq]<=0.0001,colq]=0.0001
  
  # plot
  plot(log10(res[,colq]),res[,colfstat],xlim=rev(range(log10(res[,colq]))),xlab="log10(q value)",ylab=names(res[colfstat]),type="n")
  points(log10(res[non_highlight_rows,colq]),res[non_highlight_rows,colfstat],pch=19,cex=size)
  
  if (name_highlighted) {
    if (length(highlight_rows)>0) {
      text(log10(res[highlight_rows,colq]),res[highlight_rows,colfstat],row.names(res[highlight_rows,]),col="red",cex=size*1.2,font=2)
    }
  }
  else {
    points(log10(res[highlight_rows,colq]),res[highlight_rows,colfstat],col="red",pch=19,cex=size)
    # add names of loci over p and vertical line
    if (ok_outliers & add_text) {
      text(log10(res[res[,colq]<=FDR,][,colq])+pos*(round(runif(nrow(res[res[,colq]<=FDR,]),1,2))*2-3),res[res[,colq]<=FDR,][,colfstat],row.names(res[res[,colq]<=FDR,]),cex=size)
    }
  }
  lines(c(log10(FDR),log10(FDR)),c(-1,1),lwd=2)
  
  return(list("outliers"=outliers,"nb_outliers"=length(outliers)))
}
args <- commandArgs(trailingOnly = TRUE)
# Get arguments
#sourcefile <- args[1]
g_fst.txt <- args[1]
g.sel <- args[2]
output_namefst <- args[3]
output_namedensity <- args[4]
#Source plot_R.r
#source(sourcefile::plot_bayescan(res, FDR = 0.05, size = 1, pos = 0.35, highlight = NULL, name_highlighted = F, add_text = T))
#Plot FST
png(paste(output_namefst, ".png", sep=""))
plot_bayescan(g_fst.txt,0,FDR=0.05)
dev.off()
mydata=read.table(g.sel,colClasses="numeric")
parameter="Fst1"
#Plot density of loci with a given selection coefficient.
png(paste(output_namedensity, ".png", sep=""))
plot(density(mydata[[parameter]]), xlab=parameter,
     main=paste(parameter,"posterior distribution"))
dev.off()
#install.packages("boa")
#boa.hpd(mydata[[parameter]],0.05)
