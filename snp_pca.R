#!/usr/bin/Rscript
# snp_pca.R performs a PCA using the SNPRelate R package using a VCF file
# and an option populations files

# Usage:
# snp_pca.R vcf_file output_file_name popupations_file[optional]

library("SNPRelate")
library("plotly")
library(htmlwidgets)
library(htmltools)

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

pc.percent<- pca$varprop * 100
print(round(pc.percent, 2))

# Open figure driver
#pdf(paste(output_name, ".pdf", sep=""))

# Plots PCA
if (!is.na(pops_file)) {
  sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
  pop_code <- read.table(pops_file, sep=",")
  sorted_pops <- pop_code$V2[order(match(pop_code$V1, sample.id))]
  tab <- data.frame(sample.id = pca$sample.id,
    pop = sorted_pops,
    EV1 = pca$eigenvect[,1],
    EV2 = pca$eigenvect[,2],
    stringsAsFactors=F)
  #save(tab, file=paste(output_name, ".Rdata", sep=""))
  p <- plot_ly(tab,  x=tab$EV1, y=tab$EV2, text=tab$sample.id, color=tab$pop, colors="Set3")
  p <- layout(p, title="PCA",
              xaxis=list(title=paste("PC 1(", round(pca$eigenval[1], d=2) , "%)")),
              yaxis=list(title=paste("PC 1(", round(pca$eigenval[2], d=2) , "%)")))
  htmlwidgets::saveWidget(as.widget(p), paste(output_name, ".html", sep=""))
} else {
  tab <- data.frame(sample.id = pca$sample.id,
    EV1 = pca$eigenvect[, 1],
    EV2 = pca$eigenvect[, 2],
    stringsAsFactors=F)
  print(pca$sample.id)
  print(tab$EV1)
  print(tab$EV2)
  p <- plot_ly(tab, x=tab$EV1, y=tab$EV2, text=tab$sample.id)
  p <- layout(p, title="PCA",
              xaxis=list(title=paste("PC 1(", round(pca$eigenval[1], d=2) , "%)"),
              yaxis=list(title=paste("PC 2(", round(pca$eigenval[2], d=2) , "%)"))))
}

p <- htmlwidgets::appendContent(p,
                                htmltools::tags$input(id='inputText',
                                                      value='', ''),
                                htmltools::tags$button(id='buttonSearch',
                                                       'Search'))

p <- htmlwidgets::appendContent(p, htmltools::tags$script(HTML(
'document.getElementById("buttonSearch").addEventListener("click", function()
  {
    var i = 0;
    var j = 0;
    var found = [];
    var myDiv = document.getElementsByClassName("js-plotly-plot")[0]
    var data = JSON.parse(document.querySelectorAll("script[type=\'application/json\']")[0].innerHTML);
    console.log(data.x.data)
    for (i = 0 ;i < data.x.data.length; i += 1) {
      for (j = 0; j < data.x.data[i].text.length; j += 1) {
        if (data.x.data[i].text[j].indexOf(document.getElementById("inputText").value) !== -1) {
          found.push({curveNumber: i, pointNumber: j});
        }
      }
    }
    Plotly.Fx.hover(myDiv, found);
  }
);')))
htmlwidgets::saveWidget(as.widget(p), paste(output_name, ".html", sep=""))

# remove temporary gds file
file.remove("temp.gds")
