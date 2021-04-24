library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)
library(reshape2)
library(gridExtra)
library(rlist)
#library(ggfortify)
library(tidyverse)

#fn<-"/home/mkirsche/1kgp_100samples.matrix.txt"
#fn <- "/home/mkirsche/eclipse-workspace/ParagraphGenotypes/testmatrix.txt"
fn <- "/home/mkirsche/jasmine_data/figures/figure5/paragraph.matrix.txt"
outdir <- "/home/mkirsche/jasmine_data/figures/figure5/"
prefix <- "paragraph_all"
plot_paragraph(fn, outdir, prefix)

fn <- "/home/mkirsche/eclipse-workspace/ParagraphGenotypes/testmatrix_ins.txt"
outdir <- "/home/mkirsche/eclipse-workspace/ParagraphGenotypes/"
prefix <- "paragraph_ins"
plot_paragraph(fn, outdir, prefix)

fn <- "/home/mkirsche/eclipse-workspace/ParagraphGenotypes/testmatrix_del.txt"
outdir <- "/home/mkirsche/eclipse-workspace/ParagraphGenotypes/"
prefix <- "paragraph_del"
plot_paragraph(fn, outdir, prefix)

plot_paragraph <- function(fn, outdir, prefix) {
  df <- read.table(fn, sep = "\t", header = T)#, nrow=2000000)
  
  popfn <- "/home/mkirsche/eclipse-workspace/ParagraphGenotypes/all.idx"
  popdf <- read.table(popfn, sep = " ", header = F)
  colnames(popdf)
  
  names = c("V2")
  l <- c(colnames(df))
  colordf <- data.frame(matrix(unlist(l), nrow=length(l), byrow=TRUE))
  colnames(colordf) <- names
  colordf
  merged <- merge.data.frame(colordf, popdf, all.x=TRUE)
  populationlist = merged$V3
  populationlist
  
  colnames(df)
  rsums <- rowSums(df)
  rsums
  ggplot() + aes(rsums)+ geom_histogram(binwidth=1, colour="black", fill="white")
  ofn <- paste(outdir, prefix, "_allelefreq.png", sep = "")
  #ofn <- "/home/mkirsche/1kgp_af.png"
  ggsave(ofn, width= 8, height = 8)
  ggplot() + aes(rsums[unlist(rsums > 0)])+ geom_histogram(binwidth=1, colour="black", fill="white") #+ xlim(1, 445)
  #ofn <- "/home/mkirsche/1kgp_af_nonzero.png"
  ofn <- paste(outdir, prefix, "_allelefreq_nonzero.png", sep = "")
  ggsave(ofn, width= 8, height = 8)
  #mm1 <- as.matrix(df)
  #mm2 <- matrix(mm1, ncol = ncol(df), dimnames = NULL)
  #pca <- prcomp(t(mm2), scale = FALSE)
  #pca$x
  #var_explained <- pca$sdev^2/sum(pca$sdev^2)
  #pca$x %>% as.data.frame %>% ggplot(aes(x=PC1,y=PC2, color = populationlist)) + geom_point(size=1) +
  #  theme_bw(base_size=32) + 
  #  labs(x=paste0("PC1: ",round(var_explained[1]*100,3),"%"),
  #       y=paste0("PC2: ",round(var_explained[2]*100,3),"%"))
  #ofn <- "/home/mkirsche/1kgp_100samples_pca.png"
  #ofn <- paste(outdir, prefix, "_pca.png", sep = "")
  #ggsave(ofn, width= 8, height = 8)
}

fn <- "/home/mkirsche/eclipse-workspace/ParagraphGenotypes/testmatrix_ins.txt"
outdir <- "/home/mkirsche/eclipse-workspace/ParagraphGenotypes/"
prefix <- "paragraph_ins"
plot_paragraph(fn, outdir, prefix)

fn <- "/home/mkirsche/eclipse-workspace/ParagraphGenotypes/testmatrix_del.txt"
prefix <- "paragraph_del"
plot_paragraph(fn, outdir, prefix)

fn <- "/home/mkirsche/eclipse-workspace/ParagraphGenotypes/testmatrix_del_auto.txt"
prefix <- "paragraph_del_auto"
plot_paragraph(fn, outdir, prefix)

fn <- "/home/mkirsche/eclipse-workspace/ParagraphGenotypes/testmatrix_ins_auto.txt"
prefix <- "paragraph_ins_auto"
plot_paragraph(fn, outdir, prefix)

