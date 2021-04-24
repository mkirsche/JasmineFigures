library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)
library(reshape2)
library(gridExtra)
library(rlist)
library(here)

find_denovo <- function(variants, outfile) {

  variants <- variants %>% filter(SUPP_VEC == "100000000" | SUPP_VEC == "100000" | SUPP_VEC == "100" | 
                              SUPP_VEC == "100100100" | SUPP_VEC == "100100000" | SUPP_VEC == "100000100" | SUPP_VEC == "100100")
  variants %>% count(SUPP_VEC)
  
  variants$Technology = "All"
  variants <- variants %>% mutate(Technology=ifelse(SUPP_VEC=="100000000","HiFi only", Technology))
  variants <- variants %>% mutate(Technology=ifelse(SUPP_VEC=="100","ONT only", Technology))
  variants <- variants %>% mutate(Technology=ifelse(SUPP_VEC=="100000","CLR only", Technology))
  variants <- variants %>% mutate(Technology=ifelse(SUPP_VEC=="100000100","ONT+HiFi", Technology))
  variants <- variants %>% mutate(Technology=ifelse(SUPP_VEC=="100100000","CLR+HiFi", Technology))
  variants <- variants %>% mutate(Technology=ifelse(SUPP_VEC=="100100","CLR+ONT", Technology))
  plot <- ggplot(variants, aes(x=Technology)) + geom_histogram(stat = "count", fill = "#2093c3") +
    ggtitle("Child-Only Variants") +
    ylab("Count") +
    theme(axis.text.x = element_text(size = 16, angle = 30, margin = margin(t = 21)),
          axis.text.y = element_text(size = 16),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 20),
          legend.text = element_text(size = 16),
          legend.title = element_text(size = 20),
          legend.position = c(legendx, legendy),
          plot.title = element_text(size = 24, hjust = 0.5),
    ) +
  ggsave(outfile, height = 8, width = 8)
}
projectroot <- here()
projectroot <- "/home/mkirsche/jasmine_data"
merged_table_file <- paste(projectroot, "/figures/figure4/denovo.merged.tsv", sep = "")
variants = read.table(merged_table_file, comment.char = "#", sep = "\t", header = T, stringsAsFactors=FALSE)
variants <- variants %>% filter(IS_SPECIFIC == "1" & IS_PRECISE == "1")

outfile <- paste(projectroot, "/figures/figure4/denovocandidates.png", sep = "")
find_denovo(variants, outfile)
outfile <- paste(projectroot, "/figures/figure4/denovocandidates.svg", sep = "")
find_denovo(variants, outfile)


