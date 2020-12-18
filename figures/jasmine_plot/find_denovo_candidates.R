library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)
library(reshape2)
library(gridExtra)
library(rlist)

find_denovo <- function(merged_table_file, outfile) {
  variants = read.table(merged_table_file, comment.char = "#", sep = "\t", header = T, stringsAsFactors=FALSE)
  
  variants <- variants %>% filter(IS_SPECIFIC == "1" & IS_PRECISE == "1")

  variants <- variants %>% filter(SUPP_VEC == "100000000" | SUPP_VEC == "100000" | SUPP_VEC == "100" | 
                              SUPP_VEC == "100100100" | SUPP_VEC == "100100000" | SUPP_VEC == "100000100" | SUPP_VEC == "100100")
  variants %>% count(SUPP_VEC)
  
  variants$Technology = "All"
  variants <- variants %>% mutate(Technology=ifelse(SUPP_VEC=="100000000","Hifi only", Technology))
  variants <- variants %>% mutate(Technology=ifelse(SUPP_VEC=="100","ONT only", Technology))
  variants <- variants %>% mutate(Technology=ifelse(SUPP_VEC=="100000","CLR only", Technology))
  variants <- variants %>% mutate(Technology=ifelse(SUPP_VEC=="100000100","ONT+Hifi", Technology))
  variants <- variants %>% mutate(Technology=ifelse(SUPP_VEC=="100100000","CLR+Hifi", Technology))
  variants <- variants %>% mutate(Technology=ifelse(SUPP_VEC=="100100","CLR+ONT", Technology))
  plot <- ggplot(variants, aes(x=Technology, fill = Technology)) + geom_histogram(stat = "count") +
    ggtitle("Child-Only Variants") +
    scale_fill_discrete(rainbow) +
    theme_bw() +
    theme(axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          legend.title = element_blank(),
          legend.text = element_blank(),
          legend.position = "none",
          plot.title = element_text(hjust = 0.5, size = 18))
  ggsave(outfile, device = "png", height = 8, width = 8)
}

merged_table_file <- "/home/mkirsche/jasmine_data/figures/figure4/denovo.merged.tsv"
outfile <- "/home/mkirsche/jasmine_data/figures/figure4/denovocandidates.png"
find_denovo(merged_table_file, outfile)