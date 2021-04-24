library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)
library(reshape2)
library(gridExtra)
library(rlist)

suppvec_hist <- function(df, caller, outfile) {
  df$SUPP_VEC_STRING <- str_pad(as.character(df$SUPP_VEC), 3, "left", "0")
  df$SUPP_VEC_STRING <- ifelse(df$SUPP_VEC_STRING == "100", "Child Only", df$SUPP_VEC_STRING)
  df$SUPP_VEC_STRING <- ifelse(df$SUPP_VEC_STRING == "010", "Father Only", df$SUPP_VEC_STRING)
  df$SUPP_VEC_STRING <- ifelse(df$SUPP_VEC_STRING == "001", "Mother Only", df$SUPP_VEC_STRING)
  df$SUPP_VEC_STRING <- ifelse(df$SUPP_VEC_STRING == "110", "Son/Father", df$SUPP_VEC_STRING)
  df$SUPP_VEC_STRING <- ifelse(df$SUPP_VEC_STRING == "101", "Son/Mother", df$SUPP_VEC_STRING)
  df$SUPP_VEC_STRING <- ifelse(df$SUPP_VEC_STRING == "011", "Both parents", df$SUPP_VEC_STRING)
  df$SUPP_VEC_STRING <- ifelse(df$SUPP_VEC_STRING == "111", "All three", df$SUPP_VEC_STRING)
  
  suppveccounts <- df %>% count(SUPP_VEC_STRING)
  suppveccounts
  suppveccounts$TYPE = "INS"
  ggplot(df, aes(x = SUPP_VEC_STRING, y = 1, fill = TYPE)) +
    geom_bar(position = "stack", stat = "identity") +
    labs(title = paste("SVs by Sample Presence (", caller, ")", sep = "")) +
    xlab("Samples") +
    ylab("Count") +
    scale_y_continuous(expand = expansion(mult = c(0, .1))) +
    theme(plot.title = element_text(size = 18, hjust = 0.5),
          axis.text.x = element_text(size = 12, angle = 30, margin = margin(t = 15)),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 16),
    ) +
    scale_fill_discrete(name = "SVTYPE") +
    geom_text(data = suppveccounts, aes(x = SUPP_VEC_STRING, y=n, label=n), position=position_dodge(width=0.9), vjust=-0.75)
  
  ggsave(outfile, width= 7, height = 8)
}

plot_length <- function(df, caller, outfile) 
{
  df$LenCategory = "TRA"
  
  df$ABSLEN = abs(df$LEN)
  
  df$LenCategory = ifelse(df$ABSLEN >= 0 & df$ABSLEN <= 50, "30-50", df$LenCategory)
  df$LenCategory = ifelse(df$ABSLEN > 50 & df$ABSLEN <= 100, "50-100", df$LenCategory)
  df$LenCategory = ifelse(df$ABSLEN > 100 & df$ABSLEN <= 200, "100-200", df$LenCategory)
  df$LenCategory = ifelse(df$ABSLEN > 200 & df$ABSLEN <= 300, "200-300", df$LenCategory)
  df$LenCategory = ifelse(df$ABSLEN > 300 & df$ABSLEN <= 400, "300-400", df$LenCategory)
  df$LenCategory = ifelse(df$ABSLEN > 400 & df$ABSLEN <= 500, "400-500", df$LenCategory)
  df$LenCategory = ifelse(df$ABSLEN > 500 & df$ABSLEN <= 750, "500-750", df$LenCategory)
  df$LenCategory = ifelse(df$ABSLEN > 750 & df$ABSLEN <= 1000, "750-1k", df$LenCategory)
  df$LenCategory = ifelse(df$ABSLEN > 1000 & df$ABSLEN <= 2000, "1k-2k", df$LenCategory)
  df$LenCategory = ifelse(df$ABSLEN > 2000 & df$ABSLEN <= 5000, "2k-5k", df$LenCategory)
  df$LenCategory = ifelse(df$ABSLEN > 5000 & df$ABSLEN <= 10000, "5k-10k", df$LenCategory)
  df$LenCategory = ifelse(df$ABSLEN > 10000, "10k+", df$LenCategory)
  
  df$LenCategory = ifelse(df$TYPE == "TRA", "TRA", df$LenCategory)
  
  labellist = c("30-50", "50-100", "100-200", "200-300", "300-400", "400-500", "500-750", "750-1k", "1k-2k", "2k-5k", "5k-10k", "10k+", "TRA")
  df$LenCategory <- factor(df$LenCategory,levels = labellist)
  
  df %>% group_by(df$LenCategory) %>% tally()
  
  delcounts <- df %>% filter(df$TYPE == "DEL") %>% group_by(LenCategory) %>% count()
  delcounts$TYPE = "DEL"
  
  nondelcounts <- df %>% filter(df$TYPE != "DEL") %>% group_by(LenCategory) %>% count()
  nondelcounts$TYPE = "INS"
  
  ggplot(df) +
    geom_bar(data = df %>% filter(TYPE != "DEL"), position = "stack", stat = "identity", aes(x = LenCategory, fill = TYPE, y = 1))+
    geom_bar(data = df %>% filter(TYPE == "DEL"), position = "stack", stat = "identity", aes(x = LenCategory, fill = TYPE, y = -1))+
    scale_x_discrete(labels=labellist) +
    xlab("Length") +
    ylab("Count") +
    labs(title = paste("Indels by Length (", caller, ")", sep = "")) +
    theme(plot.title = element_text(size = 24, hjust = 0.5),
          axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 16),
          text = element_text(size = 16),
    ) +
    scale_fill_discrete(name = "Type") +
    geom_text(data = nondelcounts, aes(x = LenCategory, y=n, label=n), position=position_dodge(width=0.9), vjust=-0.75) +
    geom_text(data = delcounts, aes(x = LenCategory, y=-1*n, label=n), position=position_dodge(width=0.9), vjust=1.5)
  
  ggsave(outfile, width= 12, height = 8)
}

plot_allele_frequencies <- function(variants, outfile) {
  variants$SUPP_INT = strtoi(variants$SUPP)
  ggplot(variants, aes(x = SUPP_INT, y = 1, fill = SUPP_INT)) +
    geom_bar(position = "stack", stat = "identity") +
    ggtitle("Population Allele Frequencies") +
    scale_fill_viridis_c(option = "plasma") +
    theme_bw() +
    theme(axis.text.x = element_text(size = 16),
          axis.text.y = element_text(size = 16),
          axis.title.x = element_text(size = 18),
          axis.title.y = element_text(size = 18),
          legend.title = element_blank(),
          legend.text = element_blank(),
          legend.position = "none",
          plot.title = element_text(hjust = 0.5, size = 20)) +
    xlab("Allele Frequency") +
    ylab("Number of variants") +
    
    ggsave(outfile, device = "png", height = 8, width = 8)
}

fn <- '/home/mkirsche/jasmine_data/figures/figure2/hg002_jasmine_noTRA_augmented.tsv'
df <- read.table(fn, sep = "\t", header = TRUE)
df <- df %>% filter(SPECIFIC_FLAG == 1 & PRECISE_FLAG == 1)

nrow(df)
nrow(df %>% filter(EXON_FLAG == 1))
nrow(df %>% filter(GENE_FLAG == 1))
nrow(df %>% filter(REPEAT_FLAG == 1))
nrow(df %>% filter(CENTROMERE_FLAG == 1))

colnames(df)

outfile <- "/home/mkirsche/jasmine_data/figures/figure2/anno_jasmine_all_svhist.png"
suppvec_hist(df, "JasmineAll", outfile)
outfile <- "/home/mkirsche/jasmine_data/figures/figure2/anno_jasmine_all_length.png"
plot_length(df, "JasmineAll", outfile)

outfile <- "/home/mkirsche/jasmine_data/figures/figure2/anno_jasmine_cen_svhist.png"
suppvec_hist(df %>% filter(CENTROMERE_FLAG == 1), "JasmineCentromere", outfile)
outfile <- "/home/mkirsche/jasmine_data/figures/figure2/anno_jasmine_cen_length.png"
plot_length(df %>% filter(CENTROMERE_FLAG == 1), "JasmineCentromere", outfile)

outfile <- "/home/mkirsche/jasmine_data/figures/figure2/anno_jasmine_repeat_svhist.png"
suppvec_hist(df %>% filter(REPEAT_FLAG == 1), "JasmineRepeat", outfile)
outfile <- "/home/mkirsche/jasmine_data/figures/figure2/anno_jasmine_repeat_length.png"
plot_length(df %>% filter(REPEAT_FLAG == 1), "JasmineRepeat", outfile)

outfile <- "/home/mkirsche/jasmine_data/figures/figure2/anno_jasmine_gene_svhist.png"
suppvec_hist(df %>% filter(GENE_FLAG == 1), "JasmineGene", outfile)
outfile <- "/home/mkirsche/jasmine_data/figures/figure2/anno_jasmine_gene_length.png"
plot_length(df %>% filter(GENE_FLAG == 1), "JasmineGene", outfile)

outfile <- "/home/mkirsche/jasmine_data/figures/figure2/anno_jasmine_exon_svhist.png"
suppvec_hist(df %>% filter(EXON_FLAG == 1), "JasmineExon", outfile)
outfile <- "/home/mkirsche/jasmine_data/figures/figure2/anno_jasmine_exon_length.png"
plot_length(df %>% filter(EXON_FLAG == 1), "JasmineExon", outfile)


fn <- '/home/mkirsche/jasmine_data/figures/figure5/all4_pop_jasmine_specprec_noTRA_augmented.tsv'
df <- read.table(fn, sep = "\t", header = TRUE)
df <- df %>% filter(SPECIFIC_FLAG == 1 & PRECISE_FLAG == 1)

nrow(df)
nrow(df %>% filter(EXON_FLAG == 1))
nrow(df %>% filter(GENE_FLAG == 1))
nrow(df %>% filter(REPEAT_FLAG == 1))
nrow(df %>% filter(CENTROMERE_FLAG == 1))

colnames(df)
