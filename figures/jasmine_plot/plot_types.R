library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)
library(reshape2)
library(gridExtra)
library(rlist)
library(VennDiagram)
library(Cairo)
library(ggpubr)
library(here)

plot_table <- function(df, title)
{
  df$LenCategory = "TRA"
  
  df$ABSLEN = abs(df$SVLEN)
  
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
  
  df$LenCategory = ifelse(df$SVTYPE == "TRA", "TRA", df$LenCategory)
  
  labellist = c("30-50", "50-100", "100-200", "200-300", "300-400", "400-500", "500-750", "750-1k", "1k-2k", "2k-5k", "5k-10k", "10k+", "TRA")
  df$LenCategory <- factor(df$LenCategory,levels = labellist)
  
  colorpalette <- c(DEL = '#CC6677',DUP = '#332288',INS = '#44AA99',TRA = '#88CCEE',INV = '#DDCC77')
  
  
  df %>% group_by(df$LenCategory) %>% tally()
  
  delcounts <- df %>% filter(df$SVTYPE == "DEL") %>% group_by(LenCategory) %>% count()
  delcounts$SVTYPE = "DEL"
  
  nondelcounts <- df %>% filter(df$SVTYPE != "DEL") %>% group_by(LenCategory) %>% count()
  nondelcounts$SVTYPE = "INS"
  
  plot <- ggplot(df) +
    geom_bar(data = df %>% filter(SVTYPE != "DEL"), position = "stack", stat = "identity", aes(x = LenCategory, fill = SVTYPE, y = 1))+
    geom_bar(data = df %>% filter(SVTYPE == "DEL"), position = "stack", stat = "identity", aes(x = LenCategory, fill = SVTYPE, y = -1))+
    scale_x_discrete(labels=labellist) +
    xlab("Length") +
    ylab("Count") +
    labs(title = paste("Indels by Length (", title, ")", sep = "")) +
    theme(plot.title = element_text(size = 24, hjust = 0.5),
          axis.text.x = element_text(size = 14, angle=30),
          axis.text.y = element_text(size = 14),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 16),
          text = element_text(size = 16),
          legend.position = "None",
    ) +
    scale_fill_manual(name = "SVTYPE", values = colorpalette) +
    geom_text(data = nondelcounts, aes(x = LenCategory, y=n, label=n), position=position_dodge(width=0.9), vjust=-0.75) +
    geom_text(data = delcounts, aes(x = LenCategory, y=-1*n, label=n), position=position_dodge(width=0.9), vjust=1.5)
  return(plot)
}

suppvec_hist_highlightdiscordant <- function(df, caller,div) {
  offset <- nrow(df)/div
  df$SUPP_VEC_STRING <- str_pad(as.character(df$SUPP_VEC), 3, "left", "0")
  df$SUPP_VEC_STRING <- ifelse(df$SUPP_VEC_STRING == "100", "Child Only", df$SUPP_VEC_STRING)
  df$SUPP_VEC_STRING <- ifelse(df$SUPP_VEC_STRING == "010", "Father Only", df$SUPP_VEC_STRING)
  df$SUPP_VEC_STRING <- ifelse(df$SUPP_VEC_STRING == "001", "Mother Only", df$SUPP_VEC_STRING)
  df$SUPP_VEC_STRING <- ifelse(df$SUPP_VEC_STRING == "110", "Son/Father", df$SUPP_VEC_STRING)
  df$SUPP_VEC_STRING <- ifelse(df$SUPP_VEC_STRING == "101", "Son/Mother", df$SUPP_VEC_STRING)
  df$SUPP_VEC_STRING <- ifelse(df$SUPP_VEC_STRING == "011", "Both parents", df$SUPP_VEC_STRING)
  df$SUPP_VEC_STRING <- ifelse(df$SUPP_VEC_STRING == "111", "All three", df$SUPP_VEC_STRING)
  df$SUPP_VEC_STRING <- factor(df$SUPP_VEC_STRING,levels = c("All three", "Both parents", "Father Only", 
                                                             "Mother Only", "Son/Father", "Son/Mother", "Child Only"))
  colorpalette <- c(DEL = '#CC6677',DUP = '#332288',INS = '#44AA99',TRA = '#88CCEE',INV = '#DDCC77')
  
  suppveccounts <- df %>% count(SUPP_VEC_STRING)
  suppveccounts
  suppveccounts$SVTYPE = "INS"
  
  childcounts <- df %>% filter(df$SUPP_VEC_STRING == "Child Only") %>% group_by(SUPP_VEC_STRING) %>% count()
  childcounts$SVTYPE = "INS"
  childcounts$shape = 23
  plot <- ggplot(df, aes(x = SUPP_VEC_STRING, y = 1, fill = SVTYPE)) +
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
          legend.position = "None"#c(legendx, legendy),
    ) +
    scale_fill_manual(name = "SVTYPE", values = colorpalette) +
    geom_text(data = suppveccounts, aes(x = SUPP_VEC_STRING, y=n, label=n), position=position_dodge(width=0.9), vjust=-0.75) +
    geom_point(data = childcounts, aes(x = SUPP_VEC_STRING, y=n+offset), shape = 25, fill = "darkred", color = "darkred", size = 5)
  return(plot)
}

suppvec_hist_crosstech <- function(df, caller) {
  df$SUPP_VEC_STRING <- str_pad(as.character(df$SUPP_VEC), 3, "left", "0")
  df$SUPP_VEC_STRING <- ifelse(df$SUPP_VEC_STRING == "100", "Hifi Only", df$SUPP_VEC_STRING)
  df$SUPP_VEC_STRING <- ifelse(df$SUPP_VEC_STRING == "010", "CLR Only", df$SUPP_VEC_STRING)
  df$SUPP_VEC_STRING <- ifelse(df$SUPP_VEC_STRING == "001", "ONT Only", df$SUPP_VEC_STRING)
  df$SUPP_VEC_STRING <- ifelse(df$SUPP_VEC_STRING == "110", "Hifi/CLR", df$SUPP_VEC_STRING)
  df$SUPP_VEC_STRING <- ifelse(df$SUPP_VEC_STRING == "101", "Hifi/ONT", df$SUPP_VEC_STRING)
  df$SUPP_VEC_STRING <- ifelse(df$SUPP_VEC_STRING == "011", "CLR/ONT", df$SUPP_VEC_STRING)
  df$SUPP_VEC_STRING <- ifelse(df$SUPP_VEC_STRING == "111", "All three", df$SUPP_VEC_STRING)
  colorpalette <- c(DEL = '#CC6677',DUP = '#332288',INS = '#44AA99',TRA = '#88CCEE',INV = '#DDCC77')
  
  suppveccounts <- df %>% count(SUPP_VEC_STRING)
  suppveccounts
  suppveccounts$SVTYPE = "INS"
  
  plot <- ggplot(df, aes(x = SUPP_VEC_STRING, y = 1, fill = SVTYPE)) +
    geom_bar(position = "stack", stat = "identity") +
    labs(title = paste("SVs by Sequencing Technology (", caller, ")", sep = "")) +
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
          legend.position = "None"#c(legendx, legendy),
    ) +
    scale_fill_manual(name = "SVTYPE", values = colorpalette) +
    geom_text(data = suppveccounts, aes(x = SUPP_VEC_STRING, y=n, label=n), position=position_dodge(width=0.9), vjust=-0.75)
  return(plot)
}

fn <- "/home/mkirsche/jasmine_data/beds/hg002_hifi_jasmine_annotated.tsv"
df <- read.table(fn, sep = "\t", header = TRUE)
colnames(df)
filtered <- df %>% filter(IS_SPECIFIC == 1 & IS_PRECISE == 1 & SVTYPE != "TRA")
all <- plot_table(filtered, "All Variants")
exon <- plot_table(filtered %>% filter(INTERSECTS_EXON == 1), "Exon")
gene <- plot_table(filtered %>% filter(INTERSECTS_GENE == 1), "Gene")
centromere <- plot_table(filtered %>% filter(INTERSECTS_CENTROMERES == 1), "Centromere")
line <- plot_table(filtered %>% filter(INTERSECTS_LINE == 1), "LINE")
satellite <- plot_table(filtered %>% filter(INTERSECTS_SATELLITE == 1), "Satellite")
simplerepeat <- plot_table(filtered %>% filter(INTERSECTS_SIMPLE_REPEAT == 1), "Simple Repeat")
sine <- plot_table(filtered %>% filter(INTERSECTS_SINE == 1), "SINE")
lowcomplexity <- plot_table(filtered %>% filter(INTERSECTS_LOW_COMPLEXITY == 1), "Low Complexity Repeat")
ggarrange(all,exon, gene, line, sine, satellite, simplerepeat, lowcomplexity, centromere,nrow = 3, ncol=3)
ggsave("/home/mkirsche/jasmine_data/beds/hg002_hifi_jasmine_anno.png", width=20, height=20)
div <-22
all <- suppvec_hist_highlightdiscordant(filtered, "All Variants",div)
exon <- suppvec_hist_highlightdiscordant(filtered %>% filter(INTERSECTS_EXON == 1), "Exon",div)
gene <- suppvec_hist_highlightdiscordant(filtered %>% filter(INTERSECTS_GENE == 1), "Gene",div)
centromere <- suppvec_hist_highlightdiscordant(filtered %>% filter(INTERSECTS_CENTROMERES == 1), "Centromere",div)
line <- suppvec_hist_highlightdiscordant(filtered %>% filter(INTERSECTS_LINE == 1), "LINE",div)
satellite <- suppvec_hist_highlightdiscordant(filtered %>% filter(INTERSECTS_SATELLITE == 1), "Satellite",div)
simplerepeat <- suppvec_hist_highlightdiscordant(filtered %>% filter(INTERSECTS_SIMPLE_REPEAT == 1), "Simple Repeat",div)
sine <- suppvec_hist_highlightdiscordant(filtered %>% filter(INTERSECTS_SINE == 1), "SINE",div)
lowcomplexity <- suppvec_hist_highlightdiscordant(filtered %>% filter(INTERSECTS_LOW_COMPLEXITY == 1), "Low Complexity Repeat",div)
ggarrange(all,exon, gene, line, sine, satellite, simplerepeat, lowcomplexity, centromere,nrow = 3, ncol=3)
ggsave("/home/mkirsche/jasmine_data/beds/hg002_hifi_jasmine_anno_suppvec.png", width=20, height=20)

fn <- "/home/mkirsche/jasmine_data/beds/crosstech_annotated.tsv"
df <- read.table(fn, sep = "\t", header = TRUE)
colnames(df)
filtered <- df %>% filter(IS_SPECIFIC == 1 & IS_PRECISE == 1 & SVTYPE != "TRA")
all <- suppvec_hist_crosstech(filtered, "All Variants")
exon <- suppvec_hist_crosstech(filtered %>% filter(INTERSECTS_EXON == 1), "Exon")
gene <- suppvec_hist_crosstech(filtered %>% filter(INTERSECTS_GENE == 1), "Gene")
centromere <- suppvec_hist_crosstech(filtered %>% filter(INTERSECTS_CENTROMERES == 1), "Centromere")
line <- suppvec_hist_crosstech(filtered %>% filter(INTERSECTS_LINE == 1), "LINE")
satellite <- suppvec_hist_crosstech(filtered %>% filter(INTERSECTS_SATELLITE == 1), "Satellite")
simplerepeat <- suppvec_hist_crosstech(filtered %>% filter(INTERSECTS_SIMPLE_REPEAT == 1), "Simple Repeat")
sine <- suppvec_hist_crosstech(filtered %>% filter(INTERSECTS_SINE == 1), "SINE")
lowcomplexity <- suppvec_hist_crosstech(filtered %>% filter(INTERSECTS_LOW_COMPLEXITY == 1), "Low Complexity Repeat")
ggarrange(all,exon, gene, line, sine, satellite, simplerepeat, lowcomplexity, centromere,nrow = 3, ncol=3)
ggsave("/home/mkirsche/jasmine_data/beds/crosstech_jasmine_anno_suppvec.png", width=20, height=20)


fn <- "/home/mkirsche/jasmine_data/beds/highaf_annotated.tsv"
df <- read.table(fn, sep = "\t", header = TRUE)
colnames(df)
filtered <- df %>% filter(IS_SPECIFIC == 1 & IS_PRECISE == 1 & SVTYPE != "TRA")
all <- plot_table(filtered, "All Variants")
exon <- plot_table(filtered %>% filter(INTERSECTS_EXON == 1), "Exon")
gene <- plot_table(filtered %>% filter(INTERSECTS_GENE == 1), "Gene")
centromere <- plot_table(filtered %>% filter(INTERSECTS_CENTROMERES == 1), "Centromere")
line <- plot_table(filtered %>% filter(INTERSECTS_LINE == 1), "LINE")
satellite <- plot_table(filtered %>% filter(INTERSECTS_SATELLITE == 1), "Satellite")
simplerepeat <- plot_table(filtered %>% filter(INTERSECTS_SIMPLE_REPEAT == 1), "Simple Repeat")
sine <- plot_table(filtered %>% filter(INTERSECTS_SINE == 1), "SINE")
lowcomplexity <- plot_table(filtered %>% filter(INTERSECTS_LOW_COMPLEXITY == 1), "Low Complexity Repeat")
ggarrange(all,exon, gene, line, sine, satellite, simplerepeat, lowcomplexity, centromere,nrow = 3, ncol=3)
nrow(filtered)
nrow(filtered %>% filter(INTERSECTS_CENTROMERES == 1))
ggsave("/home/mkirsche/jasmine_data/beds/highaf_anno.png", width=20, height=20)
