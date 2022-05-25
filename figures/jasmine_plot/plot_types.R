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
  summarized <- df %>% group_by(SVTYPE, LenCategory) %>% summarise(counts=n())
  plot <- ggplot(summarized) +
    geom_bar(data = summarized %>% filter(SVTYPE != "DEL"), position = "stack", stat = "identity", aes(x = LenCategory, fill = SVTYPE, y = counts))+
    geom_bar(data = summarized %>% filter(SVTYPE == "DEL"), position = "stack", stat = "identity", aes(x = LenCategory, fill = SVTYPE, y = -counts))+
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
    guides(fill=guide_legend(title="Type")) +
    geom_text(data = nondelcounts, aes(x = LenCategory, y=n, label=n), position=position_dodge(width=0.9), vjust=-0.75) +
    geom_text(data = delcounts, aes(x = LenCategory, y=-1*n, label=n), position=position_dodge(width=0.9), vjust=1.5)
  return(plot)
}

suppvec_hist_highlightdiscordant <- function(df, caller,div, lengthfilter) {
  offset <- nrow(df)/div
  df$SUPP_VEC_STRING <- str_pad(as.character(df$SUPP_VEC), 3, "left", "0")
  df$SUPP_VEC_STRING <- ifelse(df$SUPP_VEC_STRING == "100", "HG002 Only", df$SUPP_VEC_STRING)
  df$SUPP_VEC_STRING <- ifelse(df$SUPP_VEC_STRING == "010", "HG003 Only", df$SUPP_VEC_STRING)
  df$SUPP_VEC_STRING <- ifelse(df$SUPP_VEC_STRING == "001", "HG004 Only", df$SUPP_VEC_STRING)
  df$SUPP_VEC_STRING <- ifelse(df$SUPP_VEC_STRING == "110", "HG002/HG003", df$SUPP_VEC_STRING)
  df$SUPP_VEC_STRING <- ifelse(df$SUPP_VEC_STRING == "101", "HG002/HG004", df$SUPP_VEC_STRING)
  df$SUPP_VEC_STRING <- ifelse(df$SUPP_VEC_STRING == "011", "HG003/HG004", df$SUPP_VEC_STRING)
  df$SUPP_VEC_STRING <- ifelse(df$SUPP_VEC_STRING == "111", "All three", df$SUPP_VEC_STRING)
  df$SUPP_VEC_STRING <- factor(df$SUPP_VEC_STRING,levels = c("All three", "HG003/HG004", "HG003 Only", 
                                                             "HG004 Only", "HG002/HG003", "HG002/HG004", "HG002 Only"))
  colorpalette <- c(DEL = '#CC6677',DUP = '#332288',INS = '#44AA99',TRA = '#88CCEE',INV = '#DDCC77')
  titlestart <- "Variants by Sample Presence"
  if(lengthfilter)
  {
    titlestart <- "SVs by Sample Presence"
    df <- df %>% filter(SVLEN == 0 | abs(SVLEN) >= 50 | SVTYPE == "TRA")
    
  }
  suppveccounts <- df %>% count(SUPP_VEC_STRING)
  suppveccounts
  suppveccounts$SVTYPE = "INS"
  
  childcounts <- df %>% filter(df$SUPP_VEC_STRING == "HG002 Only") %>% group_by(SUPP_VEC_STRING) %>% count()
  childcounts$SVTYPE = "INS"
  childcounts$shape = 23
  
  summarized <- df %>% group_by(SVTYPE, SUPP_VEC_STRING) %>% summarise(counts=n())
  
  
  plot <- ggplot(summarized, aes(x = SUPP_VEC_STRING, y = counts, fill = SVTYPE)) +
    geom_bar(position = "stack", stat = "identity") +
    labs(title = paste(titlestart, " (", caller, ")", sep = "")) +
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
    guides(fill=guide_legend(title="Type")) +
    geom_text(data = suppveccounts, aes(x = SUPP_VEC_STRING, y=n, label=n), position=position_dodge(width=0.9), vjust=-0.75) +
    geom_point(data = childcounts, aes(x = SUPP_VEC_STRING, y=n+offset), shape = 25, fill = "darkred", color = "darkred", size = 5)
  return(plot)
}

suppvec_hist_crosstech <- function(df, caller, lengthfilter) {
  df$SUPP_VEC_STRING <- str_pad(as.character(df$SUPP_VEC), 3, "left", "0")
  df$SUPP_VEC_STRING <- ifelse(df$SUPP_VEC_STRING == "100", "Hifi Only", df$SUPP_VEC_STRING)
  df$SUPP_VEC_STRING <- ifelse(df$SUPP_VEC_STRING == "010", "CLR Only", df$SUPP_VEC_STRING)
  df$SUPP_VEC_STRING <- ifelse(df$SUPP_VEC_STRING == "001", "ONT Only", df$SUPP_VEC_STRING)
  df$SUPP_VEC_STRING <- ifelse(df$SUPP_VEC_STRING == "110", "Hifi/CLR", df$SUPP_VEC_STRING)
  df$SUPP_VEC_STRING <- ifelse(df$SUPP_VEC_STRING == "101", "Hifi/ONT", df$SUPP_VEC_STRING)
  df$SUPP_VEC_STRING <- ifelse(df$SUPP_VEC_STRING == "011", "CLR/ONT", df$SUPP_VEC_STRING)
  df$SUPP_VEC_STRING <- ifelse(df$SUPP_VEC_STRING == "111", "All three", df$SUPP_VEC_STRING)
  colorpalette <- c(DEL = '#CC6677',DUP = '#332288',INS = '#44AA99',TRA = '#88CCEE',INV = '#DDCC77')
  titlestart <- "Variants by Technology"
  if(lengthfilter)
  {
    titlestart <- "SVs by Technology"
    df <- df %>% filter(SVLEN == 0 | abs(SVLEN) >= 50 | SVTYPE == "TRA")
  }
  suppveccounts <- df %>% count(SUPP_VEC_STRING)
  suppveccounts
  suppveccounts$SVTYPE = "INS"
  summarized <- df %>% group_by(SVTYPE, SUPP_VEC_STRING) %>% summarise(counts=n())
  
  plot <- ggplot(summarized, aes(x = SUPP_VEC_STRING, y = counts, fill = SVTYPE)) +
    geom_bar(position = "stack", stat = "identity") +
    labs(title = paste(titlestart, " (", caller, ")", sep = "")) +
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
    guides(fill=guide_legend(title="Type")) +
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
lowcomplexity <- plot_table(filtered %>% filter(INTERSECTS_LOW_COMPLEXITY == 1), "Low Complexity")
ggarrange(all,exon, gene, line, sine, satellite, simplerepeat, lowcomplexity, centromere,nrow = 3, ncol=3)
ggsave("/home/mkirsche/jasmine_data/beds/hg002_hifi_jasmine_anno.png", width=20, height=20)
ggarrange(all,exon, gene, line, sine, satellite, simplerepeat, lowcomplexity, centromere,nrow = 3, ncol=3)
ggsave("/home/mkirsche/jasmine_data/beds/hg002_hifi_jasmine_anno.svg", width=20, height=20)
div <-22
all <- suppvec_hist_highlightdiscordant(filtered, "All Variants",div, FALSE)
exon <- suppvec_hist_highlightdiscordant(filtered %>% filter(INTERSECTS_EXON == 1), "Exon",div, FALSE)
gene <- suppvec_hist_highlightdiscordant(filtered %>% filter(INTERSECTS_GENE == 1), "Gene",div, FALSE)
centromere <- suppvec_hist_highlightdiscordant(filtered %>% filter(INTERSECTS_CENTROMERES == 1), "Centromere",div,FALSE)
line <- suppvec_hist_highlightdiscordant(filtered %>% filter(INTERSECTS_LINE == 1), "LINE",div, FALSE)
satellite <- suppvec_hist_highlightdiscordant(filtered %>% filter(INTERSECTS_SATELLITE == 1), "Satellite",div,FALSE)
simplerepeat <- suppvec_hist_highlightdiscordant(filtered %>% filter(INTERSECTS_SIMPLE_REPEAT == 1), "Simple Repeat",div,FALSE)
sine <- suppvec_hist_highlightdiscordant(filtered %>% filter(INTERSECTS_SINE == 1), "SINE",div, FALSE)
lowcomplexity <- suppvec_hist_highlightdiscordant(filtered %>% filter(INTERSECTS_LOW_COMPLEXITY == 1), "Low Complexity",div, FALSE)
ggarrange(all,exon, gene, line, sine, satellite, simplerepeat, lowcomplexity, centromere,nrow = 3, ncol=3)
ggsave("/home/mkirsche/jasmine_data/beds/hg002_hifi_jasmine_anno_suppvec.png", width=20, height=20)
ggarrange(all,exon, gene, line, sine, satellite, simplerepeat, lowcomplexity, centromere,nrow = 3, ncol=3)
ggsave("/home/mkirsche/jasmine_data/beds/hg002_hifi_jasmine_anno_suppvec.svg", width=20, height=20)

all <- suppvec_hist_highlightdiscordant(filtered, "All Variants",div, TRUE)
exon <- suppvec_hist_highlightdiscordant(filtered %>% filter(INTERSECTS_EXON == 1), "Exon",div, TRUE)
gene <- suppvec_hist_highlightdiscordant(filtered %>% filter(INTERSECTS_GENE == 1), "Gene",div, TRUE)
centromere <- suppvec_hist_highlightdiscordant(filtered %>% filter(INTERSECTS_CENTROMERES == 1), "Centromere",div,TRUE)
line <- suppvec_hist_highlightdiscordant(filtered %>% filter(INTERSECTS_LINE == 1), "LINE",div, TRUE)
satellite <- suppvec_hist_highlightdiscordant(filtered %>% filter(INTERSECTS_SATELLITE == 1), "Satellite",div,TRUE)
simplerepeat <- suppvec_hist_highlightdiscordant(filtered %>% filter(INTERSECTS_SIMPLE_REPEAT == 1), "Simple Repeat",div,TRUE)
sine <- suppvec_hist_highlightdiscordant(filtered %>% filter(INTERSECTS_SINE == 1), "SINE",div, TRUE)
lowcomplexity <- suppvec_hist_highlightdiscordant(filtered %>% filter(INTERSECTS_LOW_COMPLEXITY == 1), "Low Complexity",div, TRUE)
ggarrange(all,exon, gene, line, sine, satellite, simplerepeat, lowcomplexity, centromere,nrow = 3, ncol=3)
ggsave("/home/mkirsche/jasmine_data/beds/hg002_hifi_jasmine_anno_suppvec_50plus.png", width=20, height=20)
ggarrange(all,exon, gene, line, sine, satellite, simplerepeat, lowcomplexity, centromere,nrow = 3, ncol=3)
ggsave("/home/mkirsche/jasmine_data/beds/hg002_hifi_jasmine_anno_suppvec_50plus.svg", width=20, height=20)

fn <- "/home/mkirsche/jasmine_data/beds/crosstech_annotated.tsv"
df <- read.table(fn, sep = "\t", header = TRUE)
colnames(df)
filtered <- df %>% filter(IS_SPECIFIC == 1 & IS_PRECISE == 1 & SVTYPE != "TRA")
all <- suppvec_hist_crosstech(filtered, "All Variants", FALSE)
exon <- suppvec_hist_crosstech(filtered %>% filter(INTERSECTS_EXON == 1), "Exon", FALSE)
gene <- suppvec_hist_crosstech(filtered %>% filter(INTERSECTS_GENE == 1), "Gene", FALSE)
centromere <- suppvec_hist_crosstech(filtered %>% filter(INTERSECTS_CENTROMERES == 1), "Centromere", FALSE)
line <- suppvec_hist_crosstech(filtered %>% filter(INTERSECTS_LINE == 1), "LINE", FALSE)
satellite <- suppvec_hist_crosstech(filtered %>% filter(INTERSECTS_SATELLITE == 1), "Satellite",FALSE)
simplerepeat <- suppvec_hist_crosstech(filtered %>% filter(INTERSECTS_SIMPLE_REPEAT == 1), "Simple Repeat", FALSE)
sine <- suppvec_hist_crosstech(filtered %>% filter(INTERSECTS_SINE == 1), "SINE", FALSE)
lowcomplexity <- suppvec_hist_crosstech(filtered %>% filter(INTERSECTS_LOW_COMPLEXITY == 1), "Low Complexity", FALSE)
ggarrange(all,exon, gene, line, sine, satellite, simplerepeat, lowcomplexity, centromere,nrow = 3, ncol=3)
ggsave("/home/mkirsche/jasmine_data/beds/crosstech_jasmine_anno_suppvec.png", width=20, height=20)
ggarrange(all,exon, gene, line, sine, satellite, simplerepeat, lowcomplexity, centromere,nrow = 3, ncol=3)
ggsave("/home/mkirsche/jasmine_data/beds/crosstech_jasmine_anno_suppvec.svg", width=20, height=20)

all <- suppvec_hist_crosstech(filtered, "All Variants", TRUE)
exon <- suppvec_hist_crosstech(filtered %>% filter(INTERSECTS_EXON == 1), "Exon", TRUE)
gene <- suppvec_hist_crosstech(filtered %>% filter(INTERSECTS_GENE == 1), "Gene", TRUE)
centromere <- suppvec_hist_crosstech(filtered %>% filter(INTERSECTS_CENTROMERES == 1), "Centromere", TRUE)
line <- suppvec_hist_crosstech(filtered %>% filter(INTERSECTS_LINE == 1), "LINE", TRUE)
satellite <- suppvec_hist_crosstech(filtered %>% filter(INTERSECTS_SATELLITE == 1), "Satellite",TRUE)
simplerepeat <- suppvec_hist_crosstech(filtered %>% filter(INTERSECTS_SIMPLE_REPEAT == 1), "Simple Repeat", TRUE)
sine <- suppvec_hist_crosstech(filtered %>% filter(INTERSECTS_SINE == 1), "SINE", TRUE)
lowcomplexity <- suppvec_hist_crosstech(filtered %>% filter(INTERSECTS_LOW_COMPLEXITY == 1), "Low Complexity", TRUE)
ggarrange(all,exon, gene, line, sine, satellite, simplerepeat, lowcomplexity, centromere,nrow = 3, ncol=3)
ggsave("/home/mkirsche/jasmine_data/beds/crosstech_jasmine_anno_suppvec_50plus.png", width=20, height=20)
ggarrange(all,exon, gene, line, sine, satellite, simplerepeat, lowcomplexity, centromere,nrow = 3, ncol=3)
ggsave("/home/mkirsche/jasmine_data/beds/crosstech_jasmine_anno_suppvec_50plus.svg", width=20, height=20)


#fn <- "/home/mkirsche/jasmine_data/beds/highaf_annotated.tsv"
#df <- read.table(fn, sep = "\t", header = TRUE)
#colnames(df)
#filtered <- df %>% filter(IS_SPECIFIC == 1 & IS_PRECISE == 1 & SVTYPE != "TRA")
#all <- plot_table(filtered, "All Variants")
#exon <- plot_table(filtered %>% filter(INTERSECTS_EXON == 1), "Exon")
#gene <- plot_table(filtered %>% filter(INTERSECTS_GENE == 1), "Gene")
#centromere <- plot_table(filtered %>% filter(INTERSECTS_CENTROMERES == 1), "Centromere")
#line <- plot_table(filtered %>% filter(INTERSECTS_LINE == 1), "LINE")
#satellite <- plot_table(filtered %>% filter(INTERSECTS_SATELLITE == 1), "Satellite")
#simplerepeat <- plot_table(filtered %>% filter(INTERSECTS_SIMPLE_REPEAT == 1), "Simple Repeat")
#sine <- plot_table(filtered %>% filter(INTERSECTS_SINE == 1), "SINE")
#lowcomplexity <- plot_table(filtered %>% filter(INTERSECTS_LOW_COMPLEXITY == 1), "Low Complexity Repeat")
#ggarrange(all,exon, gene, line, sine, satellite, simplerepeat, lowcomplexity, centromere,nrow = 3, ncol=3)
#nrow(filtered)
#nrow(filtered %>% filter(INTERSECTS_CENTROMERES == 1))
#ggsave("/home/mkirsche/jasmine_data/beds/highaf_anno.png", width=20, height=20)
