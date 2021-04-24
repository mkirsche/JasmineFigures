library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)
library(reshape2)
library(gridExtra)
library(rlist)
library(VennDiagram)
library(Cairo)
library(here)

plot_length <- function(df, caller, outfile) 
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
  
  df %>% group_by(df$LenCategory) %>% tally()
  
  delcounts <- df %>% filter(df$SVTYPE == "DEL") %>% group_by(LenCategory) %>% count()
  delcounts$SVTYPE = "DEL"
  
  nondelcounts <- df %>% filter(df$SVTYPE != "DEL") %>% group_by(LenCategory) %>% count()
  nondelcounts$SVTYPE = "INS"
  
  ggplot(df) +
    geom_bar(data = df %>% filter(SVTYPE != "DEL"), position = "stack", stat = "identity", aes(x = LenCategory, fill = SVTYPE, y = 1))+
    geom_bar(data = df %>% filter(SVTYPE == "DEL"), position = "stack", stat = "identity", aes(x = LenCategory, fill = SVTYPE, y = -1))+
    scale_x_discrete(labels=labellist) +
    xlab("Length") +
    ylab("Count") +
    labs(title = paste("Indels by Length (", caller, ")")) +
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

indir <- "/home/mkirsche/jasmine_data/figures/investigation/crosstech"

mergedtsv <- paste(indir, "/merged.tsv", sep = "")

df <- read.table(mergedtsv, sep = "\t", header = T) 
filtered <- df %>% filter(IS_SPECIFIC == 1 & IS_PRECISE == 1)
nrow(filtered)

outfile <- paste(indir, "/agreed_length.png", sep = "")
plot_length(filtered %>% filter(SUPP_VEC=="11"), "both aligners", outfile) 
outfile <- paste(indir, "/wm_only_length.png", sep = "")
plot_length(filtered %>% filter(SUPP_VEC == "1" | SUPP_VEC == "01"), "winnowmap", outfile) 
outfile <- paste(indir, "/ngmlr_only_length.png", sep = "")
plot_length(filtered %>% filter(SUPP_VEC == "10"), "ngmlr", outfile) 

merged_hiconf <- paste(indir, "/hiconf_header.tsv", sep = "")
df <- read.table(merged_hiconf, sep = "\t", header = T) 
filtered <- df %>% filter(IS_SPECIFIC == 1 & IS_PRECISE == 1)
outfile <- paste(indir, "/hiconf_all.png", sep = "")
plot_length(filtered, "all high confidence", outfile) 

nrow(filtered %>% filter(SUPP_VEC=="11"))
nrow(filtered %>% filter(SUPP_VEC=="10"))
nrow(filtered %>% filter(SUPP_VEC == "1" | SUPP_VEC == "01"))

outfile <- paste(indir, "/hiconf_both.png", sep = "")
plot_length(filtered %>% filter(SUPP_VEC=="11"), "high confidence concordant", outfile)

outfile <- paste(indir, "/hiconf_ngmlr.png", sep = "")
plot_length(filtered %>% filter(SUPP_VEC=="10"), "high confidence ngmlr only", outfile) 

outfile <- paste(indir, "/hiconf_winnowmap.png", sep = "")
plot_length(filtered %>% filter(SUPP_VEC == "1" | SUPP_VEC == "01"), "high confidence winnowmap only", outfile) 


indir <- "/home/mkirsche/jasmine_data/figures/investigation/hifi_ont"

mergedtsv <- paste(indir, "/merged.tsv", sep = "")
df <- read.table(mergedtsv, sep = "\t", header = T) 
filtered <- df %>% filter(IS_SPECIFIC == 1 & IS_PRECISE == 1)
nrow(filtered)

outfile <- paste(indir, "/both.png", sep = "")
plot_length(filtered %>% filter(SUPP_VEC=="11"), "concordant", outfile)

outfile <- paste(indir, "/ngmlr.png", sep = "")
plot_length(filtered %>% filter(SUPP_VEC=="10"), "ngmlr only", outfile) 

outfile <- paste(indir, "/winnowmap.png", sep = "")
plot_length(filtered %>% filter(SUPP_VEC == "1" | SUPP_VEC == "01"), "winnowmap only", outfile) 

merged_hiconf <- paste(indir, "/hiconf_header.tsv", sep = "")
df <- read.table(merged_hiconf, sep = "\t", header = T) 
filtered <- df %>% filter(IS_SPECIFIC == 1 & IS_PRECISE == 1)
outfile <- paste(indir, "/hiconf_all.png", sep = "")
plot_length(filtered, "all high confidence", outfile) 

outfile <- paste(indir, "/hiconf_both.png", sep = "")
plot_length(filtered %>% filter(SUPP_VEC=="11"), "high confidence concordant", outfile)

outfile <- paste(indir, "/hiconf_ont.png", sep = "")
plot_length(filtered %>% filter(SUPP_VEC=="10"), "high confidence ONT only", outfile) 

outfile <- paste(indir, "/hiconf_hifi.png", sep = "")
plot_length(filtered %>% filter(SUPP_VEC == "1" | SUPP_VEC == "01"), "high confidence HiFi only", outfile) 


indir <- "/home/mkirsche/jasmine_data/figures/investigation/chm13_hifi"

mergedtsv <- paste(indir, "/merged.tsv", sep = "")
df <- read.table(mergedtsv, sep = "\t", header = T) 
filtered <- df %>% filter(IS_SPECIFIC == 1 & IS_PRECISE == 1)
nrow(filtered)

outfile <- paste(indir, "/agreed_length.png", sep = "")
plot_length(filtered %>% filter(SUPP_VEC=="11"), "both techs", outfile) 
outfile <- paste(indir, "/hifi_only_length.png", sep = "")
plot_length(filtered %>% filter(SUPP_VEC == "1" | SUPP_VEC == "01"), "hifi", outfile) 
outfile <- paste(indir, "/ont_only_length.png", sep = "")
plot_length(filtered %>% filter(SUPP_VEC == "10"), "ONT", outfile) 

outfile <- paste(indir, "/hiconf_both.png", sep = "")
plot_length(filtered %>% filter(SUPP_VEC=="11"), "high confidence concordant", outfile)

outfile <- paste(indir, "/hiconf_ngmlr.png", sep = "")
plot_length(filtered %>% filter(SUPP_VEC=="10"), "high confidence ngmlr only", outfile) 

outfile <- paste(indir, "/hiconf_winnowmap.png", sep = "")
plot_length(filtered %>% filter(SUPP_VEC == "1" | SUPP_VEC == "01"), "high confidence winnowmap only", outfile) 

merged_hiconf <- paste(indir, "/hiconf_header.tsv", sep = "")
df <- read.table(merged_hiconf, sep = "\t", header = T) 
filtered <- df %>% filter(IS_SPECIFIC == 1 & IS_PRECISE == 1)
outfile <- paste(indir, "/hiconf_all.png", sep = "")
plot_length(filtered, "all high confidence", outfile) 

outfile <- paste(indir, "/hiconf_both.png", sep = "")
plot_length(filtered %>% filter(SUPP_VEC=="11"), "high confidence concordant", outfile)

outfile <- paste(indir, "/hiconf_ont.png", sep = "")
plot_length(filtered %>% filter(SUPP_VEC=="10"), "high confidence ONT only", outfile) 

outfile <- paste(indir, "/hiconf_hifi.png", sep = "")
plot_length(filtered %>% filter(SUPP_VEC == "1" | SUPP_VEC == "01"), "high confidence HiFi only", outfile) 

library(circlize)

library(RCircos)
library(dplyr)
library(ggplot2)
library(Cairo)


plot_circos <- function(fn, outdir, prefix, supp)
{
  df <- read.table(fn, header = TRUE, sep = "\t")
  nrow(df)
  df <- df %>% filter(IS_PRECISE == "1" & IS_SPECIFIC == "1" & (SUPP_VEC == supp))
  nrow(df)
  colnames(df)
  unique(df$CHROM)
  goodlist <- paste('chr', c('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', 
                             '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', "X", "Y"), sep = "")
  
  goodlist
  df <- df %>% filter(CHROM %in% goodlist)
  df$CHROM <- droplevels(df$CHROM)
  df
  
  binsize <- 1000000
  df$STARTBIN <- floor(df$POS / binsize) * binsize
  groups <- df %>% group_by(CHROM, STARTBIN, SVTYPE) %>% tally()
  groups$CHROM <- droplevels(groups$CHROM)
  groups
  
  outfile <- paste(outdir, prefix, "_flatcircos.png", sep = "")
  ggplot(groups %>% filter(CHROM == "chr1"), aes(x = STARTBIN, y = n, fill = SVTYPE)) + geom_bar( stat = "identity")
  ggsave(outfile, width = 16, height = 6)
  
  insgroups <- groups %>% filter(SVTYPE == "INS")
  delgroups <- groups %>% filter(SVTYPE == "DEL")
  
  outfile <- paste(outdir, prefix, "_circos.svg", sep = "")
  #png(outfile, width = 4096, height = 4096)
  Cairo(file=outfile, 
        type="svg",
        units="in", 
        width=10, 
        height=10, 
        pointsize=12, 
        dpi=72)
  circos.initializeWithIdeogram()
  circos.track(insgroups$CHROM, x = insgroups$STARTBIN, y = insgroups$n, panel.fun = function(x, y) {
    circos.barplot(y, x, col = 3, bar_width = binsize, lwd = 0.00001, lty = 0)
  })
  circos.track(delgroups$CHROM, x = delgroups$STARTBIN, y = delgroups$n, panel.fun = function(x, y) {
    circos.barplot(y, x, col = 2, bar_width = binsize, lwd = 0.00001, lty = 0)
  })
  dev.off()
}

outdir <- "/home/mkirsche/jasmine_data/figures/investigation/"

fn <- "/home/mkirsche/jasmine_data/figures/investigation/chm13_hifi/merged.tsv"
prefix <- "chm13_wmonly"
plot_circos(fn, outdir, prefix, "1")

fn <- "/home/mkirsche/jasmine_data/figures/investigation/crosstech/merged.tsv"
prefix <- "grch38_wmonly"
plot_circos(fn, outdir, prefix, "1")

fn <- "/home/mkirsche/jasmine_data/figures/investigation/chm13_hifi/merged.tsv"
prefix <- "chm13_ngmlronly"
plot_circos(fn, outdir, prefix, "10")

fn <- "/home/mkirsche/jasmine_data/figures/investigation/crosstech/merged.tsv"
prefix <- "grch38_ngmlronly"
plot_circos(fn, outdir, prefix, "10")

fn <- "/home/mkirsche/jasmine_data/figures/investigation/chm13_hifi/merged.tsv"
prefix <- "chm13_shared"
plot_circos(fn, outdir, prefix, "11")

fn <- "/home/mkirsche/jasmine_data/figures/investigation/crosstech/merged.tsv"
prefix <- "grch38_shared"
plot_circos(fn, outdir, prefix, "11")
