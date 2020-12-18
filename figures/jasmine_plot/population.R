library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)
library(reshape2)
library(gridExtra)
library(rlist)

plot_variant_discovery <- function(variants, outfile) {
  variants$SUPP_VEC_LENGTH <- nchar(variants$SUPP_VEC)
  maxlength=max(variants$SUPP_VEC_LENGTH)
  variants$SAMPLE_DISCOVERED <- maxlength + 1 - variants$SUPP_VEC_LENGTH
  counts <- variants %>% count(SAMPLE_DISCOVERED)
  counts <- counts %>% mutate(cumsum = cumsum(n))
  counts
  ggplot(counts, aes(x = SAMPLE_DISCOVERED, y = cumsum, fill = SAMPLE_DISCOVERED)) +
    geom_bar(stat = "sum") +
  ggtitle("Variant Discovery") +
    scale_fill_viridis_c() +
    theme_bw() +
    theme(axis.text.x = element_text(size = 16),
          axis.text.y = element_text(size = 16),
          axis.title.x = element_text(size = 18),
          axis.title.y = element_text(size = 18),
          legend.title = element_blank(),
          legend.text = element_blank(),
          legend.position = "none",
          plot.title = element_text(hjust = 0.5, size = 20)) +
    xlab("Samples") +
    ylab("Variants Discovered")
  ggsave(outfile, device = "png", height = 8, width = 8)
}

plot_allele_frequencies <- function(variants, outfile) {
  ggplot(variants, aes(x = SUPP, y = 1, fill = SUPP)) +
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

plot_length <- function(df, outfile) 
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
    labs(title = "Indels by Length") +
    theme_bw() +
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

plot_shuffle_count <- function(counttablefile, outfile)
{
  shuffle_data <- read.table(counttablefile, sep = "\t", header = TRUE)
  ggplot(shuffle_data, aes(COUNT, MERGER)) + geom_boxplot() +
    xlab("Number of variants") +
    ylab("Merging software")
  ggsave(outfile, width= 8, height = 4)
}


get_merges <- function(tablefile, merger) {
  df <- read.table(tablefile, sep = "\t", header = TRUE)
  df <- df %>% filter(SPECIFIC_FLAG == 1 & PRECISE_FLAG == 1)
  df$startspan = df$MAX_START - df$MIN_START
  df$endspan = df$MAX_END - df$MIN_END
  df$MaxRange = pmax(df$startspan, df$endspan)
  df$Range <- log(df$MaxRange + 1, base = 2)
  multi <- df %>% filter(NUMVARS > 1)
  multi$Merger <- merger
  
  return(multi %>% select(Range, Merger))
}

plot_range <- function(jasminefile, survivorfile, svtoolsfile, svimmerfile, outfile) {
  names = c("Range", "Merger")
  rangedf <- data.frame(setNames(rep(list(NA), length(names)), names))
  rangedf$Range = c()
  rangedf$Merger = c()

  jasminedf <- get_merges(jasminefile, 'Jasmine')
  rangedf <- rbind(rangedf, jasminedf)
  
  survivordf <- get_merges(survivorfile, 'SURVIVOR')
  rangedf <- rbind(rangedf, survivordf)
  
  svtoolsdf <- get_merges(svtoolsfile, 'svtools')
  rangedf <- rbind(rangedf, svtoolsdf)
  
  svimmerdf <- get_merges(svimmerfile, 'svimmer')
  rangedf <- rbind(rangedf, svimmerdf)
  
  ggplot(rangedf, aes(Range, Merger)) + geom_boxplot() +
    xlab("Breakpoint Range") +
    ylab("Merging software")
  ggsave(outfile, width= 8, height = 8)
}

merged_table_file <- "/home/mkirsche/jasmine_data/figures/figure5/population.merged.tsv"
variants = read.table(merged_table_file, comment.char = "#", sep = "\t", header = T, stringsAsFactors=FALSE)
variants <- variants %>% filter(IS_SPECIFIC == "1" & IS_PRECISE == "1")

outfile <- "/home/mkirsche/jasmine_data/figures/figure5/variantdiscovery.png"
plot_variant_discovery(variants, outfile)

outfile <- "/home/mkirsche/jasmine_data/figures/figure5/allelefreqs.png"
plot_allele_frequencies(variants, outfile)

outfile <- "/home/mkirsche/jasmine_data/figures/figure5/variantlengths.png"
plot_length(variants, outfile)

counttablefile <- "/home/mkirsche/jasmine_data/figures/figure5/population.counts.txt"
outfile <- "/home/mkirsche/jasmine_data/figures/figure5/shuffle_boxplot.png"
plot_shuffle_count(counttablefile, outfile)

jasminefile <- "/home/mkirsche/jasmine_data/figures/figure5/population.jasminetable.txt"
survivorfile <- "/home/mkirsche/jasmine_data/figures/figure5/population.survivortable.txt"

outfile <- "/home/mkirsche/jasmine_data/figures/figure5/rangeplot.png"
plot_range(jasminefile, survivorfile, jasminefile, survivorfile, outfile)
