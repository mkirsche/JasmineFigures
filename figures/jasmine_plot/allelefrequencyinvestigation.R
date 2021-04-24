library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)
library(reshape2)
library(gridExtra)
library(rlist)

plot_allele_frequencies <- function(variants, outfile) {
  variants$SUPP_INT = strtoi(variants$SUPP)
  ggplot(variants, aes(x = SUPP_INT, y = 1, fill = SVPOP_SINGLE)) +
    geom_bar(position = "stack", stat = "identity") +
    ggtitle("Population Allele Frequencies") +
    #scale_fill_viridis_c(option = "plasma") +
    theme_bw() +
    theme(axis.text.x = element_text(size = 16),
          axis.text.y = element_text(size = 16),
          axis.title.x = element_text(size = 18),
          axis.title.y = element_text(size = 18),
          #legend.title = element_blank(),
          #legend.text = element_blank(),
          #legend.position = "none",
          plot.title = element_text(hjust = 0.5, size = 20)) +
    xlab("Allele Frequency") +
    ylab("Number of variants") +
    
    ggsave(outfile, device = "png", height = 8, width = 8)
}

jasminefile <- "/home/mkirsche/eclipse-workspace/AlleleFrequencyInvestigation/jasmine_svpopanno.tsv"
outfile <- "/home/mkirsche/jasmine_data/figures/figure5/jasmine_vs_svpop_af.png"
variants = read.table(jasminefile, comment.char = "#", sep = "\t", header = T, stringsAsFactors=FALSE, colClasses = 'character')
variants <- variants %>% filter(SPECIFIC_FLAG == "1" & PRECISE_FLAG == "1")
unique(variants$SVPOP_SINGLE)
plot_allele_frequencies(variants, outfile)

jasminetriosupportfile <- "/home/mkirsche/jasmine_data/figures/figure2/jasmine_support.txt"
outfile <- "/home/mkirsche/jasmine_data/figures/figure2/jasmine_re_af.png"
df = read.table(jasminetriosupportfile, comment.char = "#", sep = "\t", header = T)
df$RE_NUM <- as.numeric(df$RE)
df$AF_NUM <- as.character(df$AF)
hg002 <- df %>% filter(SAMPLE == 'HG002')
hg002 %>% group_by(RE) %>% tally()
ggplot(df %>% filter(RE_NUM < 30), aes(x = RE_NUM, y = 1, fill = AF_NUM)) + geom_bar(position = "stack", stat = "identity") +
  facet_grid(~SAMPLE) +
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        #legend.title = element_blank(),
        #legend.text = element_blank(),
        #legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 20)) +
  xlab("Read Support") +
  ylab("Number of variants")
ggsave(outfile, device = "png", height = 12, width = 16)

jasminetriosupportfile <- "/home/mkirsche/jasmine_data/figures/figure3/jasmine_crosstech_support.txt"
outfile <- "/home/mkirsche/jasmine_data/figures/figure3/jasmine_crosstech_re_af.png"
df = read.table(jasminetriosupportfile, comment.char = "#", sep = "\t", header = T)
df$RE_NUM <- as.numeric(df$RE)
df$AF_NUM <- as.character(df$AF)
df$SAMPLE_TECH <- paste(df$SAMPLE, df$TECH, sep = '_')
hg002 <- df %>% filter(SAMPLE == 'HG002')
hg002 %>% group_by(RE) %>% tally()
ggplot(df %>% filter(RE_NUM < 30), aes(x = RE_NUM, y = 1, fill = AF_NUM)) + geom_bar(position = "stack", stat = "identity") +
  facet_grid(~SAMPLE_TECH) +
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        #legend.title = element_blank(),
        #legend.text = element_blank(),
        #legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 20)) +
  xlab("Read Support") +
  ylab("Number of variants")
ggsave(outfile, device = "png", height = 12, width = 16)

# Read in sample info
fn <- "/home/mkirsche/jasmine_data/figures/figure5/filelist.txt"
samples <- read.table(fn, sep = "\t", header = T, colClasses = 'character')
samples$SAMPLE
samples$Label <- paste(samples$SAMPLE, "_", samples$TECH, sep = "")

jasminepopulationsupportfile <- "/home/mkirsche/jasmine_data/figures/figure5/all4_pop_jasmine_support_filtered.txt"
jasminepopulationsupportfile <- "/home/mkirsche/jasmine_data/figures/figure5/all4_pop_jasmine_onlyspecprec_support.txt"
outfile <- "/home/mkirsche/jasmine_data/figures/figure5/jasmine_population_re_af.png"
df = read.table(jasminepopulationsupportfile, comment.char = "#", sep = "\t", header = T)
unique(df$SAMPLE)
df$RE_NUM <- as.numeric(df$RE)
df <- df %>% filter(RE_NUM <= 30)
df$AF_NUM <- as.factor(as.character(df$AF))
levels(df$AF_NUM)
sorted_labels <- paste(sort(as.integer(levels(df$AF_NUM))))
sorted_labels
df$AF_NUM <- factor(df$AF_NUM, levels = sorted_labels)
#df <- df[order(df$AF), ]

#df <- merge(df, samples, by = SAMPLE)
for (x in unique(df$SAMPLE)) {
  print(x)
  ggplot(df %>% filter(SAMPLE == x), aes(x = RE_NUM, y = 1, fill = AF_NUM)) + geom_bar(position = "stack", stat = "identity") +
    facet_grid(~SAMPLE) +
    theme(axis.text.x = element_text(size = 8),
          axis.text.y = element_text(size = 8),
          axis.title.x = element_text(size = 18),
          axis.title.y = element_text(size = 18),
          #legend.title = element_blank(),
          #legend.text = element_blank(),
          #legend.position = "none",
          plot.title = element_text(hjust = 0.5, size = 20)) +
    xlab("Read Support") +
    ylab("Number of variants")
  ggsave(paste("/home/mkirsche/jasmine_data/figures/figure5/jasmine_population_re_af_" , x, "_.png" , sep = ""), device = "png", height = 12, width = 12)  
}

ggplot(df, aes(x = RE_NUM, y = 1, fill = AF_NUM)) + geom_bar(position = "stack", stat = "identity") +
  facet_grid(~SAMPLE) +
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        #legend.title = element_blank(),
        #legend.text = element_blank(),
        #legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 20)) +
  xlab("Read Support") +
  ylab("Number of variants")
ggsave(paste("/home/mkirsche/jasmine_data/figures/figure5/jasmine_population_re_af_" , "all", "_.png" , sep = ""), device = "png", height = 20, width = 20)
