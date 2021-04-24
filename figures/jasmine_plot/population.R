library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)
library(reshape2)
library(gridExtra)
library(rlist)
library(svglite)

plot_variant_discovery <- function(variants, outfile) {
  variants$TRIMMED_SUPP_VEC = str_remove(variants$SUPP_VEC, "^0+")
  variants$SUPP_VEC_LENGTH <- nchar(variants$TRIMMED_SUPP_VEC)
  unique(variants$SUPP_VEC_LENGTH)
  maxlength=max(variants$SUPP_VEC_LENGTH)
  maxlength
  variants$SAMPLE_DISCOVERED <- maxlength + 1 - variants$SUPP_VEC_LENGTH
  variants$SAMPLE_DISCOVERED
  counts <- variants %>% count(SAMPLE_DISCOVERED)
  counts <- counts %>% mutate(cumsum = cumsum(n))
  counts
  ggplot(counts, aes(x = SAMPLE_DISCOVERED, y = cumsum, fill = SAMPLE_DISCOVERED)) +
    geom_bar(stat = "sum", fill = "#2093c3") +
  ggtitle("Variant Discovery") +
    theme(axis.text.x = element_text(size = 16),
          axis.text.y = element_text(size = 16),
          axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20),
          legend.title = element_blank(),
          legend.text = element_blank(),
          legend.position = "none",
          plot.title = element_text(hjust = 0.5, size = 24)) +
    xlab("Samples") +
    ylab("Variants Discovered")
  ggsave(outfile, height = 8, width = 8)
}

plot_allele_frequencies <- function(variants, outfile) {
  variants$SUPP_INT = strtoi(variants$SUPP)
  ggplot(variants, aes(x = SUPP_INT, fill = SUPP_INT)) +
    geom_bar(position = "stack", stat = "count", fill = "#2093c3") +
    ggtitle("Population Allele Frequencies") +
    theme(axis.text.x = element_text(size = 16),
          axis.text.y = element_text(size = 16),
          axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20),
          legend.title = element_blank(),
          legend.text = element_blank(),
          legend.position = "none",
          plot.title = element_text(hjust = 0.5, size = 24)) +
    xlab("Allele Frequency") +
    ylab("Number of variants") +

  ggsave(outfile, height = 8, width = 8)
}

plot_length <- function(df, outfile, legendx, legendy) 
{
  df$LenCategory = "TRA"
  
  df$ABSLEN = abs(strtoi(df$SVLEN))
  
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
  
  colorpalette <- c(DEL = '#CC6677',DUP = '#332288',INS = '#44AA99',TRA = '#88CCEE',INV = "#ddaa33")#'#DDCC77')
  
  summarized <- df %>% group_by(SVTYPE, LenCategory) %>% summarise(counts=n())
  
  ggplot(summarized) +
    geom_bar(data = summarized %>% filter(SVTYPE != "DEL"), position = "stack", stat = "identity", aes(x = LenCategory, fill = SVTYPE, y = counts))+
    geom_bar(data = summarized %>% filter(SVTYPE == "DEL"), position = "stack", stat = "identity", aes(x = LenCategory, fill = SVTYPE, y = -counts))+
    scale_x_discrete(labels=labellist) +
    xlab("Length") +
    ylab("Count") +
    labs(title = "Population SV Size Distribution") +
    theme(plot.title = element_text(size = 24, hjust = 0.5),
          axis.text.x = element_text(size = 18, angle = 30, margin = margin(t = 17, r = 0, b = 0, l = 0)),
          axis.text.y = element_text(size = 18),
          axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20),
          legend.text = element_text(size = 18),
          legend.title = element_text(size = 20),
          text = element_text(size = 20),
          legend.position = c(legendx, legendy),
    ) +
    scale_fill_manual(name = "SVTYPE", values = colorpalette) +
    geom_text(data = nondelcounts, aes(x = LenCategory, y=n, label=n), position=position_dodge(width=0.9), vjust=-0.75) +
    geom_text(data = delcounts, aes(x = LenCategory, y=-1*n, label=n), position=position_dodge(width=0.9), vjust=1.5)
  
  ggsave(outfile, width= 10, height = 8)
}

plot_shuffle_count <- function(counttablefile, outfile)
{
  shuffle_data <- read.table(counttablefile, sep = "\t", header = TRUE)
  shuffle_data$order <- "1"
  shuffle_data$order <- ifelse(shuffle_data$MERGER == "Jasmine", "6", shuffle_data$order)
  shuffle_data$order <- ifelse(shuffle_data$MERGER == "dbsvmerge", "5", shuffle_data$order)
  shuffle_data$order <- ifelse(shuffle_data$MERGER == "svpop", "4", shuffle_data$order)
  shuffle_data$order <- ifelse(shuffle_data$MERGER == "SURVIVOR", "3", shuffle_data$order)
  shuffle_data$order <- ifelse(shuffle_data$MERGER == "svimmer", "2", shuffle_data$order)
  ggplot() + geom_boxplot(data = shuffle_data %>% filter(MERGER == "Jasmine" | MERGER == "svtools" | MERGER == "svimmer"), aes(x = COUNT, y = order)) + 
    geom_violin(scale = "width", draw_quantiles = c(.5), data = shuffle_data %>% filter(MERGER == "SURVIVOR" | MERGER == "svpop" | MERGER == "dbsvmerge"), aes(x = COUNT, y = order, fill = MERGER)) +
    scale_y_discrete(limits=c("1", "2", "3", "4", "5", "6"),labels=c("svtools", "svimmer", "SURVIVOR", "svpop", "dbsvmerge", "Jasmine")) +
    xlab("Number of variants across 100 permutations") +
    ylab("Merging software") + 
    theme(legend.title = element_blank(), legend.position = "none", legend.text = element_blank(), axis.title.y = element_blank(), axis.text = element_text(size = 18), axis.title.x = element_text(size = 18)) + 
    scale_fill_manual(values = c(SURVIVOR = "#33bbee", svimmer = "#ee3377", svtools = "#ee7733", Jasmine = "#009988", svpop = "#cc3311", dbsvmerge = "#0077bb"))
  ggsave(outfile, width= 8, height = 4)
}

get_merges <- function(tablefile, merger) {
  df <- read.table(tablefile, sep = "\t", header = TRUE)
  df <- df %>% filter(SPECIFIC_FLAG == 1 & PRECISE_FLAG == 1)
  df <- df %>% filter(abs(LEN) >= 50 | TYPE == 'TRA')
  df$startspan = df$MAX_START - df$MIN_START
  df$endspan = df$MAX_END - df$MIN_END
  df$MaxRange = pmax(df$startspan, df$endspan)
  df$Range <- log(df$MaxRange + 1, base = 2)
  multi <- df %>% filter(NUMVARS > 1)
  sum(df$NUMVARS)
  nrow(df)
  nrow(multi)
  multi$Merger <- merger
  
  return(multi %>% select(Range, Merger))
}

plot_range <- function(jasminefile, survivorfile, svtoolsfile, svimmerfile, svpopfile, dbsvmergefile, outfile) {
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
  
  svpopdf <- get_merges(svpopfile, 'svpop')
  rangedf <- rbind(rangedf, svpopdf)
  
  dbsvmergedf <- get_merges(dbsvmergefile, 'dbsvmerge')
  rangedf <- rbind(rangedf, dbsvmergedf)
  
  rangedf$order <- "1"
  rangedf$order <- ifelse(rangedf$Merger == "Jasmine", "6", rangedf$order)
  rangedf$order <- ifelse(rangedf$Merger == "dbsvmerge", "5", rangedf$order)
  rangedf$order <- ifelse(rangedf$Merger == "svpop", "4", rangedf$order)
  rangedf$order <- ifelse(rangedf$Merger == "SURVIVOR", "3", rangedf$order)
  rangedf$order <- ifelse(rangedf$Merger == "svimmer", "2", rangedf$order)
  
  ggplot(rangedf, aes(x = Range, y = order, fill = Merger)) + geom_violin(draw_quantiles = c(.5), scale = "width") +
    xlab("Breakpoint Range (log2)") +
    ylab("Merging software") + theme(legend.title = element_blank(), legend.position = "none", legend.text = element_blank(), axis.title.y = element_blank(), axis.text = element_text(size = 18), axis.title.x = element_text(size = 18))  +
    scale_y_discrete(limits=c("1", "2", "3", "4", "5", "6"),labels=c("svtools", "svimmer", "SURVIVOR", "svpop", "dbsvmerge", "Jasmine")) +
    scale_fill_manual(values = c(SURVIVOR = "#33bbee", svimmer = "#ee3377", svtools = "#ee7733", Jasmine = "#009988", svpop = "#cc3311", dbsvmerge = "#0077bb"))
  ggsave(outfile, width= 8, height = 4)
}

violin_plot <- function(fn, merger) {
  variants = read.table(fn, comment.char = "#", sep = "\t", header = T, stringsAsFactors=FALSE, colClasses = 'character')
  variants <- variants %>% filter(SPECIFIC_FLAG == "1" & PRECISE_FLAG == "1")
  variants$NUM_SUPP <- as.factor(as.numeric(variants$SUPP))
  variants$NUM_SPECIFIC <- as.numeric(variants$SPECIFIC_COUNT)
  
  fewplot <- ggplot(variants %>% filter(as.numeric(SUPP) <= 10), aes(x=NUM_SUPP, y=NUM_SPECIFIC, fill = NUM_SUPP)) + 
    geom_violin(draw_quantiles = c(.5)) + ylab('Number of High Confidence Samples') +
    theme(plot.title = element_text(size = 24, hjust = 0.5),
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 16),
          legend.text = element_blank(),
          text = element_text(size = 16),
          legend.position = "none"
    )
  outfile <- paste("/home/mkirsche/jasmine_data/figures/figure5/", merger, "violin_few.png", sep = "")
  ggsave(outfile, width = 20, height = 8)
  
  
  manyplot <- ggplot(variants %>% filter(as.numeric(SUPP) > 10), aes(x=NUM_SUPP, y=NUM_SPECIFIC, fill = NUM_SUPP)) + 
    geom_violin() + xlab('Number of Samples') + ylab('Number of High Confidence Samples') +
    theme(plot.title = element_text(size = 24, hjust = 0.5),
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          legend.text = element_blank(),
          text = element_text(size = 16),
          legend.position = "none"
    )
  outfile <- paste("/home/mkirsche/jasmine_data/figures/figure5/", merger, "violin_many.png", sep = "")
  ggsave(outfile, width = 20, height = 8)
  
  
  outfile <- paste("/home/mkirsche/jasmine_data/figures/figure5/", merger, "violin_all.png", sep = "")
  g <- arrangeGrob(fewplot, manyplot, nrow=2)
  ggsave(outfile, g, width = 16, height = 10)
}

plot_ism <- function(indir, outfile) {
  jasminefile <- paste(indir, "/jasmineism.txt", sep = "")
  survivorfile <- paste(indir, "/survivorism.txt", sep = "")
  svtoolsfile <- paste(indir, "/svtoolsism.txt", sep = "")
  svimmerfile <- paste(indir, "/svimmerism.txt", sep = "")
  dbsvmergefile <- paste(indir, "/dbsvmergeism.txt", sep = "")
  svpopfile <- paste(indir, "/svpopism.txt", sep = "")
  df <- data.frame(matrix(nrow = 0, ncol = 2))
  colnames(df) <- c("Intrasample", "Merger")
  
  tmpdf <- read.table(jasminefile, header = F, sep = "\t")
  colnames(tmpdf) <- "Intrasample"
  tmpdf$Merger <- "Jasmine"
  df <- rbind(df, tmpdf)
  
  tmpdf <- read.table(survivorfile, header = F, sep = "\t")
  colnames(tmpdf) <- "Intrasample"
  tmpdf$Merger <- "SURVIVOR"
  df <- rbind(df, tmpdf)
  
  tmpdf <- read.table(svtoolsfile, header = F, sep = "\t")
  colnames(tmpdf) <- "Intrasample"
  tmpdf$Merger <- "svtools"
  df <- rbind(df, tmpdf)
  
  tmpdf <- read.table(svimmerfile, header = F, sep = "\t")
  colnames(tmpdf) <- "Intrasample"
  tmpdf$Merger <- "svimmer"
  df <- rbind(df, tmpdf)
  
  tmpdf <- read.table(svpopfile, header = F, sep = "\t")
  colnames(tmpdf) <- "Intrasample"
  tmpdf$Merger <- "svpop"
  df <- rbind(df, tmpdf)
  
  tmpdf <- read.table(dbsvmergefile, header = F, sep = "\t")
  colnames(tmpdf) <- "Intrasample"
  tmpdf$Merger <- "dbsvmerge"
  df <- rbind(df, tmpdf)
  
  df$order <- "1"
  df$order <- ifelse(df$Merger == "Jasmine", "6", df$order)
  df$order <- ifelse(df$Merger == "dbsvmerge", "5", df$order)
  df$order <- ifelse(df$Merger == "svpop", "4", df$order)
  df$order <- ifelse(df$Merger == "SURVIVOR", "3", df$order)
  df$order <- ifelse(df$Merger == "svimmer", "2", df$order)
  
  ggplot() +  
    geom_boxplot(data = df %>% filter(Merger == "Jasmine" | Merger == "svpop" | Merger == "dbsvmerge"), aes(x = log2(Intrasample + 1), y = order)) +
    geom_violin(scale = "width", draw_quantiles = c(.5), data = df %>% filter(Merger == "SURVIVOR" | Merger == "svimmer" | Merger == "svtools"), aes(x = log2(Intrasample + 1), y = order, fill = Merger)) + 
    theme(legend.title = element_blank(), legend.position = "none", legend.text = element_blank(), axis.title.y = element_blank(), axis.text = element_text(size = 18), axis.title.x = element_text(size = 18)) + xlab("Intrasample Merges per Variant (log2)") + ylab("Merging software") +
    scale_y_discrete(limits=c("1", "2", "3", "4", "5", "6"),labels=c("svtools", "svimmer", "SURVIVOR", "svpop", "dbsvmerge", "Jasmine")) +
    scale_fill_manual(values = c(SURVIVOR = "#33bbee", svimmer = "#ee3377", svtools = "#ee7733", Jasmine = "#009988", svpop = "#cc3311", dbsvmerge = "#0077bb"))
  ggsave(outfile, height = 4, width = 8)
} 
indir <- "/home/mkirsche/jasmine_data/figures/figure5"
outfile <- "/home/mkirsche/jasmine_data/figures/figure5/ism_counts.png"
plot_ism(indir, outfile)
outfile <- "/home/mkirsche/jasmine_data/figures/figure5/ism_counts.svg"
plot_ism(indir, outfile)

merged_table_file <- "/home/mkirsche/jasmine_data/figures/figure5/population_final_specprec.merged.tsv"
variants = read.table(merged_table_file, comment.char = "#", sep = "\t", header = T, stringsAsFactors=FALSE, colClasses = 'character')
variants <- variants %>% filter(IS_SPECIFIC == "1" & IS_PRECISE == "1")

outfile <- "/home/mkirsche/jasmine_data/figures/figure5/variantdiscovery.png"
plot_variant_discovery(variants, outfile)
outfile <- "/home/mkirsche/jasmine_data/figures/figure5/variantdiscovery.svg"
plot_variant_discovery(variants, outfile)

outfile <- "/home/mkirsche/jasmine_data/figures/figure5/allelefreqs.png"
plot_allele_frequencies(variants, outfile)
outfile <- "/home/mkirsche/jasmine_data/figures/figure5/allelefreqs.svg"
plot_allele_frequencies(variants, outfile)

outfile <- "/home/mkirsche/jasmine_data/figures/figure5/variantlengths.png"
plot_length(variants, outfile, .92, .84)
outfile <- "/home/mkirsche/jasmine_data/figures/figure5/variantlengths.svg"
plot_length(variants, outfile, .92, .84)

counttablefile <- "/home/mkirsche/jasmine_data/figures/figure5/population.newcounts.txt"
outfile <- "/home/mkirsche/jasmine_data/figures/figure5/shuffle_boxplot.png"
plot_shuffle_count(counttablefile, outfile)
outfile <- "/home/mkirsche/jasmine_data/figures/figure5/shuffle_boxplot.svg"
plot_shuffle_count(counttablefile, outfile)

#jasminefile <- "/home/mkirsche/jasmine_data/figures/figure5/population.jasminetable.txt"
#survivorfile <- "/home/mkirsche/jasmine_data/figures/figure5/population.survivortable.txt"

jasminefile <- "/home/mkirsche/jasmine_data/figures/figure5/jasmine_pop_final_specprec.txt"
alljasminefile <- "/home/mkirsche/jasmine_data/figures/figure5/all4_pop.jasmine_augmented.txt"
allsvpopfile <- "/home/mkirsche/jasmine_data/figures/figure5/all4_pop.svpop_augmented.txt"
survivorfile <- "/home/mkirsche/jasmine_data/figures/figure5/all4_pop.survivor.specprec.txt"
svtoolsfile <- "/home/mkirsche/jasmine_data/figures/figure5/all4_pop.svtools.specprec.txt"
svimmerfile <- "/home/mkirsche/jasmine_data/figures/figure5/all4_pop.svimmer.specprec.txt"
svpopfile <- "/home/mkirsche/jasmine_data/figures/figure5/all4_pop.svpop.specprec.txt"
outfile <- "/home/mkirsche/jasmine_data/figures/figure5/rangeplot.png"
jasminespec10file <- "/home/mkirsche/jasmine_data/figures/allelefreq/spec10_filtered.txt"
jasminespec5file <- "/home/mkirsche/jasmine_data/figures/allelefreq/spec5_filtered.txt"
jasminehififile <- "/home/mkirsche/jasmine_data/figures/figure5/all4_pop_jasminehifi_augmented.txt"
jasminenonhififivefile <- "/home/mkirsche/jasmine_data/figures/figure5/all4_pop_nonhififiveplus_augmented.txt"
jasminenonhifithreefile <- "/home/mkirsche/jasmine_data/figures/figure5/all4_pop_nonhifithreeplus_augmented.txt"
jasminenonhififourpercentfile <- "/home/mkirsche/jasmine_data/figures/figure5/all4_pop_nonhififourpercentplus_augmented.txt"
jasminebigfile <- "/home/mkirsche/jasmine_data/figures/figure5/filter50_augmented.txt"
jasminelineardistfile <- "/home/mkirsche/jasmine_data/figures/figure5/lineardist_augmented_specprec.txt"
jasmine50_200lineardistfile <- "/home/mkirsche/jasmine_data/figures/figure5/lineardist50_200_specprec.txt"
jasmine100_200lineardistfile <- "/home/mkirsche/jasmine_data/figures/figure5/lineardist100_200_specprec.txt"
dbsvmergefile <- "/home/mkirsche/jasmine_data/figures/figure5/all4_pop.dbsvmerge.specprec.txt"
jasmine20_200lineardistfile <- "/home/mkirsche/jasmine_data/figures/figure5/lineardist20_200_specprec.txt"
jasmine50_50lineardistfile <- "/home/mkirsche/jasmine_data/figures/figure5/lineardist50_50_specprec.txt"

plot_range(jasminefile, survivorfile, svtoolsfile, svimmerfile, svpopfile, dbsvmergefile, outfile)
outfile <- "/home/mkirsche/jasmine_data/figures/figure5/rangeplot.svg"

plot_range(jasminefile, survivorfile, svtoolsfile, svimmerfile, svpopfile, dbsvmergefile, outfile)



outfile <- "/home/mkirsche/jasmine_data/figures/figure5/jasmineallelefreqs.png"
variants = read.table(jasminefile, comment.char = "#", sep = "\t", header = T, stringsAsFactors=FALSE, colClasses = 'character')
variants <- variants %>% filter(SPECIFIC_FLAG == "1" & PRECISE_FLAG == "1")
plot_allele_frequencies(variants, outfile)
outfile <- "/home/mkirsche/jasmine_data/figures/figure5/jasmineallelefreqs.svg"
plot_allele_frequencies(variants, outfile)


outfile <- "/home/mkirsche/jasmine_data/figures/figure5/jasminehifiaf.png"
variants = read.table(jasminehififile, comment.char = "#", sep = "\t", header = T, stringsAsFactors=FALSE, colClasses = 'character')
variants <- variants %>% filter(SPECIFIC_FLAG == "1" & PRECISE_FLAG == "1")
plot_allele_frequencies(variants, outfile)

#outfile <- "/home/mkirsche/jasmine_data/figures/figure5/all_jasmineallelefreqs.png"
#variants = read.table(alljasminefile, comment.char = "#", sep = "\t", header = T, stringsAsFactors=FALSE, colClasses = 'character')
#plot_allele_frequencies(variants, outfile)

#outfile <- "/home/mkirsche/jasmine_data/figures/figure5/all_svpopallelefreqs.png"
#variants = read.table(allsvpopfile, comment.char = "#", sep = "\t", header = T, stringsAsFactors=FALSE, colClasses = 'character')
#plot_allele_frequencies(variants, outfile)

outfile <- "/home/mkirsche/jasmine_data/figures/figure5/jasmineallelefreqs_long.png"
variants = read.table(jasminefile, comment.char = "#", sep = "\t", header = T, stringsAsFactors=FALSE, colClasses = 'character')
variants <- variants %>% filter(SPECIFIC_FLAG == "1" & PRECISE_FLAG == "1")
variants <- variants%>% filter(abs(as.numeric(LEN)) >= 50 | TYPE == "TRA")
plot_allele_frequencies(variants, outfile)
outfile <- "/home/mkirsche/jasmine_data/figures/figure5/jasmineallelefreqs_long.svg"
plot_allele_frequencies(variants, outfile)

outfile <- "/home/mkirsche/jasmine_data/figures/figure5/survivorallelefreqs.png"
variants = read.table(survivorfile, comment.char = "#", sep = "\t", header = T, stringsAsFactors=FALSE, colClasses = 'character')
variants <- variants %>% filter(SPECIFIC_FLAG == "1" & PRECISE_FLAG == "1")
plot_allele_frequencies(variants, outfile)
outfile <- "/home/mkirsche/jasmine_data/figures/figure5/survivorallelefreqs.svg"
plot_allele_frequencies(variants, outfile)
variants <- variants%>% filter(abs(as.numeric(LEN)) >= 50 | TYPE == "TRA")
outfile <- "/home/mkirsche/jasmine_data/figures/figure5/survivorallelefreqs_long.png"
plot_allele_frequencies(variants, outfile)
outfile <- "/home/mkirsche/jasmine_data/figures/figure5/survivorallelefreqs_long.svg"
plot_allele_frequencies(variants, outfile) 

outfile <- "/home/mkirsche/jasmine_data/figures/figure5/svtoolsallelefreqs.png"
variants = read.table(svtoolsfile, comment.char = "#", sep = "\t", header = T, stringsAsFactors=FALSE, colClasses = 'character')
variants <- variants %>% filter(SPECIFIC_FLAG == "1" & PRECISE_FLAG == "1")
plot_allele_frequencies(variants, outfile)
outfile <- "/home/mkirsche/jasmine_data/figures/figure5/svtoolsallelefreqs.svg"
plot_allele_frequencies(variants, outfile)
variants <- variants%>% filter(abs(as.numeric(LEN)) >= 50 | TYPE == "TRA")
outfile <- "/home/mkirsche/jasmine_data/figures/figure5/svtoolsallelefreqs_long.png"
plot_allele_frequencies(variants, outfile)
outfile <- "/home/mkirsche/jasmine_data/figures/figure5/svtoolsallelefreqs_long.svg"
plot_allele_frequencies(variants, outfile) 

outfile <- "/home/mkirsche/jasmine_data/figures/figure5/svimmerallelefreqs.png"
variants = read.table(svimmerfile, comment.char = "#", sep = "\t", header = T, stringsAsFactors=FALSE, colClasses = 'character')
variants <- variants %>% filter(SPECIFIC_FLAG == "1" & PRECISE_FLAG == "1")
plot_allele_frequencies(variants, outfile)
outfile <- "/home/mkirsche/jasmine_data/figures/figure5/svimmerallelefreqs.svg"
plot_allele_frequencies(variants, outfile)
variants <- variants%>% filter(abs(as.numeric(LEN)) >= 50 | TYPE == "TRA")
outfile <- "/home/mkirsche/jasmine_data/figures/figure5/svimmerallelefreqs_long.png"
plot_allele_frequencies(variants, outfile)
outfile <- "/home/mkirsche/jasmine_data/figures/figure5/svimmerallelefreqs_long.svg"
plot_allele_frequencies(variants, outfile) 

outfile <- "/home/mkirsche/jasmine_data/figures/figure5/dbsvmergeallelefreqs.png"
variants = read.table(dbsvmergefile, comment.char = "#", sep = "\t", header = T, stringsAsFactors=FALSE, colClasses = 'character')
variants <- variants %>% filter(SPECIFIC_FLAG == "1" & PRECISE_FLAG == "1")
plot_allele_frequencies(variants, outfile)
outfile <- "/home/mkirsche/jasmine_data/figures/figure5/dbsvmergeallelefreqs.svg"
plot_allele_frequencies(variants, outfile)
variants <- variants%>% filter(abs(as.numeric(LEN)) >= 50 | TYPE == "TRA")
outfile <- "/home/mkirsche/jasmine_data/figures/figure5/dbsvmergeallelefreqs_long.png"
plot_allele_frequencies(variants, outfile)
outfile <- "/home/mkirsche/jasmine_data/figures/figure5/dbsvmergeallelefreqs_long.svg"
plot_allele_frequencies(variants, outfile) 

outfile <- "/home/mkirsche/jasmine_data/figures/figure5/svpopallelefreqs.png"
variants = read.table(svpopfile, comment.char = "#", sep = "\t", header = T, stringsAsFactors=FALSE, colClasses = 'character')
variants <- variants %>% filter(SPECIFIC_FLAG == "1" & PRECISE_FLAG == "1")
plot_allele_frequencies(variants, outfile)
outfile <- "/home/mkirsche/jasmine_data/figures/figure5/svpopallelefreqs.svg"
plot_allele_frequencies(variants, outfile)

get_pairwise <- function(fn) {
  variants = read.table(fn, comment.char = "#", sep = "\t", header = T, stringsAsFactors=FALSE, colClasses = 'character')
  variants <- variants %>% filter(SPECIFIC_FLAG == "1" & PRECISE_FLAG == "1")
  nrow(variants)
  length <- max(nchar(variants$SUPP_VEC))
  mat <- array(dim = c(length, length))
  for(i in 1:length) {
    for(j in 1:length) {
      #mat[i, j] = nrow(variants %>% filter(charAt(as.character(filter$SUPP_VEC, i)) == '1'))
      mat[i, j] = sum(substr(variants$SUPP_VEC,i,i) == "1" & substr(variants$SUPP_VEC,j,j) == "1")
      if(i == j) {
        #mat[i, j] = 0
      }
    } 
  }
  return(mat)
}
mat <- get_pairwise(jasminefile)

#mat_clr100 <- get_pairwise("/home/mkirsche/jasmine_data/figures/figure5/clr_100_specprec.txt")
mat
library(gplots)

# Read in sample info
fn <- "/home/mkirsche/jasmine_data/figures/figure5/filelist.txt"
samples <- read.table(fn, sep = "\t", header = T, colClasses = 'character')
samples$SAMPLE
samples$Label <- paste(samples$SAMPLE, " (", samples$TECH, ")", sep = "")

rownames(mat) <- samples$Label
colnames(mat) <- samples$Label
#rownames(mat_clr100) <- samples$Label
#colnames(mat_clr100) <- samples$Label

samples$COVERAGE_NUM <- as.numeric(samples$COVERAGE)
samples$NUMVARS <- as.numeric(samples$NUMVARS)

outfile <- "/home/mkirsche/jasmine_data/figures/figure5/numvars_singlesample.png"
ggplot(samples, aes(x = COVERAGE_NUM, y = NUMVARS, label = SAMPLE, color = TECH)) + geom_point() + xlab("Coverage") + ylab("Number of Variants") +
  theme(plot.title = element_text(size = 24, hjust = 0.5),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        #legend.text = element_blank(),
        text = element_text(size = 8),
        #legend.position = "none"
  ) + geom_text(nudge_y = -400, size = 4) + guides(color=guide_legend(title="Tech"))
ggsave(outfile, width = 12, height = 8)

outfile <- "/home/mkirsche/jasmine_data/figures/figure5/numvars_singlesample_superpop.png"
ggplot(samples, aes(x = COVERAGE_NUM, y = NUMVARS, label = SAMPLE, color = SUPERPOPULATION)) + geom_point() + xlab("Coverage") + ylab("Number of Variants") +
  theme(plot.title = element_text(size = 24, hjust = 0.5),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        #legend.text = element_blank(),
        text = element_text(size = 8),
        #legend.position = "none"
  ) + geom_text(nudge_y = -400, size = 4) + guides(color=guide_legend(title="Superpopulation"))
ggsave(outfile, width = 12, height = 8)

outfile <- "/home/mkirsche/jasmine_data/figures/figure5/correlation.png"
png(file=outfile, width = 1024, height = 1024)
heatmap.2(mat, col =cm.colors(256), 
          scale = "column", trace = "none")
dev.off()
library(pheatmap)
library(grid)
library(RColorBrewer)
mat
mat_colors <- list(Superpopulation = brewer.pal(6, "Set1"))
names(mat_colors$Superpopulation) <- unique(samples$SUPERPOPULATION)
mat_colors
mat_col <- data.frame(Superpopulation = samples$SUPERPOPULATION)
rownames(mat_col) <- rownames(mat)
mat_col
outfile <- "/home/mkirsche/jasmine_data/figures/figure5/correlation_corrplot.png"
png(file=outfile, width = 1024, height = 1024)
pheatmap(mat,show_rownames=TRUE,show_colnames=TRUE, annotation_colors = mat_colors, annotation_col = mat_col, annotation_row = mat_col, drop_levels = TRUE,
         scale = "none",clustering_method="ward.D2",
         clustering_distance_cols="euclidean", fontsize = 16)
dev.off()

outfile <- "/home/mkirsche/jasmine_data/figures/figure5/correlation_corrplot.svg"
#svg(file=outfile, width = 1024, height = 1024)
Cairo(1024, 1024, file = outfile, type = "svg", bg = "white", dpi = 80)
pheatmap(mat,show_rownames=TRUE,show_colnames=TRUE, annotation_colors = mat_colors, annotation_col = mat_col, annotation_row = mat_col, drop_levels = TRUE,
         scale = "none",clustering_method="ward.D2",
         clustering_distance_cols="euclidean", fontsize = 12)
dev.off()
#ggplot(data = melt(mat), aes(x=Var1, y=Var2, fill=value)) + ggcorrplot()
#ggsave(outfile, width = 10, height = 10)

# Make Jasmine PCA plots
tmp = read.table(jasminefile, comment.char = "#", sep = "\t", header = T, stringsAsFactors=FALSE, colClasses = 'character')
nrow(tmp)
head(tmp$LEN)
variants <- tmp %>% filter(SPECIFIC_FLAG == "1" & PRECISE_FLAG == "1") #%>% filter(abs(as.numeric(LEN)) >= 50 | TYPE == "TRA")
nrow(variants)

make_pca_scatter <- function(variants, pref, samples)
{

  variants$SEP <- sub("\\s+$", "", gsub('(.{1})', '\\1 ', variants$SUPP_VEC))
  samples$SAMPLE
  sepdf <- variants %>% separate(SEP, unlist(samples$SAMPLE), sep=" ")
  #colnames(sepdf)
  #nrow(sepdf)
  #length(unique(variants$X.scratch4.mschatz1.mkirsche.JasmineFigures.figures.figure5.NA19238vGRCh38_wm_50md_PBCCS_sniffles.s2l20.refined.nSVtypes.ism.vcf))
  #nrow(sepdf %>% filter(NA19238 == "1" & X.scratch4.mschatz1.mkirsche.JasmineFigures.figures.figure5.NA19238vGRCh38_wm_50md_PBCCS_sniffles.s2l20.refined.nSVtypes.ism.vcf == "."))
  
  #test <- (sepdf %>% filter(NA19238 == "1"))
  #nrow(test)
  #head(test %>% filter(as.numeric(LEN) < 20 & as.numeric(SUPP) > 20 & TYPE == "INS"))
  #bad <- (test %>% filter(as.numeric(LEN) >= 20 & as.numeric(LEN) < 25 & as.numeric(SUPP) > 30 & TYPE == "INS"))
  #ggplot(bad, aes(x = as.numeric(SPECIFIC_COUNT))) + geom_bar(stat = "count")
  #ggsave("/home/mkirsche/jasmine_data/figures/figure5/na19238_badspecificcount.png")
  
  #nrow(bad)
  #head(bad$X.scratch4.mschatz1.mkirsche.JasmineFigures.figures.figure5.NA19238vGRCh38_wm_50md_PBCCS_sniffles.s2l20.refined.nSVtypes.ism.vcf)
  #ggplot(test %>% filter(TYPE == "INS" & as.numeric(LEN) < 50), aes(x = as.numeric(LEN))) + geom_bar(stat = "count")
  #ggsave("/home/mkirsche/jasmine_data/figures/figure5/na19238_smallinslen.png")
  #ggplot(test, aes(x = as.numeric(SUPP))) + geom_bar(stat = "count")
  #ggsave("/home/mkirsche/jasmine_data/figures/figure5/na19238_supp.png")
  #ggplot(test%>% filter(abs(as.numeric(LEN)) >= 50 | TYPE == "TRA"), aes(x = as.numeric(SUPP))) + geom_bar(stat = "count")
  #ggsave("/home/mkirsche/jasmine_data/figures/figure5/na19238_supp_big.png")
  #ins <- test %>% filter(TYPE == "INS")%>% filter(abs(as.numeric(LEN)) < 1000000)%>% filter(abs(as.numeric(LEN)) < 50)
  #del <- test %>% filter(TYPE == "DEL") %>% filter(abs(as.numeric(LEN)) < 1000000)%>% filter(abs(as.numeric(LEN)) < 50)
  #sum(as.numeric(ins$LEN))
  #sum(as.numeric(del$LEN))
  #plot_length2(test, "/home/mkirsche/jasmine_data/figures/figure5/na19238.png", .88, .88)

  #test <- (sepdf %>% filter(HG006 == "1"))
  #plot_length2(test, "/home/mkirsche/jasmine_data/figures/figure5/hg006.png", .88, .88)

  #test <- (sepdf %>% filter(HG02106 == "1"))
  #plot_length2(test, "/home/mkirsche/jasmine_data/figures/figure5/hg02106.png", .88, .88)
  
  #test <- (sepdf %>% filter(HG00732 == "1"))
  #plot_length2(test, "/home/mkirsche/jasmine_data/figures/figure5/hg00732.png", .88, .88)
  
  #test <- (sepdf %>% filter(HG02723 == "1"))
  #plot_length2(test, "/home/mkirsche/jasmine_data/figures/figure5/hg02723.png", .88, .88)

  
  pcadf <- sepdf[,unlist(samples$SAMPLE)]
  colnames(pcadf)
  pcamatrix<-t(data.matrix(pcadf))
  pca <- prcomp(pcamatrix, center = TRUE)
  df_out <- as.data.frame(pca$x)
  df_out$SAMPLE <- rownames(df_out)
  df_out <- merge(df_out, samples, by = "SAMPLE")
  ggplot(df_out, aes(x = PC1, y = PC2, label = SAMPLE, color = POPULATION)) + geom_point() + geom_text(nudge_y = -1.5, size=2)
  outfile <- paste("/home/mkirsche/jasmine_data/figures/figure5/", pref, "_", "pcapop.png", sep = "")
  ggsave(outfile, width = 8, height = 8)
  ggplot(df_out, aes(x = PC1, y = PC2, label = SAMPLE, color = SUPERPOPULATION)) + geom_point() + geom_text(nudge_y = -1.5, size=2)
  outfile <- paste("/home/mkirsche/jasmine_data/figures/figure5/", pref, "_", "pcasuperpop.png", sep = "")
  ggsave(outfile, width = 8, height = 8)
  
  # Make svpop PCA plots
  #variants = read.table(svpopfile, comment.char = "#", sep = "\t", header = T, stringsAsFactors=FALSE, colClasses = 'character')
  #variants <- variants %>% filter(SPECIFIC_FLAG == "1" & PRECISE_FLAG == "1")
  #variants$SEP <- sub("\\s+$", "", gsub('(.{1})', '\\1 ', variants$SUPP_VEC))
  #samples$SAMPLE
  #sepdf <- variants %>% separate(SEP, unlist(samples$SAMPLE), sep=" ")
  #colnames(sepdf)
  #pcadf <- sepdf[,unlist(samples$SAMPLE)]
  #colnames(pcadf)
  #pcamatrix<-t(data.matrix(pcadf))
  #pca <- prcomp(pcamatrix, center = TRUE)
  #df_out <- as.data.frame(pca$x)
  #df_out$SAMPLE <- rownames(df_out)
  #df_out <- merge(df_out, samples, by = "SAMPLE")
  #ggplot(df_out, aes(x = PC1, y = PC2, label = SAMPLE, color = POPULATION)) + geom_point() + geom_text(nudge_y = -1.5, size=2)
  #outfile <- "/home/mkirsche/jasmine_data/figures/figure5/svpoppcapop.png"
  #ggsave(outfile, width = 8, height = 8)
  #ggplot(df_out, aes(x = PC1, y = PC2, label = SAMPLE, color = SUPERPOPULATION)) + geom_point() + geom_text(nudge_y = -1.5, size=2)
  #outfile <- "/home/mkirsche/jasmine_data/figures/figure5/svpoppcasuperpop.png"
  #ggsave(outfile, width = 8, height = 8)
  
  # Get the number of SV calls in each sample
  sampletally <- sepdf %>% gather(x, value, unlist(samples$SAMPLE)) %>% group_by(x) %>% tally(value==1)
  names(sampletally)[names(sampletally) == 'x'] <- 'SAMPLE'
  colnames(sampletally)
  colnames(samples)
  sampletally <- merge(sampletally, samples, by = "SAMPLE")
  sampletally$TECH
  ggplot(sampletally, aes(x = SAMPLE, y = n, fill = TECH)) + geom_bar(stat = "identity") + xlab("Sample") + ylab("Number of Variants") +
    theme(plot.title = element_text(size = 24, hjust = 0.5),
          axis.text.x = element_text(size = 10, angle = 40),
          axis.text.y = element_text(size = 10),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          #legend.text = element_blank(),
          text = element_text(size = 16),
          #legend.position = "none"
    ) + guides(fill=guide_legend(title="Tech"))
  outfile <- paste("/home/mkirsche/jasmine_data/figures/figure5/", pref, "_", "samplecountsbar.png", sep = "")
  ggsave(outfile, width = 12, height = 8)
  sampletally
  sampletally$numcov <- as.numeric(sampletally$COVERAGE)
  ggplot(sampletally, aes(x = numcov, y = n, label = SAMPLE, color = TECH)) + geom_point() + xlab("Coverage") + ylab("Number of Variants") +
    theme(plot.title = element_text(size = 24, hjust = 0.5),
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          #legend.text = element_blank(),
          text = element_text(size = 8),
          #legend.position = "none"
    ) + geom_text(nudge_y = -400, size = 4) + guides(color=guide_legend(title="Tech"))
  outfile <- paste("/home/mkirsche/jasmine_data/figures/figure5/", pref, "_", "samplecountsscatter.png", sep = "")
  ggsave(outfile, width = 12, height = 8)
  
  ggplot(sampletally, aes(x = numcov, y = n, label = SAMPLE, color = TECH)) + geom_point() + xlab("Coverage") + ylab("Number of Variants") +
    theme(plot.title = element_text(size = 24, hjust = 0.5),
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          #legend.text = element_blank(),
          text = element_text(size = 8),
          #legend.position = "none"
    ) + geom_text(nudge_y = -400, size = 4) + guides(color=guide_legend(title="Tech"))
  outfile <- paste("/home/mkirsche/jasmine_data/figures/figure5/", pref, "_", "samplecountsscatter.svg", sep = "")
  ggsave(outfile, width = 12, height = 8)
}
make_pca_scatter(variants, "big", samples)
make_pca_scatter(tmp, "allvars", samples)

#violin_plot(jasminefile, "jasmine")
#violin_plot(svpopfile, "svpop")

plot_length2 <- function(df, outfile, legendx, legendy) 
{
  df$LenCategory = "TRA"
  
  df$ABSLEN = abs(strtoi(df$LEN))
  
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
  
  colorpalette <- c(DEL = '#CC6677',DUP = '#332288',INS = '#44AA99',TRA = '#88CCEE',INV = "#ddaa33")#'#DDCC77')
  
  summarized <- df %>% group_by(TYPE, LenCategory) %>% summarise(counts=n())
  
  ggplot(summarized) +
    geom_bar(data = summarized %>% filter(TYPE != "DEL"), position = "stack", stat = "identity", aes(x = LenCategory, fill = TYPE, y = counts))+
    geom_bar(data = summarized %>% filter(TYPE == "DEL"), position = "stack", stat = "identity", aes(x = LenCategory, fill = TYPE, y = -counts))+
    scale_x_discrete(labels=labellist) +
    xlab("Length") +
    ylab("Count") +
    labs(title = "Population SV Size Distribution") +
    theme(plot.title = element_text(size = 24, hjust = 0.5),
          axis.text.x = element_text(size = 18, angle = 30, margin = margin(t = 17, r = 0, b = 0, l = 0)),
          axis.text.y = element_text(size = 18),
          axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20),
          legend.text = element_text(size = 18),
          legend.title = element_text(size = 20),
          text = element_text(size = 20),
          legend.position = c(legendx, legendy),
    ) +
    scale_fill_manual(name = "TYPE", values = colorpalette) +
    geom_text(data = nondelcounts, aes(x = LenCategory, y=n, label=n), position=position_dodge(width=0.9), vjust=-0.75) +
    geom_text(data = delcounts, aes(x = LenCategory, y=-1*n, label=n), position=position_dodge(width=0.9), vjust=1.5)
  
  ggsave(outfile, width= 12, height = 8)
}
