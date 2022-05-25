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
library(svglite)

fn <- '/home/mkirsche/jasmine_data/figures/supplement/jasmine_dist/jasmine_md_clr50.merged.tsv'
df <- read.table(fn, sep = "\t", header = TRUE)
filtered = df %>% filter(IS_SPECIFIC == 1 & IS_PRECISE == 1)
nrow(filtered)
nrow(filtered %>% filter(SUPP_VEC == "100"))
filtered <- filtered %>% filter(SVLEN == 0 | SVLEN >= 50 | SVLEN <= -50 | SVTYPE == "TRA")
nrow(filtered)
nrow(filtered %>% filter(SUPP_VEC == "100"))

fn <- '/home/mkirsche/jasmine_data/figures/supplement/jasmine_dist/jasmine_md_ont50.merged.tsv'
df <- read.table(fn, sep = "\t", header = TRUE)
filtered = df %>% filter(IS_SPECIFIC == 1 & IS_PRECISE == 1)
nrow(filtered)
nrow(filtered %>% filter(SUPP_VEC == "100"))
filtered <- filtered %>% filter(SVLEN == 0 | abs(SVLEN) >= 50 | SVTYPE == "TRA")
nrow(filtered)
unique(filtered$SUPP_VEC)
filtered %>% count(SUPP_VEC)
nrow(filtered %>% filter(SUPP_VEC == "100"))

fn <- '/home/mkirsche/jasmine_data/figures/figure2/jasmine_md50.merged.tsv'
df <- read.table(fn, sep = "\t", header = TRUE)
filtered = df %>% filter(IS_SPECIFIC == 1 & IS_PRECISE == 1)
nrow(filtered)
nrow(filtered %>% filter(SUPP_VEC == "100"))
filtered <- filtered %>% filter(SVLEN == 0 | SVLEN >= 50 | SVLEN <= -50 | SVTYPE == "TRA")
nrow(filtered)
unique(filtered$SUPP_VEC)
filtered %>% count(SUPP_VEC)
nrow(filtered %>% filter(SUPP_VEC == "100"))

plot_table <- function(fn, outdir, prefix, caller, discvec, lengthfilter){
  
  df <- read.table(fn, sep = "\t", header = TRUE)
  filtered = df %>% filter(SPECIFIC_FLAG == 1 & PRECISE_FLAG == 1)
  if (lengthfilter) {
    filtered <- filtered %>% filter(LEN == 0 | LEN >= 50 | LEN <= -50 | TYPE == "TRA")
    prefix <- paste(prefix, "_50plus", sep = "")
  }
  filtered$diff <- filtered$NUMVARS - filtered$SUPP
  diffcount <- filtered %>% group_by(diff) %>% summarise(counts=log2(n()))
  filtered <- arrange(filtered, SUPP)
  ggplot(filtered, aes(x = diff, y = 1, fill = as.character(SUPP))) + geom_bar(position = "stack", stat = "identity") +
    scale_fill_discrete(labels=c("1", "2", "3"), name = "Number of Samples") +
    theme_bw() +
    labs(title = paste("Intrasample Merging (", caller, ")", sep = "")) +
    theme(axis.text.x = element_text(size = 16),
          axis.text.y = element_text(size = 16),
          axis.title.x = element_text(size = 18),
          axis.title.y = element_text(size = 18),
          plot.title = element_text(hjust = 0.5, size = 20)) +
    xlab("Number of excess variants merged") +
    ylab("Count")
  #ggsave(paste(outdir, "/", prefix, "_", caller, ".numvars.png", sep = ""), width = 8, height = 8)
  ggplot(diffcount, aes(x = diff, y = counts, fill = diff)) + geom_bar(stat = "identity") +
    theme_bw() +
    scale_fill_distiller(palette = "Set2") +
    scale_x_continuous(limits = c(-5, 150)) +
    labs(title = paste("Intrasample Merging (", caller, ")", sep = "")) +
    theme(axis.text.x = element_text(size = 16),
          axis.text.y = element_text(size = 16),
          axis.title.x = element_text(size = 18),
          axis.title.y = element_text(size = 18),
          legend.position = "none",
          plot.title = element_text(hjust = 0.5, size = 20)) +
    xlab("Number of excess variants merged") +
    ylab("log2(Count)")
  #ggsave(paste(outdir, "/", prefix, "_", caller, ".log.numvars.png", sep = ""), width = 8, height = 8)
  filtered = filtered %>% filter(SUPP_VEC == discvec)
  filtered$caller <- caller
  filtered$Intrasample = 'Allows Intrasample Merging'
  filtered$Intrasample <- ifelse(filtered$caller == "Jasmine", "No Intrasample Merging", as.character(filtered$Intrasample))
  filtered$Intrasample <- ifelse(filtered$caller == "dbsvmerge", "No Intrasample Merging", as.character(filtered$Intrasample))
  filtered$Intrasample <- ifelse(filtered$caller == "svpop", "No Intrasample Merging", as.character(filtered$Intrasample))
  filtered$Intrasample <- ifelse(filtered$caller == "svmerger", "No Intrasample Merging", as.character(filtered$Intrasample))
  
  
  specprec <- df %>% filter(SPECIFIC_FLAG == 1 & PRECISE_FLAG == 1)
  if (lengthfilter) {
    specprec <- specprec %>% filter(LEN == 0 | LEN >= 50 | LEN <= -50 | TYPE == "TRA")
  }
  totalcount <- sum(specprec$NUMVARS)
  mergedcount <- nrow(specprec)
  return (list(filtered, totalcount, mergedcount))
}

plot_all_callers <- function(outdir, prefix, discvec, legendx, legendy, width, height, lengthfilter){
  jasmineresults <- plot_table(paste(outdir, prefix, ".jasmine_augmented.txt", sep = ""), outdir, prefix, "Jasmine", discvec, lengthfilter)
  jasmineresults[[3]]
  survivorresults <- plot_table(paste(outdir, prefix, ".survivor_augmented.txt", sep = ""), outdir, prefix, "Survivor", discvec, lengthfilter)
  svtoolsresults <- plot_table(paste(outdir, prefix, ".svtools_augmented.txt", sep = ""), outdir, prefix, "svtools", discvec, lengthfilter)
  svimmerresults <- plot_table(paste(outdir, prefix, ".svimmer_augmented.txt", sep = ""), outdir, prefix, "svimmer", discvec, lengthfilter)
  jasmineintraresults <- plot_table(paste(outdir, prefix, ".jasmineintra_augmented.txt", sep = ""), outdir, prefix, "Jasmine_Intrasample", discvec, lengthfilter)
  svpopresults <- plot_table(paste(outdir, prefix, ".svpop_augmented.txt", sep = ""), outdir, prefix, "svpop", discvec, lengthfilter)
  dbsvmergeresults <- plot_table(paste(outdir, prefix, ".dbsvmerge_augmented.txt", sep = ""), outdir, prefix, "dbsvmerge", discvec, lengthfilter)
  svmergerresults <- plot_table(paste(outdir, prefix, ".svmerger_augmented.txt", sep = ""), outdir, prefix, "svmerger", discvec, lengthfilter)
 
  ratetitle <- "Discordant and Invalid Variants in Child"
  fnprefix <- prefix
  if (lengthfilter) {
    fnprefix <- paste(fnprefix, "_50plus", sep = "")
    legendx = .67
    ratetitle <- "Discordant and Invalid SVs in Child"
  }
  discordant <- data.frame()
  discordant <- rbind(discordant, jasmineresults[[1]])
  colnames(discordant)
  colnames(survivorresults[[1]])
  jasmineresults[[1]]
  discordant <- rbind(discordant, survivorresults[[1]])
  discordant <- rbind(discordant, svtoolsresults[[1]])
  discordant <- rbind(discordant, svimmerresults[[1]])
  discordant <- rbind(discordant, jasmineintraresults[[1]])
  discordant <- rbind(discordant, svpopresults[[1]])
  discordant <- rbind(discordant, dbsvmergeresults[[1]])
  colnames(svmergerresults[[1]])[2] = colnames(discordant)[2]
  colnames(svmergerresults[[1]])[3] = colnames(discordant)[3]
  colnames(svmergerresults[[1]])[4] = colnames(discordant)[4]
  discordant <- rbind(discordant, svmergerresults[[1]])
  nrow(discordant %>% filter(caller == "Jasmine"))
  nrow(jasmineresults[[1]])
  
  errors <- read.table(paste(outdir, prefix, ".errors.txt", sep = ""), sep="\t", header = TRUE)
  errors$caller <- errors$SOFTWARE
  errors$caller <- ifelse(errors$caller == "survivor", "Survivor", as.character(errors$caller))
  errors$caller <- ifelse(errors$caller == "jasmine", "Jasmine", as.character(errors$caller))
  errors$caller <- ifelse(errors$caller == "jasmineintra", "Jasmine_Intrasample", as.character(errors$caller))
  errors$totaldisc <- errors$DISC_MIXED_STRAND_AND_TYPE + errors$DISC_MIXED_STRAND_ONLY + errors$DISC_MIXED_TYPE_ONLY
  errors$mixedtype <- errors$DISC_MIXED_STRAND_AND_TYPE + errors$DISC_MIXED_TYPE_ONLY
  errors$mixedstrand <- errors$DISC_MIXED_STRAND_ONLY
  errors$ism <- errors$DISC_COMBINED_ISM
  errors
  
  #totals <- c(jasmineresults[[2]], survivorresults[[2]], svtoolsresults[[2]], svimmerresults[[2]], jasmineintraresults[[2]],
  #            jasmineresults[[3]], survivorresults[[3]], svtoolsresults[[3]], svimmerresults[[3]], jasmineintraresults[[3]])
  discordant
  disccounts <- discordant %>% group_by(caller, Intrasample) %>% summarise(counts=n(), sums=sum(NUMVARS))
  disccounts
  jasmineresults[[2]]
  jasmineresults[[3]]
  totaldf <- data.frame(matrix(ncol = 2, nrow = 8))
  colnames(totaldf) <- c("caller", "total_vars")
  totallist <- c(jasmineresults[[3]], survivorresults[[3]], svtoolsresults[[3]], svimmerresults[[3]], 
                 jasmineintraresults[[3]], svpopresults[[3]], dbsvmergeresults[[3]], svmergerresults[[3]])
  callerlist <- c("Jasmine", "Survivor", "svtools", "svimmer", "Jasmine_Intrasample", "svpop", "dbsvmerge", "svmerger")
  totaldf$caller <- callerlist
  totaldf$total_vars <- totallist
  totaldf
  
  disccounts <- inner_join(disccounts, errors, by = c("caller"))
  disccounts <- inner_join(disccounts, totaldf, by = c("caller"))
  disccounts
  colnames(disccounts)
  disccounts$total_vars <- disccounts$total_vars + disccounts$ism
  mdfr <- melt(disccounts, id.vars = c("caller", "Intrasample", "total_vars"))
  mdfr
  ggplot(disccounts, aes(x = caller, y = counts, fill = caller)) +
    geom_bar(position = "stack", stat = "identity") +
    labs(title = paste("Discordant Variants (", prefix, ")")) +
    xlab("Merging Software") +
    ylab("Discordant Count") +
    theme(plot.title = element_text(size = 22, hjust = 0.5),
          axis.text.x = element_text(size = 16),
          axis.text.y = element_text(size = 16),
          axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20),
    ) + theme(legend.position = "none") + facet_grid(cols = vars(Intrasample), scales = "free")
  outfile <- paste(outdir, fnprefix, "_discordant.png", sep="")
  ggsave(outfile, width= 8, height = 8)
  
  mdfr <- mdfr %>% filter(caller != "Jasmine_Intrasample")
  mdfr

  ggplot(mdfr %>% filter(variable == "counts" | variable == "mixedtype" | variable == "mixedstrand" | variable == "ism"), 
         aes(x = factor(caller, levels = c("Jasmine", "svpop", "dbsvmerge", "Survivor", "svimmer", "svmerger", "svtools")), y = as.numeric(value), fill = variable)) +
    geom_bar(stat = "identity") +
    labs(title = paste(ratetitle, " (", prefix, ")")) +
    xlab("Merging Software") +
    ylab("Count") +
    theme(plot.title = element_text(size = 18, hjust = 0.5),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          legend.position = c(legendx, legendy),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 14),
    ) +
    scale_fill_brewer(name = "Type", labels = c("Discordant", "Mixed Type", "Mixed Strand", "Discordant (Intrasample)"), palette = "Set2") + 
    facet_grid(cols = vars(factor(Intrasample,labels = c("No Intrasample Merging", "Allows Intrasample Merging"))), scales = "free")
  outfile <- paste(outdir, fnprefix, "_discordant_errors.png", sep="")
  ggsave(outfile, width= 10, height = 8)
  
  colorpalette <- c(counts = "#882255", mixedtype = "#AA4499", mixedstrand = "#332288", ism = "#88CCEE")
  ggplot(mdfr %>% filter(variable == "counts" | variable == "mixedtype" | variable == "mixedstrand" | variable == "ism"), 
         aes(x = factor(caller, levels = c("Jasmine", "svpop", "dbsvmerge", "svmerger", "svtools", "svimmer", "Survivor")), y = as.numeric(value) / as.numeric(total_vars), fill = variable)) +
    geom_bar(stat = "identity") +
    labs(title = paste(ratetitle, " (", "HG002 HiFi", ")", sep = "")) +
    xlab("Merging Software") +
    ylab("Proportion") +
    theme(plot.title = element_text(size = 24, hjust = 0.5),
          axis.text.x = element_text(size = 22,angle = 15),
          axis.text.y = element_text(size = 22),
          axis.title.x = element_text(size = 22),
          axis.title.y = element_text(size = 22),
          legend.position = c(legendx, legendy),
          legend.text = element_text(size = 18),
          legend.title = element_blank(),
          strip.text.x = element_text(size = 24),
          strip.text.y = element_text(size = 24)
    ) +
    scale_fill_manual(name = "Type", labels = c("Discordant", "Mixed Type", "Mixed Strand", "Discordant (Intrasample)"), values = colorpalette) + 
    facet_grid(cols = vars(factor(Intrasample,levels = c("No Intrasample Merging", "Allows Intrasample Merging"))), scales = "free")
  outfile <- paste(outdir, fnprefix, "_discordant_error_rate.png", sep="")
  ggsave(outfile, width= width, height = height)
  outfile <- paste(outdir, fnprefix, "_discordant_error_rate.svg", sep="")
  ggsave(outfile, width= width, height = height)
}

plot_length <- function(df, caller, outfile, filter, legendx, legendy) 
{
  if (filter)
  {
    df <- df %>% filter(IS_SPECIFIC == 1 & IS_PRECISE == 1)
  }
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
  
  colorpalette <- c(DEL = '#CC6677',DUP = '#332288',INS = '#44AA99',TRA = '#88CCEE',INV = '#DDCC77')
  
  summarized <- df %>% group_by(SVTYPE, LenCategory) %>% summarise(counts=n())
  
  ggplot(summarized, aes(fill = SVTYPE, x = LenCategory)) +
    geom_bar(data = summarized %>% filter(SVTYPE != "DEL"), stat = "identity", aes(y = counts), position = "stack")+
    geom_bar(data = summarized %>% filter(SVTYPE == "DEL"), stat = "identity", aes(y = -counts), position = "stack")+
    scale_x_discrete(labels=labellist) +
    xlab("Length") +
    ylab("Count") +
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
    ggtitle('HG002 Trio Variant Size Distribution (Jasmine)') +
    #scale_fill_brewer(name = "Type", palette = "Set2") + 
    scale_fill_manual(name = "SVTYPE", values = colorpalette) +
    guides(fill=guide_legend(title="Type")) +
    geom_text(data = nondelcounts, aes(x = LenCategory, y=n, label=n), position=position_dodge(width=0.9), vjust=-0.75) +
    geom_text(data = delcounts, aes(x = LenCategory, y=-1*n, label=n), position=position_dodge(width=0.9), vjust=1.5)
  
  ggsave(outfile, width= 10, height = 8)
}

plot_length_line <- function(df, caller, outfile, filter, legendx, legendy) 
{
  if (filter)
  {
    df <- df %>% filter(IS_SPECIFIC == 1 & IS_PRECISE == 1)
  }
  df$ABSLEN = abs(as.numeric(df$SVLEN))
  
  colorpalette <- c(DEL = '#CC6677',DUP = '#332288',INS = '#44AA99',TRA = '#88CCEE',INV = '#DDCC77')
  
  #summarized <- df %>% group_by(SVTYPE) %>% summarise(counts=n())
  
  ggplot(df%>% filter(ABSLEN <=10000 & ABSLEN >= 100), aes(x = ABSLEN, color = SVTYPE)) +
    geom_density()+
    xlab("Length") +
    ylab("Density") +
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
    ggtitle('HG002 Trio Variant Size Distribution (Jasmine)') +
    #scale_fill_brewer(name = "Type", palette = "Set2") + 
    scale_color_manual(name = "SVTYPE", values = colorpalette) +
    guides(color=guide_legend(title="Type"))
    #geom_text(data = nondelcounts, aes(x = LenCategory, y=n, label=n), position=position_dodge(width=0.9), vjust=-0.75) +
    #geom_text(data = delcounts, aes(x = LenCategory, y=-1*n, label=n), position=position_dodge(width=0.9), vjust=1.5)
  
  ggsave(outfile, width= 10, height = 8)
}

suppvec_hist <- function(df, caller, outfile, filter, lengthfilter) {
  if (filter)
  {
    df <- df %>% filter(IS_SPECIFIC == 1 & IS_PRECISE == 1)
  }
  
  if (lengthfilter) {
    df <- df %>% filter(SVLEN == 0 | abs(SVLEN) >= 50 | SVTYPE == "TRA")
  }
  
  df$SUPP_VEC_STRING <- str_pad(as.character(df$SUPP_VEC), 3, "left", "0")
  df$SUPP_VEC_STRING <- ifelse(df$SUPP_VEC_STRING == "100", "HG002 Only", df$SUPP_VEC_STRING)
  df$SUPP_VEC_STRING <- ifelse(df$SUPP_VEC_STRING == "010", "HG003 Only", df$SUPP_VEC_STRING)
  df$SUPP_VEC_STRING <- ifelse(df$SUPP_VEC_STRING == "001", "HG004 Only", df$SUPP_VEC_STRING)
  df$SUPP_VEC_STRING <- ifelse(df$SUPP_VEC_STRING == "110", "HG002/HG003", df$SUPP_VEC_STRING)
  df$SUPP_VEC_STRING <- ifelse(df$SUPP_VEC_STRING == "101", "HG002/HG004", df$SUPP_VEC_STRING)
  df$SUPP_VEC_STRING <- ifelse(df$SUPP_VEC_STRING == "011", "HG003/HG004", df$SUPP_VEC_STRING)
  df$SUPP_VEC_STRING <- ifelse(df$SUPP_VEC_STRING == "111", "All three", df$SUPP_VEC_STRING)
  
  suppveccounts <- df %>% count(SUPP_VEC_STRING)
  suppveccounts
  suppveccounts$SVTYPE = "INS"
  ggplot(df, aes(x = SUPP_VEC_STRING, y = 1, fill = SVTYPE)) +
    geom_bar(position = "stack", stat = "identity") +
    labs(title = paste("Variants by Sample Presence (", caller, ")", sep = "")) +
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
    guides(fill=guide_legend(title="Type")) +
    geom_text(data = suppveccounts, aes(x = SUPP_VEC_STRING, y=n, label=n), position=position_dodge(width=0.9), vjust=-0.75)
  
  ggsave(outfile, width= 7, height = 8)
}

suppvec_hist_highlightdiscordant <- function(df, caller, outfile, filter, offset, legendx, legendy, title, titlesize, lengthfilter) {
  if (filter)
  {
    df <- df %>% filter(IS_SPECIFIC == 1 & IS_PRECISE == 1)
  }
  if (lengthfilter) {
    df <- df %>% filter(SVLEN == 0 | abs(SVLEN) >= 50 | SVTYPE == "TRA")
  }
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
  
  suppveccounts <- df %>% count(SUPP_VEC_STRING)
  suppveccounts
  suppveccounts$SVTYPE = "INS"
  
  colorpalette <- c(DEL = '#CC6677',DUP = '#332288',INS = '#44AA99',TRA = '#88CCEE',INV = '#DDCC77')
  
  childcounts <- df %>% filter(df$SUPP_VEC_STRING == "HG002 Only") %>% group_by(SUPP_VEC_STRING) %>% count()
  childcounts$SVTYPE = "INS"
  childcounts$shape = 23
  ggplot(df, aes(x = SUPP_VEC_STRING, fill = SVTYPE)) +
    geom_bar(position = "stack", stat = "count") +
    xlab("Samples") +
    ylab("Count") +
    scale_y_continuous(expand = expansion(mult = c(0, .1))) +
    theme(axis.text.x = element_text(size = 14, angle = 30, margin = margin(t = 21)),
          axis.text.y = element_text(size = 16),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 20),
          legend.text = element_text(size = 16),
          legend.title = element_text(size = 20),
          legend.position = c(legendx, legendy),
          plot.title = element_text(size = titlesize, hjust = 0.5),
    ) +
    ggtitle(paste("Mendelian Discordance (", title, ")", sep = "")) +
    scale_fill_manual(name = "SVTYPE", values = colorpalette) +
    guides(fill=guide_legend(title="Type")) +
    geom_text(data = suppveccounts, aes(x = SUPP_VEC_STRING, y=n, label=n), position=position_dodge(width=0.9), vjust=-0.75) +
    geom_point(data = childcounts, aes(x = SUPP_VEC_STRING, y=n+offset), shape = 25, fill = "darkred", color = "darkred", size = 5)
  ggsave(outfile, width= 7, height = 8)
}

# Plot discordance across different Sniffles max_dist parameters
plotdisc <- function(outdir, indir, wildcard) {
  
  md_files <- Sys.glob(paste(indir, wildcard, sep = ""))
  
  names = c("File", "Total", "Discordant")
  md_data <- data.frame(setNames(rep(list(NA), length(names)), names))
  md_data$File = c()
  md_data$Total = c()
  md_data$Discordant = c()
  
  for(md_file in md_files) {
    md_file
    md_df = read.table(md_file, comment.char = "#", sep = "\t", header = T)
    md_df <- md_df %>% filter(CHROM != "MT" & SVTYPE != "=SR" 
                              & SUPP_VEC != "" & IS_SPECIFIC == "1" & IS_PRECISE == "1")

    nrow(md_df)
    
    md_df_filter <- md_df %>% filter(SUPP_VEC == "100")
    
    disccount <- nrow(md_df_filter)
    disccount
    totalcount <- nrow(md_df)
    
    cur <- data.frame(md_file, totalcount, disccount)
    cur
    names(cur) <- c("File", "Total", "Discordant")
    md_data <- rbind(md_data, cur)
  }
  
  md_data
  md_data$md = as.numeric(gsub(".merged.tsv", '', substring(md_data$File, str_length(indir) + 11)))
  #md_data$md = as.numeric(gsub(".merged.vcf", '', substring(md_data$File, 37)))
  md_data$md
  
  md_data$DiscordanceRate = md_data$Discordant / md_data$Total
  
  filterpref = ""
  
  discordantplot <- ggplot(data = md_data, aes(x = md, y = Discordant)) +
    geom_point(size = 1.6) +
    ggtitle(paste("Discordant Variant Calls", sep = "")) +
    ylab("Discordant Calls") +
    theme(axis.title.x = element_blank(), axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), plot.title = element_text(size = 16), axis.title.y = element_blank())
  
  outfile <- paste(outdir, filterpref, "mddisc.png",sep="")
  ggsave(outfile, device = "png", height = 8, width = 8)
  
  totalplot <- ggplot(data = md_data, aes(x = md, y = Total)) +
    geom_point(size = 1.6) +
    ggtitle(paste("Total Variant Calls", sep = "")) +
    ylab("Total Calls") +
    theme(axis.title.x = element_blank(), axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), axis.title.y = element_blank(), plot.title = element_text(size = 16))
  
  outfile <- paste(outdir, filterpref, "mdtotal.png",sep="")
  ggsave(outfile, device = "png", height = 8, width = 8)
  
  discrateplot <- ggplot(data = md_data, aes(x = md, y = DiscordanceRate)) +
    geom_point(size = 1.6) +
    ggtitle(paste("Discordance Rate", sep = "")) +
    xlab("Sniffles max_dist parameter") + theme(axis.title.x = element_text(size = 14), axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), axis.title.y = element_blank(), plot.title = element_text(size = 16))
  
  outfile <- paste(outdir, filterpref, "mddiscrate.png",sep="")
  ggsave(outfile, device = "png", height = 8, width = 8)
  ggarrange(arrangeGrob(discordantplot, totalplot, discrateplot, ncol = 1)) %>% ggexport(filename = paste(outdir, filterpref, "mdfullplot.png",sep = ""))
  ggarrange(arrangeGrob(discordantplot, totalplot, discrateplot, ncol = 1))
  ggsave(paste(outdir, filterpref, "mdfullplot.svg",sep = ""), width = 6, height = 6)
}

plotspec <- function(infile, outdir, legendx, legendy) {
  
  colorpalette <- c("#DC050C", "#4EB265", "#5289C7", "#882E72")
  
  df <- read.table(infile, sep = "\t", header = TRUE)
  df
  df$LENGTH
  df <- df %>% filter(LENGTH < 100 & READ_SUPPORT < 50)
  #df <- df[sample(nrow(df), 2000), ]
  df$Discordance <- ifelse(df$DISCORDANT == 0, "Not discordant", "Discordant")
  df$RESCUEDFROMABSENCE
  df$Discordance <- ifelse(df$RESCUEDFROMABSENCE == 1, "Rescued from absence", df$Discordance)
  df$Discordance <- ifelse(df$RESCUEDFROMDISCORDANCE == 1, "Rescued from discordance", df$Discordance)
  ggplot(df, mapping = aes(x = READ_SUPPORT, fill = Discordance)) + 
    scale_fill_manual(values = colorpalette) +
    xlab("Read Support") +
    ylab ("Count") +
    theme(plot.title = element_text(size = 22, hjust = 0.5),
          axis.text.x = element_text(size = 16),
          axis.text.y = element_text(size = 16),
          axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20),
          legend.text = element_text(size = 16),
          legend.title = element_text(size = 20),
          legend.position = c(legendx, legendy),
    ) +
    geom_bar(position = "stack", stat = "count") +
    scale_shape_identity()
  outfile <- paste(outdir, "/", "specificity_readsupp.png", sep = "")
  ggsave(outfile, width= 8, height = 8)
  outfile <- paste(outdir, "/", "specificity_readsupp.svg", sep = "")
  ggsave(outfile, width= 8, height = 8)
  
  ggplot(df, mapping = aes(x = LENGTH, fill = Discordance)) + 
    scale_fill_manual(values = colorpalette) +
    xlab("Length") +
    ylab ("Count") +
    theme(plot.title = element_text(size = 22, hjust = 0.5),
          axis.text.x = element_text(size = 16),
          axis.text.y = element_text(size = 16),
          axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20),
          legend.text = element_text(size = 16),
          legend.title = element_text(size = 20),
          legend.position = c(legendx, legendy),
    ) +
    geom_bar(position = "stack", stat = "count") +
    scale_shape_identity()
  outfile <- paste(outdir, "/", "specificity_length.png", sep = "")
  ggsave(outfile, width= 8, height = 8)
  outfile <- paste(outdir, "/", "specificity_length.svg", sep = "")
  ggsave(outfile, width= 8, height = 8)
  
  ggplot(df %>% filter(df$RESCUEDFROMDISCORDANCE == 1), mapping = aes(x=LENGTH, y= READ_SUPPORT)) +
    scale_fill_gradient(low="lightblue1", high = "darkblue",trans = "log2") +
    geom_bin2d() +
    ggtitle("Formerly Discordant Variants in HG002") +
    xlab("Variant Length") +
    ylab("Read Support") +
    theme(plot.title = element_text(size = 18, hjust = 0.5),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16),
    )
  outfile <- paste(outdir, "/", "rescued_from_discordance.png", sep = "")
  ggsave(outfile, width= 8, height = 8)
  outfile <- paste(outdir, "/", "rescued_from_discordance.svg", sep = "")
  ggsave(outfile, width= 8, height = 8)
  ggplot(df %>% filter(df$RESCUEDFROMABSENCE == 1), mapping = aes(x=LENGTH, y= READ_SUPPORT)) +
    scale_fill_gradient(low="lightblue1", high = "darkblue",trans = "log2") +
    geom_bin2d() +
    ggtitle("Formerly Absent Variants in HG002") +
    xlab("Variant Length") +
    ylab("Read Support") +
    theme(plot.title = element_text(size = 18, hjust = 0.5),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16)
    )
  outfile <- paste(outdir, "/", "rescued_from_absence.png", sep = "")
  ggsave(outfile, width= 8, height = 8)
  outfile <- paste(outdir, "/", "rescued_from_absence.svg", sep = "")
  ggsave(outfile, width= 8, height = 8)
  
  df$ParentSupport <- max(df$PARENT1SUPPORT, df$PARENT2SUPPORT)
  df$ParentSupport <- ifelse(df$PARENT1SUPPORT == -1, df$PARENT2SUPPORT, df$ParentSupport)
  df$ParentSupport <- ifelse(df$PARENT2SUPPORT == -1, df$PARENT1SUPPORT, df$ParentSupport)
  ggplot(df%>%filter(df$Discordance != "Discordant"), mapping = aes(x=ParentSupport, y = 1, fill = ParentSupport)) + 
    scale_fill_viridis_c(option = 'viridis') +
    geom_bar(stat = "identity") +
    xlim(0, 50) +
    ggtitle("Parent Read Support (All Variants)") +
    xlab("Parent Read Support") +
    ylab("Number of Variants") +
    theme(plot.title = element_text(size = 18, hjust = 0.5),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          legend.position = "none",
    )
  outfile <- paste(outdir, "/", "all_vars_parent_support.png", sep = "")
  ggsave(outfile, width= 8, height = 8)
  outfile <- paste(outdir, "/", "all_vars_parent_support.svg", sep = "")
  ggsave(outfile, width= 8, height = 8)
  bad <- df%>%filter(RESCUEDFROMABSENCE == 1 & ParentSupport < 3)
  head(bad, 3)
  ggplot(df%>%filter(df$RESCUEDFROMABSENCE == 1), mapping = aes(x=ParentSupport,y=1, fill = ParentSupport)) + 
    geom_bar(stat = "identity") +
    xlim(0, 50) +
    scale_fill_viridis_c(option = 'viridis') +
    ggtitle("Parent Read Support (Formerly Absent Variants)") +
    xlab("Parent Read Support") +
    ylab("Number of Variants") +
    theme(plot.title = element_text(size = 18, hjust = 0.5),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          legend.position = "none",
    )
  outfile <- paste(outdir, "/", "formerly_absent_parent_support.png", sep = "")
  ggsave(outfile, width= 8, height = 8)
  outfile <- paste(outdir, "/", "formerly_absent_parent_support.svg", sep = "")
  ggsave(outfile, width= 8, height = 8)
}

projectroot <- here()
projectroot <- "/home/mkirsche/jasmine_data"
wildcard <- "jasmine_md*.merged.tsv"
indir <- paste(projectroot, "/figures/figure2/", sep = '')
outdir <- indir
plotdisc(outdir, indir, wildcard)

# Read in Jasmine HG002 Trio HiFi merging results
infile <- paste(projectroot, "/figures/figure2/jasmine_md50.merged.tsv", sep = '')
df <- read.table(infile, sep = "\t", header = TRUE)

# Produce sample presence histograms without length filter
outfile <- paste(projectroot, "/figures/figure2/suppvechist.png", sep = '')
suppvec_hist_highlightdiscordant(df, "Jasmine", outfile, TRUE, 1700, .92, .85, 'Jasmine', 20, FALSE)
outfile <- paste(projectroot, "/figures/figure2/suppvechist.svg", sep = '')
suppvec_hist_highlightdiscordant(df, "Jasmine", outfile, TRUE, 1700, .86, .85, 'Jasmine', 20, FALSE)

# Produce sample presence histograms with length filter
outfile <- paste(projectroot, "/figures/figure2/suppvechist_50plus.png", sep = '')
suppvec_hist_highlightdiscordant(df, "Jasmine", outfile, TRUE, 1700, .92, .85, 'Jasmine', 20, TRUE)
outfile <- paste(projectroot, "/figures/figure2/suppvechist_50plus.svg", sep = '')
suppvec_hist_highlightdiscordant(df, "Jasmine", outfile, TRUE, 1700, .86, .85, 'Jasmine', 20, TRUE)

# Plot length distribution as a histogram
outfile <- paste(projectroot, "/figures/figure2/indels.png", sep = '')
plot_length(df, "Jasmine", outfile, TRUE, .95, .86)
outfile <- paste(projectroot, "/figures/figure2/indels.svg", sep = '')
plot_length(df, "Jasmine", outfile, TRUE, .92, .85)

# Plot length distribution as a line plot
outfile <- paste(projectroot, "/figures/figure2/indels_line.svg", sep = '')
plot_length_line(df, "Jasmine", outfile, TRUE, .92, .85)

# Read in default pipeline HG002 Trio HiFi merging results
infile <- paste(projectroot, "/default/hg002_trio_default.merged.tsv", sep = '')
df <- read.table(infile, sep = "\t", header = TRUE)

outfile <- paste(projectroot, "/figures/figure2/defaultdiscordance.png", sep = '')
suppvec_hist_highlightdiscordant(df, "Default Pipeline", outfile, FALSE, 1100, .9, .85, 'Prior Methods', 20, FALSE)
outfile <- paste(projectroot, "/figures/figure2/defaultdiscordance.svg", sep = '')
suppvec_hist_highlightdiscordant(df, "Default Pipeline", outfile, FALSE, 1100, .88, .85, 'Prior Methods', 20, FALSE)

outfile <- paste(projectroot, "/figures/figure2/defaultdiscordance_50plus.png", sep = '')
suppvec_hist_highlightdiscordant(df, "Default Pipeline", outfile, FALSE, 1100, .9, .85, 'Prior Methods', 20, TRUE)
outfile <- paste(projectroot, "/figures/figure2/defaultdiscordance_50plus.svg", sep = '')
suppvec_hist_highlightdiscordant(df, "Default Pipeline", outfile, FALSE, 1100, .88, .85, 'Prior Methods', 20, TRUE)

indir <- "/home/mkirsche/jasmine_data/figures/figure2/"
prefix <- "hg002_hifi"
discvec <- "100"
outdir <- indir
legendx <- .81
legendy <- .9
width <- 14
height <-10
plot_all_callers(indir, "hg002_hifi", "100", legendx,legendy, width, height, FALSE)
plot_all_callers(indir, "hg002_hifi", "100", legendx,legendy, width, height, TRUE)

infile <- paste(projectroot, "/figures/figure2/jasmine_md50.specificity.tsv", sep = '')
outdir <- paste(projectroot, "/figures/figure2", sep = '')
legendx <- .72
legendy <- .88
plotspec(infile, outdir, .72, .88)
  
