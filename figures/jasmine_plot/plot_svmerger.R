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

plot_table <- function(fn, outdir, prefix, caller, discvec){
  
  df <- read.table(fn, sep = "\t", header = TRUE)
  filtered = df %>% filter(SPECIFIC_FLAG == 1 & PRECISE_FLAG == 1)
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
  ggsave(paste(outdir, "/", prefix, "_", caller, ".numvars.png", sep = ""), width = 8, height = 8)
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
  ggsave(paste(outdir, "/", prefix, "_", caller, ".log.numvars.png", sep = ""), width = 8, height = 8)
  
  colorpalette <- c(DEL = '#CC6677',DUP = '#332288',INS = '#44AA99',TRA = '#88CCEE',INV = '#DDCC77')
  
  filtered$SUPP_VEC_STRING <- str_pad(as.character(filtered$SUPP_VEC), 3, "left", "0")
  filtered$SUPP_VEC_STRING <- ifelse(filtered$SUPP_VEC_STRING == "100", "Child Only", filtered$SUPP_VEC_STRING)
  filtered$SUPP_VEC_STRING <- ifelse(filtered$SUPP_VEC_STRING == "010", "Father Only", filtered$SUPP_VEC_STRING)
  filtered$SUPP_VEC_STRING <- ifelse(filtered$SUPP_VEC_STRING == "001", "Mother Only", filtered$SUPP_VEC_STRING)
  filtered$SUPP_VEC_STRING <- ifelse(filtered$SUPP_VEC_STRING == "110", "Son/Father", filtered$SUPP_VEC_STRING)
  filtered$SUPP_VEC_STRING <- ifelse(filtered$SUPP_VEC_STRING == "101", "Son/Mother", filtered$SUPP_VEC_STRING)
  filtered$SUPP_VEC_STRING <- ifelse(filtered$SUPP_VEC_STRING == "011", "Both parents", filtered$SUPP_VEC_STRING)
  filtered$SUPP_VEC_STRING <- ifelse(filtered$SUPP_VEC_STRING == "111", "All three", filtered$SUPP_VEC_STRING)
  
  suppveccounts <- filtered %>% count(SUPP_VEC_STRING)
  suppveccounts
  suppveccounts$TYPE = "INS"
  
  ggplot(filtered, aes(x = as.character(SUPP_VEC_STRING), y = 1, fill = TYPE)) +
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
    scale_fill_manual(name = "SVTYPE", values = colorpalette) +
    geom_text(data = suppveccounts, aes(x = SUPP_VEC_STRING, y=n, label=n), position=position_dodge(width=0.9), vjust=-0.75)
  
  ggsave(paste(outdir, "/", prefix, "_", caller, ".suppvechist.png", sep = ""), width= 7, height = 8)
  
  
  
  filtered = filtered %>% filter(SUPP_VEC == discvec)
  filtered$caller <- caller
  filtered$Intrasample = 'Allows Intrasample Merging'
  filtered$Intrasample <- ifelse(filtered$caller == "Jasmine", "No Intrasample Merging", as.character(filtered$Intrasample))
  filtered$Intrasample <- ifelse(filtered$caller == "dbsvmerge", "No Intrasample Merging", as.character(filtered$Intrasample))
  filtered$Intrasample <- ifelse(filtered$caller == "svpop", "No Intrasample Merging", as.character(filtered$Intrasample))
  
  
  
  
  totalcount <- sum(df$NUMVARS)
  mergedcount <- nrow(df)
  return (list(filtered, totalcount, mergedcount))
}

fn <- "/home/mkirsche/git/sv-merger/trio_svmerger_augmented.tsv"
results <- plot_table(fn, "/home/mkirsche/git/sv-merger", "hg002_hifi_trio", "svmerger", 100)
discordant <- results[[1]]
disccounts <- discordant %>% group_by(caller, Intrasample) %>% summarise(counts=n(), sums=sum(NUMVARS))
disccounts
