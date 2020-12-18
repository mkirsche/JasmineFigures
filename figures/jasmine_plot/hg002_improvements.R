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

plot_table <- function(fn, outdir, prefix, caller, discvec){
  
  df <- read.table(fn, sep = "\t", header = TRUE)
  filtered = df %>% filter(SPECIFIC_FLAG == 1 & PRECISE_FLAG == 1 & SUPP_VEC == discvec)
  filtered$caller <- caller
  
  totalcount <- sum(df$NUMVARS)
  mergedcount <- nrow(df)
  return (list(filtered, totalcount, mergedcount))
}

plot_all_callers <- function(outdir, prefix, discvec){
  
  jasmineresults <- plot_table(paste(outdir, prefix, ".jasmine_augmented.txt", sep = ""), outdir, prefix, "Jasmine", discvec)
  survivorresults <- plot_table(paste(outdir, prefix, ".survivor_augmented.txt", sep = ""), outdir, prefix, "Survivor", discvec)
  svtoolsresults <- plot_table(paste(outdir, prefix, ".svtools_augmented.txt", sep = ""), outdir, prefix, "svtools", discvec)
  svimmerresults <- plot_table(paste(outdir, prefix, ".svimmer_augmented.txt", sep = ""), outdir, prefix, "svimmer", discvec)
  jasmineintraresults <- plot_table(paste(outdir, prefix, ".jasmineintra_augmented.txt", sep = ""), outdir, prefix, "Jasmine_Intrasample", discvec)
  
  discordant <- data.frame()
  discordant <- rbind(discordant, jasmineresults[[1]])
  discordant <- rbind(discordant, survivorresults[[1]])
  discordant <- rbind(discordant, svtoolsresults[[1]])
  discordant <- rbind(discordant, svimmerresults[[1]])
  discordant <- rbind(discordant, jasmineintraresults[[1]])
  
  errors <- read.table(paste(outdir, prefix, ".errors.txt", sep = ""), sep="\t", header = TRUE)
  errors$caller <- errors$SOFTWARE
  errors$caller <- ifelse(errors$caller == "survivor", "Survivor", as.character(errors$caller))
  errors$caller <- ifelse(errors$caller == "jasmine", "Jasmine", as.character(errors$caller))
  errors$caller <- ifelse(errors$caller == "jasmineintra", "Jasmine_Intrasample", as.character(errors$caller))
  errors$totaldisc <- errors$DISC_MIXED_STRAND_AND_TYPE + errors$DISC_MIXED_STRAND_ONLY + errors$DISC_MIXED_TYPE_ONLY
  errors$mixedtype <- errors$DISC_MIXED_STRAND_AND_TYPE + errors$DISC_MIXED_TYPE_ONLY
  errors$mixedstrand <- errors$DISC_MIXED_STRAND_ONLY
  
  #totals <- c(jasmineresults[[2]], survivorresults[[2]], svtoolsresults[[2]], svimmerresults[[2]], jasmineintraresults[[2]],
  #            jasmineresults[[3]], survivorresults[[3]], svtoolsresults[[3]], svimmerresults[[3]], jasmineintraresults[[3]])
  
  disccounts <- discordant %>% group_by(caller) %>% summarise(counts=n(), sums=sum(NUMVARS))
  disccounts
  
  disccounts <- inner_join(disccounts, errors, by = "caller")
  mdfr <- melt(disccounts, id.vars = "caller")
  
  ggplot(disccounts, aes(x = caller, y = counts, fill = caller)) +
    geom_bar(position = "stack", stat = "identity") +
    labs(title = paste("Discordant Variants (", prefix, ")")) +
    xlab("Merging Software") +
    ylab("Discordant Count") +
    theme(plot.title = element_text(size = 18, hjust = 0.5),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16),
    ) + theme(legend.position = "none")
  outfile <- paste(outdir, prefix, "_discordant.png", sep="")
  ggsave(outfile, width= 8, height = 8)
  
  ggplot(mdfr %>% filter(variable == "counts" | variable == "mixedtype" | variable == "mixedstrand"), aes(x = caller, y = as.numeric(value), fill = variable)) +
    geom_bar(stat = "identity") +
    labs(title = paste("Discordant and Invalid Variants in Child (", prefix, ")")) +
    xlab("Merging Software") +
    ylab("Count") +
    theme(plot.title = element_text(size = 18, hjust = 0.5),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16),
    ) +
    scale_fill_brewer(name = "Type", labels = c("Discordant", "Mixed Type", "Mixed Strand"), palette = "Set2")
  outfile <- paste(outdir, prefix, "_discordant_errors.png", sep="")
  ggsave(outfile, width= 8, height = 8)
}

plot_length <- function(df, caller, outfile, filter) 
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

suppvec_hist <- function(df, caller, outfile, filter) {
  if (filter)
  {
    df <- df %>% filter(IS_SPECIFIC == 1 & IS_PRECISE == 1)
  }
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
  suppveccounts$SVTYPE = "INS"
  ggplot(df, aes(x = SUPP_VEC_STRING, y = 1, fill = SVTYPE)) +
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
    theme(axis.title.x = element_blank())
  
  outfile <- paste(outdir, filterpref, "mddisc.png",sep="")
  ggsave(outfile, device = "png", height = 8, width = 8)
  
  totalplot <- ggplot(data = md_data, aes(x = md, y = Total)) +
    geom_point(size = 1.6) +
    ggtitle(paste("Total Variant Calls", sep = "")) +
    ylab("Total Calls") +
    theme(axis.title.x = element_blank())
  
  outfile <- paste(outdir, filterpref, "mdtotal.png",sep="")
  ggsave(outfile, device = "png", height = 8, width = 8)
  
  discrateplot <- ggplot(data = md_data, aes(x = md, y = DiscordanceRate)) +
    geom_point(size = 1.6) +
    ggtitle(paste("Discordance Rate", sep = "")) +
    xlab("Sniffles max_dist parameter")
  
  outfile <- paste(outdir, filterpref, "mddiscrate.png",sep="")
  ggsave(outfile, device = "png", height = 8, width = 8)
  ggarrange(arrangeGrob(discordantplot, totalplot, discrateplot, ncol = 1)) %>% ggexport(filename = paste(outdir, filterpref, "mdfullplot.png",sep = ""))
  
}

projectroot <- here()
wildcard <- "jasmine_md*.merged.tsv"
indir <- paste(projectroot, "/figures/figure2/", sep = '')
outdir <- indir
plotdisc(outdir, indir, wildcard)

infile <- paste(projectroot, "/figures/figure2/jasmine_md50.merged.tsv", sep = '')
outfile <- paste(projectroot, "/figures/figure2/suppvechist.png", sep = '')
df <- read.table(infile, sep = "\t", header = TRUE)
suppvec_hist(df, "Jasmine", outfile, TRUE)

outfile <- paste(projectroot, "/figures/figure2/indels.png", sep = '')
plot_length(df, "Hifi", outfile, TRUE)
