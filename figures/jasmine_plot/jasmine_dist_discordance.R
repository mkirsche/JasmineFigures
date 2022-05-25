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

plotdisc <- function(outdir, indir, wildcard) {
  paste(indir, wildcard, sep = "")
  md_files <- Sys.glob(paste(indir, wildcard, sep = ""))
  indir
  wildcard
    md_files
  
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
  md_data$md = as.numeric(gsub(".merged.tsv", '', substring(md_data$File, str_length(indir) + 13)))
  md_data$md
  
  md_data$DiscordanceRate = md_data$Discordant / md_data$Total
  
  filterpref = ""
  
  discordantplot <- ggplot(data = md_data, aes(x = md, y = Discordant)) +
    geom_point(size = 1.6) +
    ggtitle(paste("Discordant Variant Calls", sep = "")) +
    ylab("Discordant Calls") +
    theme(axis.title.x = element_blank())
  
  outfile <- paste(outdir, filterpref, "distdisc.png",sep="")
  ggsave(outfile, device = "png", height = 8, width = 8)
  
  totalplot <- ggplot(data = md_data, aes(x = md, y = Total)) +
    geom_point(size = 1.6) +
    ggtitle(paste("Total Variant Calls", sep = "")) +
    ylab("Total Calls") +
    theme(axis.title.x = element_blank())
  
  outfile <- paste(outdir, filterpref, "disttotal.png",sep="")
  ggsave(outfile, device = "png", height = 8, width = 8)
  
  discrateplot <- ggplot(data = md_data, aes(x = md, y = DiscordanceRate)) +
    geom_point(size = 1.6) +
    ggtitle(paste("Discordance Rate", sep = "")) +
    xlab("Jasmine min_dist parameter")
  
  outfile <- paste(outdir, filterpref, "distdiscrate.png",sep="")
  ggsave(outfile, device = "png", height = 8, width = 8)
  ggarrange(arrangeGrob(discordantplot, totalplot, discrateplot, ncol = 1)) %>% ggexport(filename = paste(outdir, filterpref, "distfullplot.png",sep = ""))
  
  ggarrange(arrangeGrob(discordantplot, totalplot, discrateplot, ncol = 1))
  ggsave(paste(outdir, filterpref, "distfullplot.svg", sep = ""), width = 8, height = 8)
  
}

plotdiscmd <- function(outdir, indir, wildcard, tech) {
  
  md_files <- Sys.glob(paste(indir, wildcard, sep = ""))
  names = c("File", "Total", "Discordant")
  md_data <- data.frame(setNames(rep(list(NA), length(names)), names))
  md_data$File = c()
  md_data$Total = c()
  md_data$Discordant = c()
  md_files
  md_fie <- "/home/mkirsche/jasmine_data/figures/supplement/jasmine_dist/jasmine_md_clr1000.merged.tsv"
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
  md_data$md = as.numeric(gsub(".merged.tsv", '', substring(md_data$File, str_length(indir) + 12 + str_length(tech))))
  md_data$md
  
  md_data$DiscordanceRate = md_data$Discordant / md_data$Total
  
  filterpref = ""
  
  discordantplot <- ggplot(data = md_data, aes(x = md, y = Discordant)) +
    geom_point(size = 1.6) +
    ggtitle(paste("Discordant Variant Calls", sep = "")) +
    ylab("Discordant Calls") +
    theme(axis.title.x = element_blank())
  
  outfile <- paste(outdir, filterpref, tech, "mddisc.png",sep="")
  ggsave(outfile, device = "png", height = 8, width = 8)
  
  totalplot <- ggplot(data = md_data, aes(x = md, y = Total)) +
    geom_point(size = 1.6) +
    ggtitle(paste("Total Variant Calls", sep = "")) +
    ylab("Total Calls") +
    theme(axis.title.x = element_blank())
  
  outfile <- paste(outdir, filterpref, tech, "mdtotal.png",sep="")
  ggsave(outfile, device = "png", height = 8, width = 8)
  
  discrateplot <- ggplot(data = md_data, aes(x = md, y = DiscordanceRate)) +
    geom_point(size = 1.6) +
    ggtitle(paste("Discordance Rate", sep = "")) +
    xlab("Sniffles max_dist parameter")
  
  outfile <- paste(outdir, filterpref, tech, "mddiscrate.png",sep="")
  ggsave(outfile, device = "png", height = 8, width = 8)
  ggarrange(arrangeGrob(discordantplot, totalplot, discrateplot, ncol = 1)) %>% ggexport(filename = paste(outdir, filterpref, tech, "mdfullplot.png",sep = ""))
  
  ggarrange(arrangeGrob(discordantplot, totalplot, discrateplot, ncol = 1))
  ggsave(paste(outdir, filterpref, tech, "mdfullplot.svg", sep = ""), width = 8, height = 8)
}
projectroot <- here()
projectroot <- '/home/mkirsche/jasmine_data'
wildcard <- "jasmine_dist*.merged.tsv"
indir <- paste(projectroot, "/figures/supplement/jasmine_dist/", sep = '')
outdir <- indir
plotdisc(outdir, indir, wildcard)

indir <- "/home/mkirsche/jasmine_data/figures/supplement/jasmine_dist/"
outdir <- indir
wildcard <- "jasmine_md_clr*.merged.tsv"
plotdiscmd(outdir, indir, wildcard, 'clr')

wildcard <- "jasmine_md_ont*.merged.tsv"
tech <- 'ont'
plotdiscmd(outdir, indir, wildcard, 'ont')
