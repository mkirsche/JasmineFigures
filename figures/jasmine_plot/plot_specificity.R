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

plotspec <- function(infile, outdir, legendx, legendy) {
  df <- read.table(infile, sep = "\t", header = TRUE)
  df
  df$LENGTH
  df <- df %>% filter(LENGTH < 100 & READ_SUPPORT < 50)
  #df <- df[sample(nrow(df), 2000), ]
  df$Discordance <- ifelse(df$DISCORDANT == 0, "Not discordant", "Discordant")
  df$RESCUEDFROMABSENCE
  df$Discordance <- ifelse(df$RESCUEDFROMABSENCE == 1, "Rescued from absence", df$Discordance)
  df$Discordance <- ifelse(df$RESCUEDFROMDISCORDANCE == 1, "Rescued from discordance", df$Discordance)
  ggplot(df, mapping = aes(x = READ_SUPPORT, y = 1, fill = Discordance)) + 
    scale_fill_discrete() +
    xlab("Read Support") +
    ylab ("Number of Variants") +
    ggtitle("Effects of Double Threshold") +
    theme(plot.title = element_text(size = 18, hjust = 0.5),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          legend.position = c(legendx, legendy),
    ) +
    geom_bar(position = "stack", stat = "identity") +
    scale_shape_identity()
  outfile <- paste(outdir, "/", "specificity_readsupp.png", sep = "")
  ggsave(outfile, width= 8, height = 8)
  
  ggplot(df, mapping = aes(x = LENGTH, y = 1, fill = Discordance)) + 
    scale_fill_discrete() +
    xlab("Variant Length") +
    ylab ("Number of Variants") +
    ggtitle("Effects of Double Threshold") +
    theme(plot.title = element_text(size = 18, hjust = 0.5),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          legend.position = c(legendx, legendy),
    ) +
    geom_bar(position = "stack", stat = "identity") +
    scale_shape_identity()
  outfile <- paste(outdir, "/", "specificity_length.png", sep = "")
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
  
  df$ParentSupport <- max(df$PARENT1SUPPORT, df$PARENT2SUPPORT)
  df$ParentSupport <- ifelse(df$PARENT1SUPPORT == -1, df$PARENT2SUPPORT, df$ParentSupport)
  df$ParentSupport <- ifelse(df$PARENT2SUPPORT == -1, df$PARENT1SUPPORT, df$ParentSupport)
  ggplot(df%>%filter(df$Discordance != "Discordant"), mapping = aes(x=ParentSupport, y = 1)) + 
    #scale_fill_viridis_c(option = 'viridis') +
    geom_bar(stat = "identity", fill = "#2093c3") +
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
  
  ggplot(df%>%filter(df$Discordance != "Discordant"), mapping = aes(x=ParentSupport, y = 1)) + 
    #scale_fill_viridis_c(option = 'viridis') +
    geom_bar(stat = "identity", fill = "#2093c3") +
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
  outfile <- paste(outdir, "/", "all_vars_parent_support.svg", sep = "")
  ggsave(outfile, width= 8, height = 8)
  
  
  bad <- df%>%filter(RESCUEDFROMABSENCE == 1 & ParentSupport < 3)
  head(bad, 3)
  ggplot(df%>%filter(df$RESCUEDFROMABSENCE == 1), mapping = aes(x=ParentSupport,y=1)) + 
    geom_bar(stat = "identity", fill = "#2093c3") +
    xlim(0, 50) +
    #scale_fill_viridis_c(option = 'viridis') +
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
  
  ggplot(df%>%filter(df$RESCUEDFROMABSENCE == 1), mapping = aes(x=ParentSupport,y=1)) + 
    geom_bar(stat = "identity", fill = "#2093c3") +
    xlim(0, 50) +
    #scale_fill_viridis_c(option = 'viridis') +
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
  outfile <- paste(outdir, "/", "formerly_absent_parent_support.svg", sep = "")
  ggsave(outfile, width= 8, height = 8)
}
projectroot <- "/home/mkirsche/jasmine_data"
infile <- paste(projectroot, "/figures/figure2/jasmine_md50.specificity.tsv", sep = '')
outdir <- paste(projectroot, "/figures/figure2", sep = '')
plotspec(infile, outdir, .86, .88)
