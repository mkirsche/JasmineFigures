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

plot_iris_scores <- function(smallold, smallnew, mediumold, mediumnew, largeold, largenew, outfile) {
  smallolddf <- read.table(smallold)
  smallolddf$Refined <- 'Unrefined Variant Calls'
  smallolddf$Size <- 'Small (50 - 200 bp)'
  
  smallnewdf <- read.table(smallnew)
  smallnewdf$Refined <- 'Refined Variant Calls'
  smallnewdf$Size <- 'Small (50 - 200 bp)'
  
  mediumolddf <- read.table(mediumold)
  mediumolddf$Refined <- 'Unrefined Variant Calls'
  mediumolddf$Size <- 'Medium (900 - 1100 bp)'
  
  mediumnewdf <- read.table(mediumnew)
  mediumnewdf$Refined <- 'Refined Variant Calls'
  mediumnewdf$Size <- 'Medium (900 - 1100 bp)'
  
  largeolddf <- read.table(largeold)
  largeolddf$Refined <- 'Unrefined Variant Calls'
  largeolddf$Size <- 'Large (4000 - 6000 bp)'
  
  largenewdf <- read.table(largenew)
  largenewdf$Refined <- 'Refined Variant Calls'
  largenewdf$Size <- 'Large (4000 - 6000 bp)'
  
  combineddf <- rbind(smallolddf, smallnewdf, mediumolddf, mediumnewdf, largeolddf, largenewdf)
  combineddf$Score <- combineddf$V1
  colnames(combineddf)
  
  combineddf$Refined = factor(combineddf$Refined, levels = c("Unrefined Variant Calls", "Refined Variant Calls"))
  combineddf$Size = factor(combineddf$Size, levels = c("Small (50 - 200 bp)", "Medium (900 - 1100 bp)", "Large (4000 - 6000 bp)"))

  combinedplot <- ggplot(combineddf, aes(x = Score)) + geom_histogram(breaks = seq(0.75, 1.0, 0.01), color="black",fill="lightblue2") +
    xlab('Insertion Sequence Accuracy') +
    ylab('Number of Variants') +
    theme(plot.title = element_text(size = 18, hjust = 0.5),
          axis.text.x = element_text(size = 12, angle = 30, margin = margin(t = 15)),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16), 
          strip.text.x = element_text(
            size = 18, face = "bold.italic"
          ),
          strip.text.y = element_text(
            size = 14, face = "bold.italic"
          )
    ) +
    facet_grid(cols = vars(Refined), rows = vars(Size))
  ggsave(outfile, plot = combinedplot, width= 10, height = 8)
  
  combinedplot <- ggplot(combineddf, aes(x = Score)) + geom_histogram(breaks = seq(0.75, 1.0, 0.01), color="black",fill="lightblue2") +
    xlab('Insertion Sequence Accuracy') +
    ylab('Number of Variants') +
    theme(plot.title = element_text(size = 18, hjust = 0.5),
          axis.text.x = element_text(size = 12, angle = 30, margin = margin(t = 15)),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16), 
          strip.text.x = element_text(
            size = 18, face = "bold.italic"
          ),
          strip.text.y = element_text(
            size = 14, face = "bold.italic"
          )
    ) +
    facet_grid(cols = vars(Refined))
  ggsave(outfile, plot = combinedplot, width= 10, height = 8)

}


indir <- "/home/mkirsche/jasmine_data/figures/supplement/iris"
smallold <- paste(indir,"/iris_sim_small.scores.txt", sep = "")
smallnew <- paste(indir, "/iris_sim_small.refined.scores.txt", sep = "")
mediumold <- paste(indir,"/iris_sim_medium.scores.txt", sep = "")
mediumnew <- paste(indir, "/iris_sim_medium.refined.scores.txt", sep = "")
largeold <- paste(indir,"/iris_sim_large.scores.txt", sep = "")
largenew <- paste(indir, "/iris_sim_large.refined.scores.txt", sep = "")

outfile <- "/home/mkirsche/jasmine_data/figures/supplement/iris/iris_sim_hist.png"
plot_iris_scores(smallold, smallnew, mediumold, mediumnew, largeold, largenew, outfile)


refinedscores <- "/home/mkirsche/jasmine_data/figures/supplement/iris/hg002_refined_scores.txt"
refined <- read.table(refinedscores, sep = '\t', header = TRUE)
colnames(refined)
nrow(refined)
nrow(refined %>% filter(REFINED == 1))
nrow(refined %>% filter(REFINED == 0))
refined <- refined[order(REFINED)]
refined$REFINED
refined$NUMSEQID <- as.numeric(refined$SEQ_IDENTITY)
refined$REFINEDCHAR <- as.character(refined$REFINED)
colnames(refined)
ggplot(refined) + geom_histogram(position = "stack", binwidth = 0.01, aes(x = SEQ_IDENTITY, fill = REFINEDCHAR))
ggsave("/home/mkirsche/jasmine_data/figures/supplement/iris/hg002_refined_scores.png", width = 8, height = 8)

unrefinedscores <- "/home/mkirsche/jasmine_data/figures/supplement/iris/hg002_unrefined_scores.txt"
unrefined <- read.table(unrefinedscores, sep = '\t', header = TRUE)
nrow(unrefined)
mutual_set_refined <- refined %>% filter(ID %in% unrefined$ID)
mutual_set_unrefined <- unrefined %>% filter(ID %in% refined$ID)
nrow(mutual_set_refined)
nrow(mutual_set_unrefined)
mean(mutual_set_refined$SEQ_IDENTITY)
mean(mutual_set_unrefined$SEQ_IDENTITY)

good <- mutual_set_refined %>% filter(REFINED == 1)
bad <- mutual_set_refined %>% filter(REFINED == 0)
mean(good$NUMSEQID)
mean(bad$NUMSEQID)

mutual_set_unrefined_changed <- unrefined %>% filter(ID %in% good$ID)
nrow(unrefined)
nrow(mutual_set_unrefined_changed)
mean(good$NUMSEQID)
mean(mutual_set_unrefined_changed$SEQ_IDENTITY)

goodplot <- ggplot(good) + geom_histogram(position = "stack", binwidth = 0.005, fill = "lightblue", aes(x = SEQ_IDENTITY)) + xlab("Sequence Identity") + ylab("Iris Refined") + 
  theme(axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16), axis.text = element_text(size = 12))
unrefinedgoodplot <- ggplot(mutual_set_unrefined_changed) + geom_histogram(position = "stack", fill = "lightblue", binwidth = 0.005, aes(x = SEQ_IDENTITY)) + ylab("Raw SV Calls") + 
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 16), axis.text = element_text(size = 12))
good$Refined <- "Refined Variant Calls"
mutual_set_unrefined_changed$Refined <- "Unrefined Variant Calls"
ncol(good)
colnames(good)
colnames(mutual_set_unrefined_changed)
ncol(mutual_set_unrefined_changed)
mutual_set_unrefined_changed$NUMSEQID <- "0"
mutual_set_unrefined_changed$REFINEDCHAR <- "0"

all <- rbind(good, mutual_set_unrefined_changed)
all$Refined = factor(all$Refined, levels = c("Unrefined Variant Calls", "Refined Variant Calls"))
combinedplot <- ggplot(all, aes(x = SEQ_IDENTITY)) + geom_histogram(breaks = seq(0.50, 1.0, 0.01), color="black",fill="lightblue2") +
  xlab('Insertion Sequence Accuracy') +
  ylab('Number of Variants') +
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        axis.text.x = element_text(size = 12, angle = 30, margin = margin(t = 15)),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16), 
        strip.text.x = element_text(
          size = 18, face = "bold.italic"
        ),
        strip.text.y = element_text(
          size = 14, face = "bold.italic"
        )
  ) +
  facet_grid(cols = vars(Refined))
ggsave("/home/mkirsche/jasmine_data/figures/supplement/iris/irisrealscores.svg", plot = combinedplot, width= 10, height = 8)
ggsave("/home/mkirsche/jasmine_data/figures/supplement/iris/irisrealscores.png", plot = combinedplot, width= 10, height = 8)

#ggarrange(arrangeGrob(unrefinedgoodplot, goodplot), ncol = 1, nrow = 2)
#ggsave(filename = )
