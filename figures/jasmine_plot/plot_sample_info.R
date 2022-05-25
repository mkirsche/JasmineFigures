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
library(ggpattern)

fn <- "/home/mkirsche/jasmine_data/figures/figure5/filelist.txt"
df = read.table(fn, comment.char = "#", sep = "\t", header = T, stringsAsFactors=FALSE, colClasses = 'character')
colnames(df)
df$TECH <- ifelse(df$TECH == "PBCCS", "HiFi", df$TECH)
data <- split(df,f=df$SUPERPOPULATION)
data$AFR
p1 <- ggplot(df, aes(x = SAMPLE, y = as.numeric(COVERAGE), fill = TECH)) + 
  geom_bar(stat = "identity",
           position = "dodge",
           color = "black", 
           pattern_fill = "black",
           pattern_angle = 0,
           pattern_density = 0.05,
           pattern_spacing = 0.01,
           pattern_key_scale_factor = 0.6
  ) + 
  facet_grid(cols = vars(SUPERPOPULATION), scales = "free_x",space = "free_x") +
  guides(fill = guide_legend(override.aes = list(pattern = "none")), pattern = guide_legend(override.aes = list(fill = "lightgray"))) +
  xlab("Sample") +
  ylab("Coverage") +
  theme(axis.text.x = element_text(size = 14, angle = 30, vjust=.65),
        axis.text.y = element_text(size = 18),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 22),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 22),
        #legend.position = "none",
        strip.text.x = element_text(size = 18),
        strip.text.y = element_text(size = 18),
  ) +
  scale_fill_discrete(name = "Technology")

#p2 <- p1 %+% data$AMR
#p3 <- p1 %+% data$ASH
#p4 <- p1 %+% data$EAS
#p5 <- p1 %+% data$SAS
#ggarrange(p1, p2, p3, p4, p5)

ggsave("/home/mkirsche/jasmine_data/figures/figure5/datasummary.png", width= 14, height = 5)                            
ggsave("/home/mkirsche/jasmine_data/figures/figure5//datasummary.svg", width= 14, height = 5)                            

