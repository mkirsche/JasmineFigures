library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)
library(reshape2)
library(gridExtra)
library(rlist)

fn <- "/home/mkirsche/jasmine_data/figures/crosstech/table.txt"
data = read.table(fn, comment.char = "#", sep = "\t", header = T, stringsAsFactors=FALSE)
data
data$DISC <- data$CLR_ONLY / (data$CLR_ONLY + data$CLR_ONLY)
ggplot(data, aes(x = MERGE_DIST, y = CLR_ONLY)) + geom_point()
ggsave("/home/mkirsche/jasmine_data/figures/crosstech/clr_scatter.png")
ggplot(data, aes(x = MERGE_DIST, y = HIFI_ONLY)) + geom_point()
ggsave("/home/mkirsche/jasmine_data/figures/crosstech/hifi_scatter.png")
ggplot(data, aes(x = MERGE_DIST, y = CLR_ONLY)) + geom_point()

ggplot(data %>% gather(TYPE, VALUE, -MERGE_DIST, - DISC), aes(x = MERGE_DIST, y = VALUE, fill = TYPE)) + 
  geom_bar(stat = "identity") +
  theme_bw()
ggsave("/home/mkirsche/jasmine_data/figures/crosstech/supp_vec_barplot.png")
