library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)
library(reshape2)
library(gridExtra)
library(rlist)
library(VennDiagram)
library(Cairo)
library(here)

samples <- c(5,
             10,
             20,
             30,
             40,
             50,
             100,
             250,
             500,
             750,
             1000,
             1500,
             2000,
             2504)

times<- c(5.201,
          7.245,
          13.398,
          18.532,
          24.710,
          30.860,
          80.084,
          350.832,
          1122.906,
          2866.614,
          4420.475,
          11306.869,
          18655.551,
          30792.497)

mems <- c(1.363,
          1.372,
          1.405,
          1.454,
          1.857,
          1.548,
          2.062,
          2.588,
          4.004,
          6.574,
          6.709,
          6.832,
          8.187,
          10.280
)

df <- data.frame(matrix(nrow = 14, ncol = 3))
df$Samples <- samples
df$Time <- times
df$Memory<- mems

ggplot(df, aes(x = Samples, y = Time)) + geom_line() + geom_point() + ylab("Time (s)") + theme(
  axis.text = element_text(size = 18), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20)
)
ggsave("/home/mkirsche/jasmine_data/revisions/simtime.png", width = 8, height = 8)
ggsave("/home/mkirsche/jasmine_data/revisions/simtime.svg", width = 8, height = 8)


ggplot(df, aes(x = Samples, y = Memory)) + geom_line() + geom_point() + ylab("Memory (GB)") + theme(
  axis.text = element_text(size = 18), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20)
)
ggsave("/home/mkirsche/jasmine_data/revisions/simmemory.png", width = 8, height = 8)
ggsave("/home/mkirsche/jasmine_data/revisions/simmemory.svg", width = 8, height = 8)
