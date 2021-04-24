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

projectroot <- here()

directory = paste(projectroot, "/figures/supplement/jasmine_dist/", sep = "")
fn1 <- paste(directory, "jasmine_manhattan.merged.tsv", sep = "")
fn2 <- paste(directory, "jasmine_euclidean.merged.tsv", sep = "")
outfile <- paste(directory, "distance_function.png", sep = "")
df1 <- read.table(fn1, sep = "\t", header = T) 
df2 <- read.table(fn2, sep = "\t", header = T)
df1 <- df1 %>% filter(IS_SPECIFIC == "1" & IS_PRECISE == "1")
df2 <- df2 %>% filter(IS_SPECIFIC == "1" & IS_PRECISE == "1")
nrow(df1)
nrow(df2)

colnames(df1)

merged <- merge(df1, df2, by = c("SUPP_VEC", "IDLIST"), all = TRUE)
nrow(merged)
euclideanonly <- nrow(merged %>% filter(is.na(CHROM.x)))
manhattanonly <- nrow(merged %>% filter(is.na(CHROM.y)))
shared <- nrow(merged %>% filter(!is.na(CHROM.x) & !is.na(CHROM.y)))

outfile
Cairo(600, 600, file = outfile, type = "png", bg = "white")
grid.newpage()
venn.plot <- draw.pairwise.venn(area1           = euclideanonly + shared,
                              area2           = manhattanonly + shared,
                              cross.area      = shared,
                              category        = c('Euclidean', 'Manhattan'),
                              fill            = c('dodgerblue', 'orchid3'),
                              cat.col         = c('dodgerblue', 'orchid3'),
                              cex             = 2,
                              cat.cex         = 2,
                              cat.pos = c(-15, 15),
                              cat.dist = 0.05,
                              #cat.dist= c(.1, .1),
                              euler           = TRUE,
                              scaled          = FALSE
)
dev.off()
disc <- merged %>% filter(SUPP_VEC == "100")
disc1 <- disc %>% filter(is.na(CHROM.y))
disc2 <- disc %>% filter(is.na(CHROM.x))
nrow(disc)
nrow(disc1)
nrow(disc2)
disc1ids <- disc1 %>% select(CHROM.x, POS.x, ID.x)
disc2ids <- disc2 %>% select(CHROM.y, POS.y, ID.y)

outfile <- paste(directory, "discordant_euclidean.tsv", sep = "")
write.table(disc1ids, outfile, sep = "\t", quote = FALSE)

outfile <- paste(directory, "discordant_manhattan.tsv", sep = "")
write.table(disc2ids, outfile, sep = "\t", quote = FALSE)

