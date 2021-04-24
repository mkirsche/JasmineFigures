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
#projectroot <- "/home/mkirsche/jasmine_data"

directory = paste(projectroot, "/figures/supplement/jasmine_dist/", sep = "")
fn1 <- paste(directory, "jasmine_centroid.merged.tsv", sep = "")
fn2 <- paste(directory, "jasmine_clique.merged.tsv", sep = "")
fn3 <- paste(directory, "jasmine_edge.merged.tsv", sep = "")
df1 <- read.table(fn1, sep = "\t", header = T) 
df2 <- read.table(fn2, sep = "\t", header = T)
df3 <- read.table(fn3, sep = "\t", header = T)
df1 <- df1 %>% filter(IS_SPECIFIC == "1" & IS_PRECISE == "1")
df2 <- df2 %>% filter(IS_SPECIFIC == "1" & IS_PRECISE == "1")
df3 <- df3 %>% filter(IS_SPECIFIC == "1" & IS_PRECISE == "1")
nrow(df1)
nrow(df2)
nrow(df3)

colnames(df1)

merged <- merge(df1, df2, by = c("SUPP_VEC", "IDLIST"), all = TRUE)
merged <- merge(merged, df3, by = c("SUPP_VEC", "IDLIST"), all = TRUE)

#merged <- join_all(list(df1, df2, df3), by = c("SUPP_VEC", "IDLIST"))
nrow(merged)
colnames(merged)
centroidonly <- nrow(merged %>% filter(!is.na(CHROM.x) & is.na(CHROM.y) & is.na(CHROM)))
cliqueonly <- nrow(merged %>% filter(is.na(CHROM.x) & !is.na(CHROM.y) & is.na(CHROM)))
edgeonly <- nrow(merged %>% filter(is.na(CHROM.x) & is.na(CHROM.y) & !is.na(CHROM)))

allbutcentroid <- nrow(merged %>% filter(is.na(CHROM.x) & !is.na(CHROM.y) & !is.na(CHROM)))
allbutclique <- nrow(merged %>% filter(!is.na(CHROM.x) & is.na(CHROM.y) & !is.na(CHROM)))
allbutedge <- nrow(merged %>% filter(!is.na(CHROM.x) & !is.na(CHROM.y) & is.na(CHROM)))

allthree <- nrow(merged %>% filter(!is.na(CHROM.x) & !is.na(CHROM.y) & !is.na(CHROM)))

outfile <- paste(directory, "merging_requirements.png", sep = "")
outfile
Cairo(600, 600, file = outfile, type = "png", bg = "white")
grid.newpage()
venn.plot <- draw.triple.venn(area1           = centroidonly + allbutclique + allbutedge + allthree,
                                area2           = cliqueonly + allbutcentroid + allbutedge + allthree,
                                area3           = edgeonly + allbutcentroid + allbutclique + allthree,
                                n12             = allbutedge + allthree,
                                n23             = allbutcentroid + allthree,
                                n13             = allbutclique + allthree,
                                n123            = allthree,
                                category        = c('Centroid', 'Clique', 'Edge'),
                              fill            = c('dodgerblue', 'seagreen3', 'orchid3'),
                              cat.col         = c('dodgerblue', 'seagreen3', 'orchid3'),
                              cex             = 1,
                              cat.cex         = 1,
                              euler.d = FALSE,
                              scaled = FALSE
)
dev.off()
disc <- merged %>% filter(SUPP_VEC == "100")
disc1 <- disc %>% filter(is.na(CHROM.y) & is.na(CHROM))
disc2 <- disc %>% filter(is.na(CHROM.x) & is.na(CHROM))
disc3 <- disc %>% filter(is.na(CHROM.x) & is.na(CHROM.y))
disc12 <- disc %>% filter(!is.na(CHROM.x) & !is.na(CHROM.y) & is.na(CHROM))
disc13 <- disc %>% filter(!is.na(CHROM.x) & is.na(CHROM.y) & !is.na(CHROM))
disc23 <- disc %>% filter(is.na(CHROM.x) & !is.na(CHROM.y) & !is.na(CHROM))


nrow(disc)
nrow(disc1)
nrow(disc2)
nrow(disc3)
nrow(disc12)
nrow(disc13)
nrow(disc23)

disc1ids <- disc1 %>% select(CHROM.x, POS.x, ID.x)
disc2ids <- disc2 %>% select(CHROM.y, POS.y, ID.y)
disc3ids <- disc3 %>% select(CHROM, POS, ID)

outfile <- paste(directory, "discordant_centroid.tsv", sep = "")
write.table(disc1ids, outfile, sep = "\t", quote = FALSE)

outfile <- paste(directory, "discordant_clique.tsv", sep = "")
write.table(disc2ids, outfile, sep = "\t", quote = FALSE)

outfile <- paste(directory, "discordant_edge.tsv", sep = "")
write.table(disc3ids, outfile, sep = "\t", quote = FALSE)

