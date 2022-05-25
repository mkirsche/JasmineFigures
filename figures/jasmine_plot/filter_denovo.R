library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)
library(reshape2)
library(gridExtra)
library(rlist)
#library(ggfortify)
library(tidyverse)
library(plotly)

projectroot <- '/home/mkirsche/jasmine_data'#here()
wildcard <- "jasmine_md*.merged.tsv"
indir <- paste(projectroot, "/figures/figure4/", sep = '')
infile <- paste(indir, "allcandidates.tsv", sep = '')

df <- read.table(infile, sep = "\t", header = TRUE)

df$BREAKPOINT_STD <- (df$STD_quant_start + df$STD_quant_stop) / 2.0
df$ABS_LEN <- log(abs(df$SVLEN), base = 2)
df$DENOVO = 0
df$DENOVO <- ifelse(df$POS == 53340465 & df$CHROM == 'chr17', 1, df$DENOVO)
df$DENOVO <- ifelse(df$POS == 62805217 & df$CHROM == 'chr18', 1, df$DENOVO)
df$DENOVO <- ifelse(df$POS == 85552367 & df$CHROM == 'chr3', 1, df$DENOVO)
df$DENOVO <- ifelse(df$POS == 97089276 & df$CHROM == 'chr5', 1, df$DENOVO)
df$DENOVO <- ifelse(df$POS == 125785998 & df$CHROM == 'chr8', 1, df$DENOVO)
df$DENOVO <- ifelse(df$POS == 23280711 & df$CHROM == 'chr14', 1, df$DENOVO)
df$DENOVO
df$shp <- 'circle'
df$shp <- ifelse(df$DENOVO == 1, 'diamond', df$shp)
df$shp
fontsize <- 20
df$shp <- ifelse(df$POS == 23280711 & df$CHROM == 'chr14', 'square', df$shp)
#df$shp <- "o"
df$BreakpointVariance <- df$BREAKPOINT_STD*df$BREAKPOINT_STD
df$LOGRE = log(df$RE, base = 2)
ggplot(df, aes(x = ABS_LEN, y = LOGRE, color = as.character(DENOVO), size = BreakpointVariance)) + geom_point() + scale_color_manual(name = "De novo", values = c('#CC6677', 'darkgreen')) + xlab("Variant Length (log2)") + ylab("Read Support (log2)")
ggsave(paste(indir, "denovoscatter.svg", sep = ''), width = 6, height = 6)

fig <- plot_ly() %>% add_trace(df, x = df$RE, y = df$BREAKPOINT_STD, z = df$ABS_LEN, color = df$DENOVO, mode = "markers", colors = c('#CC6677', 'green'), marker = list(symbol = df$shp))
fig <- fig %>% layout(
  scene = list(
    xaxis = list(title = list(text = "<b>Read Support</b>", font = list(size = fontsize)), tickfont = list(size = fontsize - 4)),
    yaxis = list(autorange = "reversed", title = list(text = "<b>Breakpoint StDev</b>", font = list(size = fontsize)), tickfont = list(size = fontsize - 4)),
    zaxis = list(title = list(text = "<b>Length (log2)</b>", font = list(size = fontsize)), tickfont = list(size = fontsize - 4))
  ),
  
)
fig


