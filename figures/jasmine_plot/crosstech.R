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

suppvec_hist <- function(df, caller, outfile, filter, legendx, legendy, lengthfilter) {
  if (filter)
  {
    df <- df %>% filter(IS_SPECIFIC == 1 & IS_PRECISE == 1)
  }
  crosstechtitle <- "Variants by Sequencing Technology"
  if (lengthfilter) {
    df <- df %>% filter(SVLEN == 0 | abs(SVLEN) >= 50 | SVTYPE == "TRA")
    crosstechtitle <- "SVs by Sequencing Technology"
  }
  colorpalette <- c(DEL = '#CC6677',DUP = '#332288',INS = '#44AA99',TRA = '#88CCEE',INV = '#DDCC77')
  
  df$SUPP_VEC_STRING <- str_pad(as.character(df$SUPP_VEC), 3, "left", "0")
  df$SUPP_VEC_STRING <- ifelse(df$SUPP_VEC_STRING == "100", "Hi-Fi Only", df$SUPP_VEC_STRING)
  df$SUPP_VEC_STRING <- ifelse(df$SUPP_VEC_STRING == "010", "CLR Only", df$SUPP_VEC_STRING)
  df$SUPP_VEC_STRING <- ifelse(df$SUPP_VEC_STRING == "001", "ONT Only", df$SUPP_VEC_STRING)
  df$SUPP_VEC_STRING <- ifelse(df$SUPP_VEC_STRING == "110", "HiFi/CLR", df$SUPP_VEC_STRING)
  df$SUPP_VEC_STRING <- ifelse(df$SUPP_VEC_STRING == "101", "HiFi/ONT", df$SUPP_VEC_STRING)
  df$SUPP_VEC_STRING <- ifelse(df$SUPP_VEC_STRING == "011", "CLR/ONT", df$SUPP_VEC_STRING)
  df$SUPP_VEC_STRING <- ifelse(df$SUPP_VEC_STRING == "111", "All three", df$SUPP_VEC_STRING)
  
  suppveccounts <- df %>% count(SUPP_VEC_STRING)
  suppveccounts
  suppveccounts$SVTYPE = "INS"
  ggplot(df, aes(x = SUPP_VEC_STRING, fill = SVTYPE)) +
    geom_bar(position = "stack", stat = "count") +
    labs(title = paste(crosstechtitle, sep = "")) +
    xlab("Technologies") +
    ylab("Count") +
    scale_y_continuous(expand = expansion(mult = c(0, .1))) +
    theme(axis.text.x = element_text(size = 16, angle = 30, margin = margin(t = 21)),
          axis.text.y = element_text(size = 16),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 20),
          legend.text = element_text(size = 16),
          legend.title = element_text(size = 20),
          legend.position = c(legendx, legendy),
          plot.title = element_text(size = 22, hjust = .5)
    ) +
    scale_fill_manual(name = "SVTYPE", values = colorpalette) +
    guides(fill=guide_legend(title="Type")) +
    geom_text(data = suppveccounts, aes(x = SUPP_VEC_STRING, y=n, label=n), position=position_dodge(width=0.9), vjust=-0.75)
  
  ggsave(outfile, width= 7, height = 8)
}

tech_specific <- function(df, caller, outfile, filter, legendx, legendy) {
  if (filter)
  {
    df <- df %>% filter(IS_SPECIFIC == 1 & IS_PRECISE == 1)
  }
  typecounts <- df %>% count(SVTYPE)
  typecounts
  ggplot(df, aes(x = SVTYPE, y = 1)) +
    geom_bar(position = "stack", stat = "identity") +
    labs(title = paste("SVs by Type (", caller, ")")) +
    xlab("SV Type") +
    ylab("Count") +
    theme(plot.title = element_text(size = 18, hjust = 0.5),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 16),
          legend.position = c(legendx, legendy),
    ) +
    geom_text(data = typecounts, aes(x = SVTYPE, y=n, label=n), position=position_dodge(width=0.9), vjust=-0.75)
  
    ggsave(outfile, width= 6, height = 8)
}

tech_specific_length <- function(df, caller, outfile, filter, sz) 
{
  if (filter)
  {
    df <- df %>% filter(IS_SPECIFIC == 1 & IS_PRECISE == 1)
  }
  df$LenCategory = "TRA"
  colorpalette <- c(DEL = '#CC6677',DUP = '#332288',INS = '#44AA99',TRA = '#88CCEE',INV = '#DDCC77')
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
  unique(df$LenCategory)
  labellist = c("30-50", "50-100", "100-200", "200-300", "300-400", "400-500", "500-750", "750-1k", "1k-2k", "2k-5k", "5k-10k", "10k+", "TRA")
  df$LenCategory <- factor(df$LenCategory,levels = labellist)
  
  df %>% group_by(df$LenCategory) %>% tally()

  delcounts <- df %>% filter(df$SVTYPE == "DEL") %>% group_by(LenCategory) %>% count()
  delcounts$SVTYPE = "DEL"
  
  nondelcounts <- df %>% filter(df$SVTYPE != "DEL") %>% group_by(LenCategory) %>% count()
  nondelcounts$SVTYPE = "INS"
  summarized <- df %>% group_by(SVTYPE, LenCategory) %>% summarise(counts=n())
  ggplot(summarized) +
        geom_bar(data = summarized %>% filter(SVTYPE != "DEL"), position = "stack", stat = "identity", aes(x = LenCategory, fill = SVTYPE, y = counts))+
        geom_bar(data = summarized %>% filter(SVTYPE == "DEL"), position = "stack", stat = "identity", aes(x = LenCategory, fill = SVTYPE, y = -counts))+
        scale_x_discrete(labels=labellist) +
        xlab("Length") +
        ylab("Count") +
        labs(title = paste(caller, " Variant Size Distribution", sep = "")) +
        theme(plot.title = element_text(size = sz + 6, hjust = 0.5),
              axis.text.x = element_text(size = sz, angle = 30, margin = margin(t = 17, r = 0, b = 0, l = 0)),
              axis.text.y = element_text(size = sz),
              axis.title.x = element_text(size = sz + 2),
              axis.title.y = element_text(size = sz + 2),
              legend.text = element_text(size = sz),
              legend.title = element_text(size = sz + 2),
              text = element_text(size = sz + 2*(sz - 16)),
              legend.position = "none",
        ) +
        scale_fill_manual(name = "SVTYPE", values = colorpalette) +
        guides(fill=guide_legend(title="Type")) +
        geom_text(data = nondelcounts, aes(x = LenCategory, y=n, label=n), position=position_dodge(width=0.9), vjust=-0.5, size = (sz-4)*5/14) +
        geom_text(data = delcounts, aes(x = LenCategory, y=-1*n, label=n), position=position_dodge(width=0.9), vjust=1.25, size = (sz-4)*5/14)

  ggsave(outfile, width= 12, height = 8)
}

plot_length_line <- function(df, outfile, legendx, legendy) 
{
  df$ABSLEN = abs(as.numeric(df$SVLEN))
  
  colorpalette <- c(DEL = '#CC6677',DUP = '#332288',INS = '#44AA99',TRA = '#88CCEE',INV = '#DDCC77')
  
  #summarized <- df %>% group_by(SVTYPE) %>% summarise(counts=n())
  
  ggplot(df%>% filter(ABSLEN <=10000 & ABSLEN >= 100), aes(x = ABSLEN, color = SVTYPE)) +
    geom_density()+
    xlab("Length") +
    ylab("Density") +
    theme(plot.title = element_text(size = 24, hjust = 0.5),
          axis.text.x = element_text(size = 18, angle = 30, margin = margin(t = 17, r = 0, b = 0, l = 0)),
          axis.text.y = element_text(size = 18),
          axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20),
          legend.text = element_text(size = 18),
          legend.title = element_text(size = 20),
          text = element_text(size = 20),
          legend.position = c(legendx, legendy),
    ) +
    ggtitle('HG002 Trio Technology-Concordant Variant Sizes') +
    #scale_fill_brewer(name = "Type", palette = "Set2") + 
    scale_color_manual(name = "SVTYPE", values = colorpalette) +
    guides(color=guide_legend(title="Type"))
  #geom_text(data = nondelcounts, aes(x = LenCategory, y=n, label=n), position=position_dodge(width=0.9), vjust=-0.75) +
  #geom_text(data = delcounts, aes(x = LenCategory, y=-1*n, label=n), position=position_dodge(width=0.9), vjust=1.5)
  
  ggsave(outfile, width= 10, height = 8)
}

projectroot <- here()
projectroot <- '/home/mkirsche/jasmine_data'

fn <- paste(projectroot, '/figures/figure3/hg002_crosstech.merged.tsv', sep = '')
df <- read.table(fn, sep = "\t", header = TRUE)
df_longspecprec <- df <- df %>% filter(SVLEN == 0 | abs(SVLEN) >= 50 | SVTYPE == "TRA") %>% filter(IS_SPECIFIC == 1 & IS_PRECISE == 1)
n
caller <- 'Jasmine CrossTech'

outfile <- paste(projectroot, '/figures/figure3/hg002_crosstech_suppvecs.png', sep = '')
suppvec_hist(df, "Jasmine CrossTech", outfile, FALSE, .88, .85, FALSE)
outfile <- paste(projectroot, '/figures/figure3/hg002_crosstech_suppvecs_spec_prec.png', sep = '')
suppvec_hist(df, "Jasmine CrossTech", outfile, TRUE, .88, .85, FALSE)
outfile <- paste(projectroot, '/figures/figure3/hg002_crosstech_suppvecs_spec_prec.svg', sep = '')
suppvec_hist(df, "Jasmine CrossTech", outfile, TRUE, .88, .85, FALSE)

outfile <- paste(projectroot, '/figures/figure3/hg002_crosstech_suppvecs_50plus.png', sep = '')
suppvec_hist(df, "Jasmine CrossTech", outfile, FALSE, .88, .85, TRUE)
outfile <- paste(projectroot, '/figures/figure3/hg002_crosstech_suppvecs_spec_prec_50plus.png', sep = '')
suppvec_hist(df, "Jasmine CrossTech", outfile, TRUE, .88, .85, TRUE)
outfile <- paste(projectroot, '/figures/figure3/hg002_crosstech_suppvecs_spec_prec_50plus.svg', sep = '')
suppvec_hist(df, "Jasmine CrossTech", outfile, TRUE, .88, .85, TRUE)

# Hifi-unique variants
#outfile <- paste(projectroot, '/figures/figure3/hg002_crosstech_hifi.png', sep = '')
#tech_specific(df %>% filter(SUPP_VEC == 100), "Hifi Only", outfile, TRUE)
outfile <- paste(projectroot, '/figures/figure3/hg002_crosstech_hifi_indels.png', sep = '')
tech_specific_length(df %>% filter(SUPP_VEC == 100), "HiFi-Only", outfile, TRUE, 22)
outfile <- paste(projectroot, '/figures/figure3/hg002_crosstech_hifi_indels.svg', sep = '')
tech_specific_length(df %>% filter(SUPP_VEC == 100), "HiFi-Only", outfile, TRUE, 22)

# CLR-unique variants
#outfile <- paste(projectroot, '/figures/figure3/hg002_crosstech_clr.png', sep = '')
#tech_specific(df %>% filter(SUPP_VEC == 010), "CLR Only", outfile, TRUE)
outfile <- paste(projectroot, '/figures/figure3/hg002_crosstech_clr_indels.png', sep = '')
tech_specific_length(df %>% filter(SUPP_VEC == 010), "CLR-Only", outfile, TRUE, 22)
outfile <- paste(projectroot, '/figures/figure3/hg002_crosstech_clr_indels.svg', sep = '')
tech_specific_length(df %>% filter(SUPP_VEC == 010), "CLR-Only", outfile, TRUE, 22)

# ONT-unique variants
#outfile <- paste(projectroot, '/figures/figure3/hg002_crosstech_ont.png', sep = '')
#tech_specific(df %>% filter(SUPP_VEC == 001), "ONT Only", outfile, TRUE)
outfile <- paste(projectroot, '/figures/figure3/hg002_crosstech_ont_indels.png', sep = '')
tech_specific_length(df %>% filter(SUPP_VEC == 001), "ONT-Only", outfile, TRUE, 22)
outfile <- paste(projectroot, '/figures/figure3/hg002_crosstech_ont_indels.svg', sep = '')
tech_specific_length(df %>% filter(SUPP_VEC == 001), "ONT-Only", outfile, TRUE, 22)

# Concordant variants
#outfile <- paste(projectroot, '/figures/figure3/hg002_crosstech_agree.png', sep = '')
#tech_specific(df %>% filter(SUPP_VEC == 111), "Concordant", outfile, TRUE)
outfile <- paste(projectroot, '/figures/figure3/hg002_crosstech_agree_indels.png', sep = '')
tech_specific_length(df %>% filter(SUPP_VEC == 111), "Concordant", outfile, TRUE, 18)
outfile <- paste(projectroot, '/figures/figure3/hg002_crosstech_agree_indels.svg', sep = '')
tech_specific_length(df %>% filter(SUPP_VEC == 111), "Concordant", outfile, TRUE, 18)
outfile <- paste(projectroot, '/figures/figure3/hg002_crosstech_agree_indels_line.svg', sep = '')
plot_length_line(df %>% filter(SUPP_VEC == 111), outfile, .88, .85)


df <- df %>% filter(IS_SPECIFIC == 1 & IS_PRECISE == 1)
pb_only <- nrow(df %>% filter(SUPP_VEC == 010))
ont_only <- nrow(df %>% filter(SUPP_VEC == 001))
hifi_only <- nrow(df %>% filter(SUPP_VEC == 100))
pb_ont <- nrow(df %>% filter(SUPP_VEC == 011))
pb_hifi <- nrow(df %>% filter(SUPP_VEC == 110))
ont_hifi <- nrow(df %>% filter(SUPP_VEC == 101))
allthree <- nrow(df %>% filter(SUPP_VEC == 111))

tmp <- df%>%filter(SUPP_VEC==001)
nrow(tmp)
nrow(tmp%>%filter(SVTYPE=="INS"))

outfile <- paste(projectroot, '/figures/figure3/hg002_crosstech_venn.png', sep = '')
Cairo(600, 600, file = outfile, type = "png", bg = "white")
grid.newpage()
venn.plot <- draw.triple.venn(area1           = pb_only + pb_ont + pb_hifi + allthree,
                              area2           = ont_only + pb_ont + ont_hifi + allthree,
                              area3           = hifi_only + pb_hifi + ont_hifi + allthree,
                              n12             = pb_ont + allthree,
                              n23             = ont_hifi + allthree,
                              n13             = pb_hifi + allthree,
                              n123            = allthree,
                              category        = c('CLR', 'ONT', 'HiFi'),
                              fill            = c('#DD3D2D', '#FEDA8B', '#C2D4EF'),
                              #cat.col         = c('#DD3D2D', '#FEDA8B', '#C2D4EF'),#c('dodgerblue', 'seagreen3', 'orchid3'),
                              cex             = 2,
                              cat.cex         = 2,
                              euler           = TRUE,
                              scaled          = FALSE
)
dev.off()

outfile <- paste(projectroot, '/figures/figure3/hg002_crosstech_venn.pdf', sep = '')
pdf(outfile, useDingbats = FALSE, width = 6, height = 6)
grid.newpage()
venn.plot <- draw.triple.venn(area1           = pb_only + pb_ont + pb_hifi + allthree,
                              area2           = ont_only + pb_ont + ont_hifi + allthree,
                              area3           = hifi_only + pb_hifi + ont_hifi + allthree,
                              n12             = pb_ont + allthree,
                              n23             = ont_hifi + allthree,
                              n13             = pb_hifi + allthree,
                              n123            = allthree,
                              category        = c('CLR', 'ONT', 'HiFi'),
                              fill            = c('#DD3D2D', '#FEDA8B', '#C2D4EF'),
                              #cat.col         = c('#DD3D2D', '#FEDA8B', '#C2D4EF'),#c('dodgerblue', 'seagreen3', 'orchid3'),
                              cex             = 2,
                              cat.cex         = 2,
                              euler           = TRUE,
                              scaled          = FALSE
)
dev.off()



df <- df %>% filter(SVLEN == 0 | abs(SVLEN) >= 50 | SVTYPE == "TRA")
pb_only <- nrow(df %>% filter(SUPP_VEC == 010))
ont_only <- nrow(df %>% filter(SUPP_VEC == 001))
hifi_only <- nrow(df %>% filter(SUPP_VEC == 100))
pb_ont <- nrow(df %>% filter(SUPP_VEC == 011))
pb_hifi <- nrow(df %>% filter(SUPP_VEC == 110))
ont_hifi <- nrow(df %>% filter(SUPP_VEC == 101))
allthree <- nrow(df %>% filter(SUPP_VEC == 111))
outfile <- paste(projectroot, '/figures/figure3/hg002_crosstech_venn_50plus.pdf', sep = '')
pdf(outfile, useDingbats = FALSE, width = 6, height = 6)
grid.newpage()
venn.plot <- draw.triple.venn(area1           = pb_only + pb_ont + pb_hifi + allthree,
                              area2           = ont_only + pb_ont + ont_hifi + allthree,
                              area3           = hifi_only + pb_hifi + ont_hifi + allthree,
                              n12             = pb_ont + allthree,
                              n23             = ont_hifi + allthree,
                              n13             = pb_hifi + allthree,
                              n123            = allthree,
                              category        = c('CLR', 'ONT', 'HiFi'),
                              fill            = c('#DD3D2D', '#FEDA8B', '#C2D4EF'),
                              #cat.col         = c('#DD3D2D', '#FEDA8B', '#C2D4EF'),#c('dodgerblue', 'seagreen3', 'orchid3'),
                              cex             = 2,
                              cat.cex         = 2,
                              euler           = TRUE,
                              scaled          = FALSE
)
dev.off()

outfile <- paste(projectroot, '/figures/figure3/hg002_crosstech_venn_50plus.png', sep = '')
Cairo(600, 600, file = outfile, type = "png", bg = "white")
grid.newpage()
venn.plot <- draw.triple.venn(area1           = pb_only + pb_ont + pb_hifi + allthree,
                              area2           = ont_only + pb_ont + ont_hifi + allthree,
                              area3           = hifi_only + pb_hifi + ont_hifi + allthree,
                              n12             = pb_ont + allthree,
                              n23             = ont_hifi + allthree,
                              n13             = pb_hifi + allthree,
                              n123            = allthree,
                              category        = c('CLR', 'ONT', 'HiFi'),
                              fill            = c('#DD3D2D', '#FEDA8B', '#C2D4EF'),
                              #cat.col         = c('#DD3D2D', '#FEDA8B', '#C2D4EF'),#c('dodgerblue', 'seagreen3', 'orchid3'),
                              cex             = 2,
                              cat.cex         = 2,
                              euler           = TRUE,
                              scaled          = FALSE
)
dev.off()

#ggsave(venn.plot, outfile, device = 'cairo_png')

fn <- paste(projectroot, '/figures/figure3/hg002_crosstech.merged.tsv', sep = '')
outfile <- paste(projectroot, '/figures/figure3/long_hg002_crosstech_suppvecs.png', sep = '')
df <- read.table(fn, sep = "\t", header = TRUE)
nrow(df)
df <- df %>% filter(abs(SVLEN) > 50 | SVTYPE == 'TRA')
nrow(df)
caller <- 'Jasmine CrossTech'
suppvec_hist(df, "Jasmine CrossTech", outfile, FALSE)
outfile <- paste(projectroot, '/figures/figure3/long_hg002_crosstech_suppvecs_spec_prec.png', sep = '')
suppvec_hist(df, "Jasmine CrossTech", outfile, TRUE)

# Hifi-unique variants
outfile <- paste(projectroot, '/figures/figure3/long_hg002_crosstech_hifi.png', sep = '')
tech_specific(df %>% filter(SUPP_VEC == 100), "HiFi-Only", outfile, TRUE)
outfile <- paste(projectroot, '/figures/figure3/long_hg002_crosstech_hifi_indels.png', sep = '')
tech_specific_length(df %>% filter(SUPP_VEC == 100), "HiFi Only", outfile, TRUE)

# CLR-unique variants
outfile <- paste(projectroot, '/figures/figure3/long_hg002_crosstech_clr.png', sep = '')
tech_specific(df %>% filter(SUPP_VEC == 010), "CLR-Only", outfile, TRUE)
outfile <- paste(projectroot, '/figures/figure3/long_hg002_crosstech_clr_indels.png', sep = '')
tech_specific_length(df %>% filter(SUPP_VEC == 010), "CLR Only", outfile, TRUE)

# ONT-unique variants
outfile <- paste(projectroot, '/figures/figure3/long_hg002_crosstech_ont.png', sep = '')
tech_specific(df %>% filter(SUPP_VEC == 001), "ONT-Only", outfile, TRUE)
outfile <- paste(projectroot, '/figures/figure3/long_hg002_crosstech_ont_indels.png', sep = '')
tech_specific_length(df %>% filter(SUPP_VEC == 001), "ONT Only", outfile, TRUE)

# Concordant variants
outfile <- paste(projectroot, '/figures/figure3/long_hg002_crosstech_agree.png', sep = '')
tech_specific(df %>% filter(SUPP_VEC == 111), "Concordant", outfile, TRUE)
outfile <- paste(projectroot, '/figures/figure3/long_hg002_crosstech_agree_indels.png', sep = '')
tech_specific_length(df %>% filter(SUPP_VEC == 111), "Concordant", outfile, TRUE)


df <- df %>% filter(IS_SPECIFIC == 1 & IS_PRECISE == 1)
pb_only <- nrow(df %>% filter(SUPP_VEC == 010))
ont_only <- nrow(df %>% filter(SUPP_VEC == 001))
hifi_only <- nrow(df %>% filter(SUPP_VEC == 100))
pb_ont <- nrow(df %>% filter(SUPP_VEC == 011))
pb_hifi <- nrow(df %>% filter(SUPP_VEC == 110))
ont_hifi <- nrow(df %>% filter(SUPP_VEC == 101))
allthree <- nrow(df %>% filter(SUPP_VEC == 111))

outfile <- paste(projectroot, '/figures/figure3/long_hg002_crosstech_venn.png', sep = '')
Cairo(600, 600, file = outfile, type = "png", bg = "white")
grid.newpage()
venn.plot <- draw.triple.venn(area1           = pb_only + pb_ont + pb_hifi + allthree,
                              area2           = ont_only + pb_ont + ont_hifi + allthree,
                              area3           = hifi_only + pb_hifi + ont_hifi + allthree,
                              n12             = pb_ont + allthree,
                              n23             = ont_hifi + allthree,
                              n13             = pb_hifi + allthree,
                              n123            = allthree,
                              category        = c('CLR', 'ONT', 'HiFi'),
                              fill            = c('dodgerblue', 'seagreen3', 'orchid3'),
                              cat.col         = c('dodgerblue', 'seagreen3', 'orchid3'),
                              cex             = 2,
                              cat.cex         = 2,
                              euler           = TRUE,
                              scaled          = FALSE
)
dev.off()

#ggsave(venn.plot, outfile, device = 'cairo_png')
