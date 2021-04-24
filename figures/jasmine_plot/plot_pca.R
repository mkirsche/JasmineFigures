library(tidyverse)

#fn <- "/home/mkirsche/plink.eigenvec"
#spofn <- "/home/mkirsche/nygc_pca_sp.png"
#ofn <- "/home/mkirsche/nygc_pca.png"

#fn <- "/home/mkirsche/jasmine_pca.eigenvec"
#spofn <- "/home/mkirsche/jasmine_paragraph_pca_sp.png"
#ofn <- "/home/mkirsche/jasmine_paragraph_pca.png"

#fn <- "/home/mkirsche/merged.eigenvec"
#spofn <- "/home/mkirsche/merged_sp.png"
#ofn <- "/home/mkirsche/merged.png"

fn <- "/home/mkirsche/jasmine_data/figures/figure5/plink.eigenvec"
spofn <- "/home/mkirsche/jasmine_data/figures/figure5/jasmine_paragraph_pca_sp.png"
ofn <- "/home/mkirsche/jasmine_data/figures/figure5/jasmine_paragraph_pca.png"


df <- read.table(fn, header = FALSE)[, c(1, 3:5)]
colnames(df) <- c("sample_id", "pc1", "pc2", "pc3")
popfn <- "/home/mkirsche/eclipse-workspace/ParagraphGenotypes/all.idx"
popdf <- read.table(popfn, sep = " ", header = F)
colnames(popdf)
popdf
names = c("V2")
l <- as.vector(df$sample_id)
l
colordf <- data.frame(matrix(l, nrow=length(l), byrow=TRUE))
colnames(colordf) <- names
colnames(colordf)
nrow(colordf)
colordf
merged <- merge.data.frame(colordf, popdf, all.x=TRUE)
merged
merged$V3 <- ifelse(is.na(merged$V3), "LR", as.character(merged$V3))
merged$SP <- merged$V3
merged$SP <- ifelse(merged$SP == "LR", "LR", as.character(merged$SP))
merged$SP <- ifelse(merged$SP == "CHB", "EAS", as.character(merged$SP))
merged$SP <- ifelse(merged$SP == "JPT", "EAS", as.character(merged$SP))
merged$SP <- ifelse(merged$SP == "CHS", "EAS", as.character(merged$SP))
merged$SP <- ifelse(merged$SP == "CDX", "EAS", as.character(merged$SP))
merged$SP <- ifelse(merged$SP == "KHV", "EAS", as.character(merged$SP))
merged$SP <- ifelse(merged$SP == "CEU", "EUR", as.character(merged$SP))
merged$SP <- ifelse(merged$SP == "TSI", "EUR", as.character(merged$SP))
merged$SP <- ifelse(merged$SP == "FIN", "EUR", as.character(merged$SP))
merged$SP <- ifelse(merged$SP == "GBR", "EUR", as.character(merged$SP))
merged$SP <- ifelse(merged$SP == "IBS", "EUR", as.character(merged$SP))
merged$SP <- ifelse(merged$SP == "YRI", "AFR", as.character(merged$SP))
merged$SP <- ifelse(merged$SP == "LWK", "AFR", as.character(merged$SP))
merged$SP <- ifelse(merged$SP == "GWD", "AFR", as.character(merged$SP))
merged$SP <- ifelse(merged$SP == "MSL", "AFR", as.character(merged$SP))
merged$SP <- ifelse(merged$SP == "ESN", "AFR", as.character(merged$SP))
merged$SP <- ifelse(merged$SP == "ASW", "AFR", as.character(merged$SP))
merged$SP <- ifelse(merged$SP == "ACB", "AFR", as.character(merged$SP))
merged$SP <- ifelse(merged$SP == "MXL", "AMR", as.character(merged$SP))
merged$SP <- ifelse(merged$SP == "PUR", "AMR", as.character(merged$SP))
merged$SP <- ifelse(merged$SP == "CLM", "AMR", as.character(merged$SP))
merged$SP <- ifelse(merged$SP == "PEL", "AMR", as.character(merged$SP))
merged$SP <- ifelse(merged$SP == "GIH", "SAS", as.character(merged$SP))
merged$SP <- ifelse(merged$SP == "PJL", "SAS", as.character(merged$SP))
merged$SP <- ifelse(merged$SP == "BEB", "SAS", as.character(merged$SP))
merged$SP <- ifelse(merged$SP == "STU", "SAS", as.character(merged$SP))
merged$SP <- ifelse(merged$SP == "ITU", "SAS", as.character(merged$SP))
merged <- merged %>% map_df(rev)
Population <- merged$V3
SuperPopulation <- merged$SP
nrow(df)
length(SuperPopulation)
df <- df %>% map_df(rev)
df$V3 <- merged$V3
df$V2 <- merged$V2
lrpoints <- subset(df, V3 == "LR")
nrow(lrpoints)
lrpoints
ggplot(data = df, aes(x = pc1, y = pc2)) +
  geom_point(aes(color = SuperPopulation)) +
  theme_bw(base_size=32) + 
  guides(color = guide_legend(override.aes = list(size=10))) +
  xlab('PC1') +
  ylab('PC2') +
  ggtitle("Genotyped SVs in 1KGP Samples (Super Population)") +
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        axis.text.x = element_text(size = 12, angle = 30, margin = margin(t = 15)),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16), 
        legend.title = element_blank()
  ) +
  geom_text(data=lrpoints, aes(pc1, pc2, label = sample_id), size = 2)
ggsave(spofn, width = 12, height = 8)

ggplot(data = df, aes(x = pc1, y = pc2)) +
  geom_point(aes(color = Population)) +
  theme_bw(base_size=32) + 
  guides(color = guide_legend(override.aes = list(size=10))) +
  xlab('PC1') +
  ylab('PC2') +
  ggtitle("Genotyped SVs in 1KGP Samples (Population)") +
theme(plot.title = element_text(size = 18, hjust = 0.5),
      axis.text.x = element_text(size = 12, angle = 30, margin = margin(t = 15)),
      axis.text.y = element_text(size = 12),
      axis.title.x = element_text(size = 16),
      axis.title.y = element_text(size = 16), 
      legend.title = element_blank()
) +
  geom_text(data=lrpoints, aes(pc1, pc2, label = sample_id), size = 2)
ggsave(ofn, width = 12, height = 8)








ggplot(data = df, aes(x = pc1, y = pc2)) +
  geom_point(aes(color = SuperPopulation)) +
  theme_bw(base_size=32) + 
  guides(color = guide_legend(override.aes = list(size=10))) +
  xlab('PC1') +
  ylab('PC2') + xlim(0, .0125) + ylim(-.05, .0125) +
  ggtitle("NYGC Small Variant Calls (SuperPopulation)") +
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        axis.text.x = element_text(size = 12, angle = 30, margin = margin(t = 15)),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16), 
        legend.title = element_blank()
  ) +
  geom_text(data=lrpoints, aes(pc1, pc2, label = sample_id))
ggsave(paste(ofn, "_2.png", sep=""), width = 12, height = 8)


fn <- "/home/mkirsche/t2t_smallvars.eigenvec"
spofn <- "/home/mkirsche/t2t_sp.png"
ofn <- "/home/mkirsche/t2t.png"

df <- read.table(fn, header = FALSE)[, c(1, 3:5)]
colnames(df) <- c("sample_id", "pc1", "pc2", "pc3")
popfn <- "/home/mkirsche/eclipse-workspace/ParagraphGenotypes/all.idx"
popdf <- read.table(popfn, sep = " ", header = F)
colnames(popdf)
popdf
names = c("V1")
l <- as.vector(df$sample_id)
l
colordf <- data.frame(matrix(l, nrow=length(l), byrow=TRUE))
colnames(colordf) <- names
colnames(colordf)
nrow(colordf)
colordf
merged <- merge.data.frame(colordf, popdf, all.x=TRUE)
merged
merged$SP <- merged$V3
merged$SP <- ifelse(merged$SP == "CHB", "EAS", as.character(merged$SP))
merged$SP <- ifelse(merged$SP == "JPT", "EAS", as.character(merged$SP))
merged$SP <- ifelse(merged$SP == "CHS", "EAS", as.character(merged$SP))
merged$SP <- ifelse(merged$SP == "CDX", "EAS", as.character(merged$SP))
merged$SP <- ifelse(merged$SP == "KHV", "EAS", as.character(merged$SP))
merged$SP <- ifelse(merged$SP == "CEU", "EUR", as.character(merged$SP))
merged$SP <- ifelse(merged$SP == "TSI", "EUR", as.character(merged$SP))
merged$SP <- ifelse(merged$SP == "FIN", "EUR", as.character(merged$SP))
merged$SP <- ifelse(merged$SP == "GBR", "EUR", as.character(merged$SP))
merged$SP <- ifelse(merged$SP == "IBS", "EUR", as.character(merged$SP))
merged$SP <- ifelse(merged$SP == "YRI", "AFR", as.character(merged$SP))
merged$SP <- ifelse(merged$SP == "LWK", "AFR", as.character(merged$SP))
merged$SP <- ifelse(merged$SP == "GWD", "AFR", as.character(merged$SP))
merged$SP <- ifelse(merged$SP == "MSL", "AFR", as.character(merged$SP))
merged$SP <- ifelse(merged$SP == "ESN", "AFR", as.character(merged$SP))
merged$SP <- ifelse(merged$SP == "ASW", "AFR", as.character(merged$SP))
merged$SP <- ifelse(merged$SP == "ACB", "AFR", as.character(merged$SP))
merged$SP <- ifelse(merged$SP == "MXL", "AMR", as.character(merged$SP))
merged$SP <- ifelse(merged$SP == "PUR", "AMR", as.character(merged$SP))
merged$SP <- ifelse(merged$SP == "CLM", "AMR", as.character(merged$SP))
merged$SP <- ifelse(merged$SP == "PEL", "AMR", as.character(merged$SP))
merged$SP <- ifelse(merged$SP == "GIH", "SAS", as.character(merged$SP))
merged$SP <- ifelse(merged$SP == "PJL", "SAS", as.character(merged$SP))
merged$SP <- ifelse(merged$SP == "BEB", "SAS", as.character(merged$SP))
merged$SP <- ifelse(merged$SP == "STU", "SAS", as.character(merged$SP))
merged$SP <- ifelse(merged$SP == "ITU", "SAS", as.character(merged$SP))
merged <- merged %>% map_df(rev)
Population <- merged$V3
SuperPopulation <- merged$SP
nrow(df)
length(SuperPopulation)
df <- df %>% map_df(rev)
df$V3 <- merged$V3
df$V2 <- merged$V2
ggplot(data = df, aes(x = pc1, y = pc2)) +
  geom_point(aes(color = SuperPopulation)) +
  theme_bw(base_size=32) + 
  guides(color = guide_legend(override.aes = list(size=10))) +
  xlab('PC1') +
  ylab('PC2') +
  ggtitle("T2T Small Variant Calls (Super Population)") +
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        axis.text.x = element_text(size = 12, angle = 30, margin = margin(t = 15)),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16), 
        legend.title = element_blank()
  )
ggsave(spofn, width = 12, height = 8)

ggplot(data = df, aes(x = pc1, y = pc2)) +
  geom_point(aes(color = Population)) +
  theme_bw(base_size=32) + 
  guides(color = guide_legend(override.aes = list(size=10))) +
  xlab('PC1') +
  ylab('PC2') +
  ggtitle("T2T Small Variant Calls (Population)") +
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        axis.text.x = element_text(size = 12, angle = 30, margin = margin(t = 15)),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16), 
        legend.title = element_blank()
  )
ggsave(ofn, width = 12, height = 8)
