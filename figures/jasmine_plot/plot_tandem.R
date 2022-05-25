library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)
library(reshape2)
library(gridExtra)
library(rlist)
#library(ggfortify)
library(tidyverse)

fn <- "/home/mkirsche/jasmine_data/revisions/hg002_trio_tandem.tsv"
df = read.table(fn, comment.char = "#", sep = "\t", header = T, stringsAsFactors=FALSE, colClasses = 'character')
nrow(df)
nrow(df %>% filter(SUPP_VEC=="100"))
df$SVLEN <- as.numeric(df$SVLEN)
lengthfiltered <- df %>% filter(SVLEN >= 50 | SVLEN <= -50 | SVLEN == 0 | SVTYPE=="TRA")
nrow(lengthfiltered)
nrow(lengthfiltered %>% filter(SUPP_VEC=="100"))


fn <- "/home/mkirsche/jasmine_data/revisions/tr_overlap_nonzero.tsv"
df = read.table(fn, comment.char = "#", sep = "\t", header = T, stringsAsFactors=FALSE, colClasses = 'character')
nrow(df)
nrow(df %>% filter(SUPP_VEC=="100"))
df$SVLEN <- as.numeric(df$SVLEN)
lengthfiltered <- df %>% filter(SVLEN >= 50 | SVLEN <= -50 | SVLEN == 0 | SVTYPE=="TRA")
nrow(lengthfiltered)
nrow(lengthfiltered %>% filter(SUPP_VEC=="100"))



fn <- "/home/mkirsche/jasmine_data/revisions/hiconf/test_all_500.vcf.tsv"#"/home/mkirsche/jasmine_data/revisions/hiconf/hg002_trio_hiconf.tsv"
df = read.table(fn, comment.char = "#", sep = "\t", header = T, stringsAsFactors=FALSE, colClasses = 'character')
nrow(df)
nrow(df %>% filter(SUPP_VEC=="100"))
df$SVLEN <- as.numeric(df$SVLEN)
lengthfiltered <- df %>% filter(SVLEN >= 50 | SVLEN <= -50 | SVLEN == 0 | SVTYPE=="TRA")
nrow(lengthfiltered)
nrow(lengthfiltered %>% filter(SUPP_VEC=="100"))

fn <- "/home/mkirsche/jasmine_data/revisions/hiconf/jasmine_cohort_hiconf.tsv"
df = read.table(fn, comment.char = "#", sep = "\t", header = T, stringsAsFactors=FALSE, colClasses = 'character')
nrow(df)
nrow(df %>% filter(SUPP_VEC=="100"))
df$SVLEN <- as.numeric(df$SVLEN)
lengthfiltered <- df %>% filter(SVLEN >= 50 | SVLEN <= -50 | SVLEN == 0 | SVTYPE=="TRA")
nrow(lengthfiltered)
nrow(lengthfiltered %>% filter(SUPP_VEC=="100"))
