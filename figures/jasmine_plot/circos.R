library(circlize)

library(RCircos)
library(dplyr)
library(ggplot2)
library(Cairo)

fn <-  "/home/mkirsche/jasmine_data/figures/figure5/population_final_specprec.merged.tsv"#"/home/mkirsche/jasmine_data/figures/figure5/all4_pop_jasmine_specprec_noTRA.tsv"

df <- read.table(fn, header = TRUE, sep = "\t")
colnames(df)
unique(df$CHROM)
goodlist <- paste('chr', c('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', 
                           '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', "X", "Y"), sep = "")

goodlist
df <- df %>% filter(CHROM %in% goodlist)
df$CHR <- droplevels(df$CHROM)

unique(df$CHROM)
binsize <- 2000000
df$STARTBIN <- floor(df$POS / binsize) * binsize
groups <- df %>% group_by(CHROM, STARTBIN, SVTYPE) %>% tally()
groups$CHROM <- droplevels(groups$CHROM)

unique(groups$CHROM)
groups$logn <- log2(groups$n)
outfile <- "/home/mkirsche/jasmine_data/figures/figure5/flatcircos.png"
ggplot(groups %>% filter(CHROM == "chr1"), aes(x = STARTBIN, y = logn, fill = SVTYPE)) + geom_bar( stat = "identity")
ggsave(outfile, width = 16, height = 6)

insgroups <- groups %>% filter(SVTYPE == "INS")
delgroups <- groups %>% filter(SVTYPE == "DEL")
invgroups <- groups %>% filter(SVTYPE == "INV")
dupgroups <- groups %>% filter(SVTYPE == "DUP")

tmp <- data.frame(dupgroups)
for (x in goodlist) {
  if((x %in% dupgroups$CHROM) == FALSE)
  {
    de<-data.frame(paste(x, "", sep = ""),0, "DUP", 0)
    names(de)<-c("CHROM","STARTBIN", "SVTYPE", "logn")
    tmp <- rbind(tmp, de)
  }
}
dupgroups <- tmp

insgroups
dupgroups

outfile <- "/home/mkirsche/jasmine_data/figures/figure5/circos.svg"
#png(outfile, width = 4096, height = 4096)
Cairo(file=outfile, 
      type="svg",
      units="in", 
      width=10, 
      height=10, 
      pointsize=12, 
      dpi=72)
circos.par(cell.padding=c(0, .1, 0, .1), track.height = .16)
colorpalette <- c(DEL = '#CC6677',DUP = '#332288',INS = '#44AA99',TRA = '#88CCEE',INV = '#DDaa33')

circos.initializeWithIdeogram()
circos.track(insgroups$CHROM, x = insgroups$STARTBIN, y = insgroups$logn, panel.fun = function(x, y) {
  circos.barplot(y, x, col = '#44AA99', bar_width = binsize, lwd = 0.00001, lty = 0)
})
circos.track(delgroups$CHROM, x = delgroups$STARTBIN, y = delgroups$logn, panel.fun = function(x, y) {
  circos.barplot(y, x, col = '#CC6677', bar_width = binsize, lwd = 0.00001, lty = 0)
})
circos.track(invgroups$CHROM, x = invgroups$STARTBIN, y = invgroups$logn, panel.fun = function(x, y) {
  circos.barplot(y, x, col = '#ddaa33', bar_width = binsize, lwd = 0.00001, lty = 0)
})
circos.track(dupgroups$CHROM, x = dupgroups$STARTBIN, y = dupgroups$logn, panel.fun = function(x, y) {
  circos.barplot(y, x, col = '#332288', bar_width = binsize, lwd = 0.00001, lty = 0)
})

#circos.text(c(20, 20, 20, 20), c(10, 10, 10, 10), c("Insertions", "Deletions", "Inversions", "Duplications"), c("chr15", "chr15", "chr15", "chr15"), c(3, 4, 5, 6))
dev.off()

outfile <- "/home/mkirsche/jasmine_data/figures/figure5/circoslegend.svg"
#png(outfile, width = 4096, height = 4096)
Cairo(file=outfile, 
      type="svg",
      units="in", 
      width=4, 
      height=6, 
      pointsize=12, 
      dpi=72)
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("topleft", 
       legend = c("Insertions", "Deletions", "Inversions", "Duplications"), 
       col = c('#44AA99', '#CC6677', '#DDAA33', '#332288'), 
       pch = 16,
       bty = "n", 
       pt.cex = 3, 
       cex = 1.5, 
       text.col = "black")
dev.off()



## Now make it again but with only SVs
fn <-  "/home/mkirsche/jasmine_data/figures/figure5/population_final_specprec.merged.tsv"#"/home/mkirsche/jasmine_data/figures/figure5/all4_pop_jasmine_specprec_noTRA.tsv"
df <- read.table(fn, header = TRUE, sep = "\t")
df <- df %>% filter(SVLEN >= 50 | SVLEN <= -50 | SVLEN == 0 | SVTYPE == "TRA")
colnames(df)
unique(df$CHROM)
goodlist <- paste('chr', c('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', 
                           '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', "X", "Y"), sep = "")

goodlist
df <- df %>% filter(CHROM %in% goodlist)
df$CHR <- droplevels(df$CHROM)

unique(df$CHROM)
binsize <- 2000000
df$STARTBIN <- floor(df$POS / binsize) * binsize
groups <- df %>% group_by(CHROM, STARTBIN, SVTYPE) %>% tally()
groups$CHROM <- droplevels(groups$CHROM)

unique(groups$CHROM)
groups$logn <- log2(groups$n)
outfile <- "/home/mkirsche/jasmine_data/figures/figure5/flatcircos_50plus.png"
ggplot(groups %>% filter(CHROM == "chr1"), aes(x = STARTBIN, y = logn, fill = SVTYPE)) + geom_bar( stat = "identity")
ggsave(outfile, width = 16, height = 6)

insgroups <- groups %>% filter(SVTYPE == "INS")
delgroups <- groups %>% filter(SVTYPE == "DEL")
invgroups <- groups %>% filter(SVTYPE == "INV")
dupgroups <- groups %>% filter(SVTYPE == "DUP")

tmp <- data.frame(dupgroups)
for (x in goodlist) {
  if((x %in% dupgroups$CHROM) == FALSE)
  {
    de<-data.frame(paste(x, "", sep = ""),0, "DUP", 0)
    names(de)<-c("CHROM","STARTBIN", "SVTYPE", "logn")
    tmp <- rbind(tmp, de)
  }
}
dupgroups <- tmp

insgroups
dupgroups

outfile <- "/home/mkirsche/jasmine_data/figures/figure5/circos_50plus.svg"
Cairo(file=outfile, 
      type="svg",
      units="in", 
      width=10, 
      height=10, 
      pointsize=12, 
      dpi=72)
circos.par(cell.padding=c(0, .1, 0, .1), track.height = .16)
colorpalette <- c(DEL = '#CC6677',DUP = '#332288',INS = '#44AA99',TRA = '#88CCEE',INV = '#DDaa33')

circos.initializeWithIdeogram()
circos.track(insgroups$CHROM, x = insgroups$STARTBIN, y = insgroups$logn, panel.fun = function(x, y) {
  circos.barplot(y, x, col = '#44AA99', bar_width = binsize, lwd = 0.00001, lty = 0)
})
circos.track(delgroups$CHROM, x = delgroups$STARTBIN, y = delgroups$logn, panel.fun = function(x, y) {
  circos.barplot(y, x, col = '#CC6677', bar_width = binsize, lwd = 0.00001, lty = 0)
})
circos.track(invgroups$CHROM, x = invgroups$STARTBIN, y = invgroups$logn, panel.fun = function(x, y) {
  circos.barplot(y, x, col = '#ddaa33', bar_width = binsize, lwd = 0.00001, lty = 0)
})
circos.track(dupgroups$CHROM, x = dupgroups$STARTBIN, y = dupgroups$logn, panel.fun = function(x, y) {
  circos.barplot(y, x, col = '#332288', bar_width = binsize, lwd = 0.00001, lty = 0)
})

dev.off()


