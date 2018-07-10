setwd("F:/1000 Brassica napus genome/玄表皮毛/xuan_five_gene_block_20180626")
library(LDheatmap)
library(genetics)
BnaC07_snp <- read.table("BnaC07_snp.txt", header = TRUE, check.names = FALSE)
BnaC07_dist <- read.table("BnaC07_dist.txt", header = FALSE)
BnaC07_dist <- as.numeric(BnaC07_dist)
for (i in 1:ncol(BnaC07_snp)){
  BnaC07_snp[,i] <- as.genotype(BnaC07_snp[,i])
}
MyHeatmap <- LDheatmap(BnaC07_snp, genetic.distances = BnaC07_dist,color = grey.colors(20))
rgb.palette <- colorRampPalette(rev(c("yellow","red")), space="rgb")
pdf("BnaC07_block.pdf")
LDheatmap(MyHeatmap,color = rgb.palette(100), flip = TRUE)
dev.off()
