setwd("F:/1000 Brassica napus genome/flower_time_GEM_association_analysis/flc_analysis/gene_flower/snp_analysis/LD_analysis")
library(LDheatmap)
library(genetics)
BnaA02_snp <- read.table("BnaA02_snp.txt", header = TRUE, check.names = FALSE)
BnaA02_dist <- read.table("BnaA02_dist.txt", header = FALSE)
BnaA02_dist <- as.numeric(BnaA02_dist)
for (i in 1:ncol(BnaA02_snp)){
  BnaA02_snp[,1] <- as.genotype(BnaA02_snp[,1])
}
rgb.palette <- colorRampPalette(rev(c("white","red")), space="rgb")
LDheatmap(MyHeatmap,color = rgb.palette(100), flip = TRUE)
