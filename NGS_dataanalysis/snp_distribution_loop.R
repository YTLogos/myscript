setwd("C:/Users/taoyan/Desktop/文章/snp_distribution/")
library(openxlsx)
library(tidyverse)
library(ggpubr)
info <- read.xlsx("other_gene_snp_info.xlsx",sheet = 10)
info$accession <- factor(info$accession)
info$type <- factor(info$type)
gene_info <- read.xlsx("candidate genes.xlsx",sheet = 3)
gene_id <- as.character(gene_info$Geneid)
gene_anno <- as.character(gene_info$info)
gene_vline_value <- as.numeric(gene_info$vline)
snp_info <- list()
for (i in 1:9){
  snp_info[[i]] <- read.xlsx("other_gene_snp_info.xlsx", sheet = i)
  snp_info[[i]]$POS <- factor(snp_info[[i]]$POS)
  snp_info_gather <- snp_info[[i]]%>%gather(key="accession",value="allele",-c(POS,CHROM))
  snp_info_gather$accession <- factor(snp_info_gather$accession)
  snp_info_final <- left_join(snp_info_gather,info, by=c("accession"="accession"))
  snp_info_final <- na.omit(snp_info_final)
  snp_info_final$allele <- factor(snp_info_final$allele)
  my_data_W <- snp_info_final%>%dplyr::filter(type=="W")
  my_data_SW <- snp_info_final%>%dplyr::filter(type=="SW")
  my_data_S <- snp_info_final%>%dplyr::filter(type=="S")
  W <- ggplot(my_data_W,aes(x=POS,y=accession,fill=allele))+
    geom_bar(stat = "identity")+
    scale_fill_manual(values = c("#EBE9CB","#CACAE7","#E93857"),
                      name="Allele code", breaks=c("0","1","2"))+
    theme(plot.background = element_blank())+
    theme(panel.background = element_blank())+
    theme(axis.text = element_blank())+
    theme(axis.line = element_blank())+
    theme(axis.title.x = element_blank())+
    theme(legend.position = "bottom")+
    theme(axis.ticks = element_blank())+
    theme(legend.text = element_text(size=10),
          legend.title = element_text(size=15))+
    theme(legend.title = element_text(size=12),
          legend.text = element_text(size=10))+
    geom_vline(xintercept = gene_vline_value[i], size=0.35,color="blue")+
    ylab("Winter\n(n=603)")+
    labs(title=gene_id[i], subtitle=gene_anno[i])+
    theme(axis.title.y = element_text(size=8,angle = 0,vjust = 0.5,hjust =0.5))
  SW <- ggplot(my_data_SW,aes(x=POS,y=accession,fill=allele))+
    geom_bar(stat = "identity")+
    scale_fill_manual(values = c("#EBE9CB","#CACAE7","#E83857"),
                      name="Allele code", breaks=c("0","1","2"))+
    theme(plot.background = element_blank())+
    theme(panel.background = element_blank())+
    theme(axis.text = element_blank())+
    theme(axis.line = element_blank())+
    theme(axis.title.x = element_blank())+
    theme(legend.position = "bottom")+
    theme(axis.ticks = element_blank())+
    theme(legend.text = element_text(size=10),
          legend.title = element_text(size=15))+
    theme(legend.title = element_text(size=12), legend.text = element_text(size=10))+
    geom_vline(xintercept = gene_vline_value[i], size=0.35,color="blue")+
    ylab("Semi-Winter\n(n=140)")+
    theme(axis.title.y = element_text(size=8,angle = 0,vjust = 0.5,hjust = 0.5))
  S <- ggplot(my_data_S,aes(x=POS,y=accession,fill=allele))+
    geom_bar(stat = "identity")+
    scale_fill_manual(values = c("#EBE9CB","#CACAE7","#E83857"),
                      name="Allele code", breaks=c("0","1","2"))+
    theme(plot.background = element_blank())+
    theme(panel.background = element_blank())+
    theme(axis.text = element_blank())+
    theme(axis.line = element_blank())+
    theme(axis.title.x = element_blank())+
    theme(legend.position = "bottom")+
    theme(axis.ticks = element_blank())+
    theme(legend.text = element_text(size=10),
          legend.title = element_text(size=15))+
    theme(legend.title = element_text(size=12), legend.text = element_text(size=10))+
    geom_vline(xintercept = gene_vline_value[i], size=0.35,color="blue")+
    labs(y="Spring\n(n=175)")+
    theme(axis.title.y = element_text(size=8,angle = 0,vjust = 0.5,hjust = 0.5))
  p <- ggarrange(W,SW,S,ncol = 1,nrow = 3,common.legend = TRUE,legend = "bottom",
            align = "v", heights = c(6.0,1.4,1.75))
  pdf(paste0(gene_id[i],"_snp_distribution.pdf"), height = 4,width = 5)
  print(p)
  dev.off()
}