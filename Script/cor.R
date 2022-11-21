rm(list = ls())
library("ggpubr")
library("phyloseq")
library("data.table")
library("tidyverse")
library("corrplot")
library("ggplot2")
library("reshape2")
library("dplyr")
setwd("~/Rdata/yujijun/Intratumoral_microbe")
#########################################
#文件都需要整理成行名为样本名（12位），数据框内都为非全零的数值型连续变量
#TMB <- read.table("~/Rdata/yujijun/TCGA_panimmune_20220917/mutation-load_updated.txt",sep = "\t",header = T,check.names = F)
# neoantigen <- read.table("~/Rdata/yujijun/TCGA_panimmune_20220917/
#                          TCGA_PCA.mc3.v0.2.8.CONTROLLED.filtered.sample_neoantigens_10062017.tsv",sep = "\t",header = T,check.names = F)
# neoantigen <- data.frame(neoantigen,row.names = 1)
# neoantigen <- neoantigen[,-ncol(neoantigen)]
# pMHCSNV <- read.table("~/Rdata/yujijun/TCGA_panimmune_20220917/
#                          TCGA_pMHC_SNV_sampleSummary_MC3_v0.2.8.CONTROLLED_170404.tsv",sep = "\t",header = T,check.names = F)
#igh <- read.table("~/Rdata/yujijun/TCGA_panimmune_20220917/tcga.pancan.igh.div.txt",sep = "\t",header = T,check.names = F)
# TCGASubtype <- read.table("~/Rdata/yujijun/TCGA_panimmune_20220917/TCGASubtype.20170308.tsv",sep = "\t",header = T,check.names = F)
# viral <- read.table("~/Rdata/yujijun/TCGA_panimmune_20220917/viral.tsv",sep = "\t",header = T,check.names = F)
# viral <- viral[,c(2,4,6:17)]%>%filter(Study=="COAD")
# viral <- viral[!duplicated(viral$ParticipantBarcode),]
# rownames(viral) <- viral[,1]
# viral <- viral[,-1]
# viral <- viral[,-1]
# viral <- viral[,colMeans(viral)>0]
TMESubtype <- read.table("~/Rdata/yujijun/TMEsubtype/TME.txt",sep = "\t",header = T,check.names = F,row.names = 1)
# TMESubtype <- TMESubtype[,c("Purity","TMB")]
TMESubtype2 <- TMESubtype[,c("MFP","Immune Subtype","P_stage")]
colnames(TMESubtype2)[2] <- "Immune_Subtype"


####read phylum
phylum <- read.table("~/Rdata/yujijun/Intratumoral_microbe/phylum.tissue.case.reads.txt",sep = "\t",header = T,
                     check.names = F,row.names = 1)
sample_info <- read.delim("~/Rdata/yujijun/Intratumoral_microbe/sample_info.tissue.case.reads.txt",
                          sep = "\t",header = T,check.names = F,row.names = 1)
sample_info <- filter(sample_info,acronym=="COAD")



#连续变量的相关性分析
cor1 <- function(immset=Cibersort,dat=phylum,name="Cibersort",sample_info=sample_info){
  # immset=viral
  # dat=phylum
  #样本取交集
  same_samples <- intersect(substr(rownames(immset),1,12),colnames(dat))
  same_samples <- intersect(same_samples,rownames(filter(sample_info,acronym=="COAD")))
  immset <- immset[same_samples,]
  dat <- dat[,same_samples]
  
  #相关性分析
  outab1 <- data.frame()
  outab2 <- data.frame()
  for(i in rownames(dat)){
    r1 <- c()
    r2 <- c()
    for(j in colnames(immset)){
      corT <- cor.test(as.numeric(dat[i,]),as.numeric(immset[,j]),method = "spearman")
      r1 <- c(r1,corT$estimate)
      r2 <- c(r1,corT$p.value)
    }
    outab1 <- rbind(outab1,r1)
    outab2 <- rbind(outab2,r1)
  }
  colnames(outab1) <- colnames(immset)
  colnames(outab2) <- colnames(immset)
  rownames(outab1) <- rownames(dat)
  rownames(outab1) <- rownames(dat)
  outab1 <- na.omit(outab1)
  outab2 <- na.omit(outab2)
  #绘图
  # 加载依赖包
  library(pheatmap)
  a <- pheatmap(t(as.matrix(outab1)),cluster_row = TRUE,cluster_col = TRUE,scale="none",
                legend = TRUE,border=FALSE,
                cellwidth = 8, cellheight = 8,
                main = paste0("Correlation between bacteria phylum and ",name))
  ggsave(paste0(name,".pdf"),plot=a,width = 10,height = 10)
}
#run function
cor1(TMESubtype,phylum,"TMESubtype",sample_info)
cor1(viral,phylum,"viral",sample_info)
cor1(neoantigen,phylum,"neoantigen",sample_info)



#分类变量的差异丰度分析
cor2 <- function(immset=Cibersort,dat=phylum,sample_info=sample_info,i="MFP",j="Actinobacteria"){
  immset=TMESubtype2
  dat=phylum
  name="Cibersort"
  sample_info=sample_info
  #样本取交集
  same_samples <- intersect(substr(rownames(immset),1,12),colnames(dat))
  same_samples <- intersect(same_samples,rownames(sample_info))
  immset <- immset[same_samples,]
  dat2 <- dat[,same_samples]%>%as.matrix()%>%t()
  dat2 <- log2(dat2+0.001)
  #合并数据
  df <- cbind(immset,dat2)
  #绘图
  a <- ggviolin(df, x=i, y=j, fill = i, 
                width = 0.5,   #小提琴宽度
                add = "boxplot", add.params = list(fill="white",size=1))+ 
    stat_compare_means(comparisons = list(c(names(table(df[,"MFP"])[1]),names(table(df[,"MFP"])[2]))), label = "p.signif")
  ggsave(paste0(i,j,".pdf"),plot=a,width = 10,height = 10)
}
cor2(TMESubtype2,phylum,sample_info,"MFP","Actinobacteria")
cor2(TMESubtype2,phylum,sample_info,"MFP","Actinobacteria")
