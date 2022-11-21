#β-多样性（样本间多样性）
rm(list = ls())
setwd("~/microbiome_project/yujijun/Intratumoral_microbe/Script")
library(dplyr)
library(phyloseq)
library(amplicon)
#读入过滤后的phyloseq
mat = readRDS("4.filt.idtrans.orig.Rds")
mat <- rarefy_even_depth(mat, replace = FALSE, rngseed = T)

metadata <- sample_data(mat)%>%as.data.frame()
metadata <- data.frame(rownames(metadata),metadata$project,row.names = rownames(metadata))
colnames(metadata) <- c("SampleID","project")

# 读取特征ASV表，多样性分析输入抽平标准化的表
otutab <- otu_table(mat)%>%as.data.frame()
# 基于特征表、原数据、距离类型、分组列名、是否添加置信椭圆，是否添加样本标签
p = beta_cpcoa(otutab, metadata, dis="bray", groupID="project", ellipse=T, label=T)
ggsave("p1.cpcoa_otutab.pdf", p, width=10, height=10, units="mm")
ggsave("p1.cpcoa_otutab.png", p, width=10, height=10, units="mm", dpi=300)

#接下来对每个癌症取平均值，然后再做CPCoA分析
#按照T/N合并样本
COAD_tumor <- subset_samples(mat, project =="COAD")%>%subset_samples(Sample =="PT")%>%
  otu_table()%>%as.data.frame()%>%apply(1,mean)
COAD_normal <- subset_samples(mat, project =="COAD")%>%subset_samples(Sample =="STN")%>%
  otu_table()%>%as.data.frame()%>%apply(1,mean)
ESCA_tumor <- subset_samples(mat, project =="ESCA")%>%subset_samples(Sample =="PT")%>%
  otu_table()%>%as.data.frame()%>%apply(1,mean)
ESCA_normal <- subset_samples(mat, project =="ESCA")%>%subset_samples(Sample =="STN")%>%
  otu_table()%>%as.data.frame()%>%apply(1,mean)
HNSC_tumor <- subset_samples(mat, project =="HNSC")%>%subset_samples(Sample =="PT")%>%
  otu_table()%>%as.data.frame()%>%apply(1,mean)
HNSC_normal <- subset_samples(mat, project =="HNSC")%>%subset_samples(Sample =="STN")%>%
  otu_table()%>%as.data.frame()%>%apply(1,mean)
READ_tumor <- subset_samples(mat, project =="READ")%>%subset_samples(Sample =="PT")%>%
  otu_table()%>%as.data.frame()%>%apply(1,mean)
READ_normal <- subset_samples(mat, project =="READ")%>%subset_samples(Sample =="STN")%>%
  otu_table()%>%as.data.frame()%>%apply(1,mean)
STAD_tumor <- subset_samples(mat, project =="STAD")%>%subset_samples(Sample =="PT")%>%
  otu_table()%>%as.data.frame()%>%apply(1,mean)
STAD_normal <- subset_samples(mat, project =="STAD")%>%subset_samples(Sample =="STN")%>%
  otu_table()%>%as.data.frame()%>%apply(1,mean)
#合并不同癌症数据
df <- cbind(COAD_tumor,COAD_normal,
            ESCA_tumor,ESCA_normal,HNSC_tumor,HNSC_normal,READ_tumor,READ_normal,STAD_tumor,STAD_normal)

metadata <- data.frame(colnames(df),row.names = colnames(df))
colnames(metadata) <- "Var1"
# 基于特征表、原数据、距离类型、分组列名、是否添加置信椭圆，是否添加样本标签
p = beta_cpcoa(df2, metadata, dis="bray", groupID="Var1", ellipse=T, label=T)
ggsave("p1.cpcoa_otutab.pdf", p, width=10, height=10, units="mm")
ggsave("p1.cpcoa_otutab.png", p, width=10, height=10, units="mm", dpi=300)