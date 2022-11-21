rm(list = ls())
setwd("~/microbiome_project/yujijun/Intratumoral_microbe/Script")
#安装包
# library(devtools)
# install_github("jokergoo/ComplexHeatmap")
#绘制复杂热图
library("ComplexHeatmap")
library(dplyr)
library(phyloseq)
mat = readRDS("4.filt.idtrans.orig.Rds")
otu <- otu_table(mat)
tax <-  tax_table(mat)
#按照指定的Top数量进行筛选
Top=50
level="species"
for (j in 1:dim(tax)[1]) {
  if (row.names(tax)[j] %in% names(sort(rowSums(otu), decreasing = TRUE)[1:Top])) {
    tax[j,level] =tax[j,level]
  } else {
    tax[j,level]= "Other"
  }
}
top_species <- tax%>%as.data.frame()%>%filter(species!="Other")

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
df2 <- df[rownames(top_species),]
#pheatmap(df2,scale = "row")
#注释列
# annotation_col = sample_data(dd)%>%as.data.frame()
# row.names(annotation_col) <- colnames(top_species)
#注释行
rownames(df2) <- top_species[,"species"]
top_species2 <- mutate(top_species,species=NULL) #删除种列
annotation_row = top_species2
row.names(annotation_row) <- rownames(top_species2)
#每画一次热图就会随机设置一种颜色，所以为了保持一致，最好设置随机种子
set.seed(123) 
B <- pheatmap(df2, scale = "row",cluster_rows = F,cluster_cols = F,border_color=NA,
              col = colorRampPalette(c("navy","white","firebrick3"))(1000),
              # annotation_col = annotation_col, #列注释信息
              annotation_row = annotation_row,#行注释信息
              row_split = annotation_row$phylum,#行截断（按照pathway，不像之前随机）
              # column_split = annotation_col$group,#列截断
              # annotation_names_row = F,#不显示行注释信息
              # annotation_names_col = F ,#不显示列注释信息
              column_title = NULL,#不显示列标题
              row_title = NULL)#不显示行标题
pdf("figure1a.pdf",width = 15,height = 10)
print(B)
dev.off()


