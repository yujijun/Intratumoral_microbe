rm(list = ls())
library("phyloseq")
library("data.table")
library("tidyverse")
library("corrplot")
library("ggplot2")
library("reshape2")
library("dplyr")
setwd("~/Rdata/yujijun/Intratumoral_microbe")
#########################################
#load phyloseq object
rt <- readRDS("~/Rdata/yujijun/Duke_TCMA_20220917/WXS/case/physeq.bacteria.unambiguous.decontam.tissue.case.reads.species.Rds")

##create new sample_table
sample_info <- rt@sam_data
sample_info <- data.frame(rt@sam_data@.Data,row.names = rt@sam_data@row.names)
colnames(sample_info) <- rt@sam_data@names
sample_info
write.table(sample_info,"sample_info.tissue.case.reads.txt",sep = "\t",quote = F,row.names = F)


#read the taxonomy file
tax <- read.delim("~/Rdata/yujijun/Duke_TCMA_20220917/metadata/taxonomy.txt")
tax <- tax[,c("tax_id","type","name")]
#界(Kingdom)、门(Phylum)、纲(Class)、目(Order)、科(Family)、属(Genus)、种(Species)
#shape the data
new_tax <- reshape2::dcast(tax,tax_id~type)
#extract tax_table
tax_mat <- as.matrix(rt@tax_table@.Data)
tax_mat <- cbind(tax_id=rownames(tax_mat),tax_mat)
tax_mat <- as.data.frame(tax_mat)
#change the tax_id in tax_table into bacteria
tax_mat$superkingdom <- new_tax$superkingdom[match(tax_mat$superkingdom, new_tax$tax_id)]
tax_mat$phylum <- new_tax$phylum[match(tax_mat$phylum, new_tax$tax_id)]
tax_mat$class <- new_tax$class[match(tax_mat$class, new_tax$tax_id)]
tax_mat$order <- new_tax$order[match(tax_mat$order, new_tax$tax_id)]
tax_mat$family <- new_tax$family[match(tax_mat$family, new_tax$tax_id)]
tax_mat$genus <- new_tax$genus[match(tax_mat$genus, new_tax$tax_id)]
tax_mat$species <- new_tax$species[match(tax_mat$species, new_tax$tax_id)]


#读入otu_table
otu_table <- rt@otu_table@.Data%>%as.data.frame()
#tax_id转菌名
otu_table <- cbind(ID=rownames(otu_table),otu_table)
otu_table$ID <- tax_mat$phylum[match(tax_mat$tax_id, otu_table$ID)] 
###利用循环将其中的重复门数据进行求和
data<-aggregate(`TCGA-2H-A9GF` ~ ID,data=otu_table,sum)
for (i in colnames(otu_table)[3:length(colnames(otu_table))]){
  data1<-aggregate(otu_table[,i]~ID,data=otu_table,sum)
  colnames(data1)[2]<-i  
  data<-merge(data,data1,by="ID")
}
df1 <- data
rownames(df1)=df1$ID
df1 <- df1[,-1]

#otu_table过滤没有表达的菌种和样本  
df2 <- df1[rowMeans(df1)>0,]
df3 <- df2[,colMeans(df2)>0]
output <- cbind(ID=rownames(df3),df3)
write.table(output,"phylum.tissue.case.reads.txt",sep = "\t",quote = F)
#########################################
#read immune subtype file
imm <- read.table("~/Rdata/yujijun/TCGA_panimmune_20220917/Scores_160_Signatures.tsv",sep = "\t",header = T,check.names = F)

Cibersort <- filter(imm,Source=="Cibersort") %>% mutate(Source= NULL) %>% 
  data.frame(row.names = 1,check.names = F) %>% as.matrix() %>% t() %>% as.data.frame()
Wolf <- filter(imm,Source=="Wolf") %>% mutate(Source= NULL) %>% 
  data.frame(row.names = 1,check.names = F) %>% as.matrix() %>% t() %>% as.data.frame()
Attractors <- filter(imm,Source=="Attractors") %>% mutate(Source= NULL) %>% 
  data.frame(row.names = 1,check.names = F) %>% as.matrix() %>% t() %>% as.data.frame()
ICR <- filter(imm,Source=="ICR") %>% mutate(Source= NULL) %>% 
  data.frame(row.names = 1,check.names = F) %>% as.matrix() %>% t() %>% as.data.frame()
c7atoms <- filter(imm,Source=="c7atoms") %>% mutate(Source= NULL) %>% 
  data.frame(row.names = 1,check.names = F) %>% as.matrix() %>% t() %>% as.data.frame()
Bindea <- filter(imm,Source=="Bindea") %>% mutate(Source= NULL) %>% 
  data.frame(row.names = 1,check.names = F) %>% as.matrix() %>% t() %>% as.data.frame()



#merge otu_table and imm_table for correlation analysis
imm_cor <- function(immset=Cibersort,dat=df3,name="Cibersort"){
  # immset=Cibersort
  # dat=df3
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
imm_cor(Cibersort,df3,"Cibersort")
imm_cor(Wolf,df3,"Wolf")
imm_cor(ICR,df3,"ICR")
imm_cor(c7atoms,df3,"c7atoms")
imm_cor(Bindea,df3,"Bindea")
