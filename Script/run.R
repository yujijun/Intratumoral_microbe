rm(list = ls())
setwd("~/Rdata/yujijun/Intratumoral_microbe/Script")
library("phyloseq")
library("data.table")
library("tidyverse")
library("corrplot")
library("ggplot2")
library("reshape2")
#load phyloseq object
rt <- readRDS("~/Rdata/yujijun/Duke_TCMA_20220917/WGS/sample/physeq.bacteria.unambiguous.decontam.tissue.sample.reads.species.Rds")
tax <- read.delim("~/Rdata/yujijun/Duke_TCMA_20220917/metadata/taxonomy.txt")

source("1.microbiome_taxIDtrans.R")
source("2.microbiome_AbundanceFilter.R")
source("3.microbiome_VarianceFilter.R")
source("4.microbiome_sampleFilter.R")
source("5.microbiome_Normalization.R")
#未完成
# source("5.microbiome_PlotTaxaAundanceBar.R")
# source("6.HeatTreePlot.R")
# source("7.microbiome_PlotAlphaData.R")

#运行函数
df <- idtrans(rt,tax)
df2 <- microbiome_AbundanceFilter(df,"prevalence",1,0.1)
df3 <- microbiome_VarianceFilter(df2,"iqr",0.1)
df4 <- microbiome_sampleFilter(df3,"mean",0)
df5 <- PerformNormalization(df4,"none","colsum","none")


#library(metacoder)
#library(viridis)
# PrepareHeatTreePlot(mbSet, "Sample", "Species", "dbgr", "dft", "PT_vs_STN", 0.05, "heat_tree_0","png",300)
# PlotAlphaData("filt","alpha_diver_0","Chao1","Sample","OTU", "default", "png")
