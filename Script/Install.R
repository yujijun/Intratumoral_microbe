# BASIC INFO ---------------------------
##
## Script name: Install.R 
##
## Purpose of script: A script which is used for installing R package in intratumoral microbe projects
##
## Author: JijunYu
##
## Date Created: 2022-09-20
## Update Date:
## 
## Copyright (c) Jijun Yu, 2022
## Email: jijunyuedu@outlook.com
##
## Notes:
##   
#--------------Main -----------------------

## R 
rm(list = ls())
## 设置镜像，国外用户请忽略
options(BioC_mirror="http://mirrors.tuna.tsinghua.edu.cn/bioconductor/")
options("repos" = c(CRAN="http://mirrors.tuna.tsinghua.edu.cn/CRAN/"))

## R packages
bioPackages <-c( 
  "phyloseq",
  "corrplot",
  "ggplot2",
  "data.table",
  "tidyverse"
)


## install BiocManager
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

## install devtools
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")

## install packages
lapply( bioPackages, 
        function( bioPackage ){
          if(!bioPackage %in% rownames(installed.packages())){
            CRANpackages <- available.packages()
            if(bioPackage %in% rownames(CRANpackages)){
              install.packages( bioPackage)
            }else{
              BiocManager::install(bioPackage, suppressUpdates=F, ask=F)
            }
          }
        })