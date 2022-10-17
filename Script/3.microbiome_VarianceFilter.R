microbiome_VarianceFilter <- function(data, filtopt, filtPerct){
  otu <- otu_table(data)%>%as.data.frame()
  tax <- tax_table(data)%>%as.data.frame()                   #接着上面过滤后的数据继续过滤
  rmn_feat <- nrow(otu)
  filter.val <- nm <- NULL
  if(filtPerct == 0) {
    remain <- rep(TRUE, rmn_feat)
  }else {
    if (filtopt == "iqr") {
      filter.val <- apply(otu, 1, IQR, na.rm = T)            #计算每一行的IQR值
      nm <- "IQR"
    }
    else if (filtopt == "sd") {                         
      filter.val <- apply(otu, 1, sd, na.rm = T)             #计算每一行的标准差
      nm <- "standard deviation"
    }
    else if (filtopt == "cov") {
      sds <- apply(otu, 1, sd, na.rm = T)
      mns <- apply(otu, 1, mean, na.rm = T)
      filter.val <- abs(sds/mns)
      nm <- "Coeffecient of variation"                       #变异系数=标准差/平均值
    }
    rk <- rank(-filter.val, ties.method = "random")          #rank()是求秩的函数，它的返回值是这个向量中对应元素的“排名”。　　"random" 是相同元素随机编排次序，避免了“先到先得”，“权重”优于“前后顺序”的机制增大了随机的程度。
    var.num <- nrow(otu)
    remain <- rk < var.num * (1 - filtPerct)                 #将菌种数量*（1-0.1），如果rk中的样本排名小于这个菌种cutoff，则返回TRUE
  }
  otu <- otu[remain, ]%>% as.matrix()%>%otu_table(taxa_are_rows=T)   #根据上述条件过滤菌种
  tax <- tax[remain, ]%>% as.matrix()%>%tax_table()
  new_phyloseq <- phyloseq(otu,tax,data@sam_data)
  write_rds(new_phyloseq,file = "3.filt.idtrans.orig.Rds")
}