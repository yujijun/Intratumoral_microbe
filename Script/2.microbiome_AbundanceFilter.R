microbiome_AbundanceFilter <- function (data, filt.opt, count, smpl.perc) 
{
  otu <- otu_table(data)%>%as.data.frame()
  tax <- tax_table(data)%>%as.data.frame()
  rmn_feat <- nrow(otu)              #data变量应该就是otu_table,这里统计行数rmn_feat，也就是一共检测到多少个菌种
  if(count == 0) {
    kept.inx <- rep(TRUE, rmn_feat)   #如果count过滤的cutoff为0，就把TRUE重复行数的数目
  }else {
    if (filt.opt == "prevalence") {           #如果函数传的第二个参数是prevalence，则大于参数3的样本量大于参数4时菌种才被保留
      rmn_feat <- nrow(otu)
      minLen <- smpl.perc * ncol(otu)        #smpl.perc=0.2表示，20%样本数量
      kept.inx <- apply(otu, MARGIN = 1, function(x) {
        sum(x >= count) >= minLen             #按行，对大于count=4的样本数量进行统计，如果该菌种样本数量大于20%的总样本量,则返回TRUE,否则这一行菌种返回为FALSE
      })
    }
    else if (filt.opt == "mean") {
      filter.val <- apply(otu, 1, mean, na.rm = T)       #如果函数传的第二个参数是mean,则平均值大于参数3的菌种才被保留
      kept.inx <- filter.val >= count
    }
    else if (filt.opt == "median") {
      filter.val <- apply(otu, 1, median, na.rm = T)     #如果函数传的第二个参数是median,则中位数大于参数3的菌种才被保留
      kept.inx <- filter.val >= count
    }
  }
  otu <- otu[kept.inx, ]%>% as.matrix()%>%otu_table(taxa_are_rows=T)
  tax <- tax[kept.inx, ]%>% as.matrix()%>%tax_table()
  new_phyloseq <- phyloseq(otu,tax,data@sam_data)
  write_rds(new_phyloseq,file = "2.filt.idtrans.orig.Rds")
}

#microbiome_AbundanceFilter(rt,"prevalence",1,0.1)
