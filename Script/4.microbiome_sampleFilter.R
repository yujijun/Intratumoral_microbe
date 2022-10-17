microbiome_sampleFilter <- function (data, filt.opt, count) 
{
  otu <- otu_table(data)%>%as.data.frame()
  if (filt.opt == "mean") {
      filter.val <- apply(otu, 2, mean, na.rm = T)       #如果函数传的第二个参数是mean,则平均值大于参数3的菌种才被保留
      kept.inx <- filter.val > count
  }
  otu <- otu[,kept.inx]%>% as.matrix()%>%otu_table(taxa_are_rows=T)
  new_phyloseq <- phyloseq(otu,data@tax_table,data@sam_data)
  write_rds(new_phyloseq,file = "4.filt.idtrans.orig.Rds")
}