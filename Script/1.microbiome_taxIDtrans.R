idtrans <- function(data,tax){
  data <- rt
  #read the taxonomy file
  tax <- tax[,c("tax_id","type","name")]
  #界(Kingdom)、门(Phylum)、纲(Class)、目(Order)、科(Family)、属(Genus)、种(Species)
  #shape the data
  new_tax <- reshape2::dcast(tax,tax_id~type)
  #extract tax_table
  tax_mat <- as.matrix(data@tax_table@.Data)
  tax_mat <- cbind(tax_id=rownames(tax_mat),tax_mat)
  tax_mat <- as.data.frame(tax_mat)
  #change the tax_id in tax_table into bacteria
  tax_mat$kingdom <- new_tax$superkingdom[match(tax_mat$superkingdom, new_tax$tax_id)]
  tax_mat$phylum <- new_tax$phylum[match(tax_mat$phylum, new_tax$tax_id)]
  tax_mat$class <- new_tax$class[match(tax_mat$class, new_tax$tax_id)]
  tax_mat$order <- new_tax$order[match(tax_mat$order, new_tax$tax_id)]
  tax_mat$family <- new_tax$family[match(tax_mat$family, new_tax$tax_id)]
  tax_mat$genus <- new_tax$genus[match(tax_mat$genus, new_tax$tax_id)]
  tax_mat$species <- new_tax$species[match(tax_mat$species, new_tax$tax_id)]
  #extract needed columns
  tax_mat <- tax_mat[,c("phylum","class","order","family","genus","species")]   #纲目科属种
  TAX <- tax_table(tax_mat %>% as.matrix())
  #重新创建新的phyloseq对象
  new_rt <- phyloseq(data@otu_table,TAX,data@sam_data)
  write_rds(new_rt,file = "1.idtrans.orig.Rds")
}