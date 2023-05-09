library(biomaRt)
#setwd("~/Desktop/SiFINeT/Result/Experimental/Monoclonal/cellcycle")

gene_name <- readRDS("../all_gene_name.rds")
mart <- useMart('ensembl', dataset="mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org")
annotation <- getBM(attributes=c("uniprot_gn_id", "uniprot_gn_symbol"), 
                    filters="uniprot_gn_symbol", values=gene_name, 
                    mart=mart)

write.table(annotation, "uniid-symbol.txt", quote = F, row.names = F, sep = "\t")
write.table(annotation[,1], "uniid.txt", quote = F, row.names = F, sep = "\t", col.names = F)
