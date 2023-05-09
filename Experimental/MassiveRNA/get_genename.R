library(restfulSE)
my10x = se1.3M()
genenames <- my10x@elementMetadata@listData$gene_name
saveRDS(genenames, "genenames.rds")
