source("../../../code/run_limma.R")

cluster_setting <- c("true", "seurat", "louvain", "cidr")

for (i in 1:100){
  set.seed(i)
	data1 <- readRDS(paste("../data/", i, "_matrix.rds", sep = ""))
  
  for (j in cluster_setting){
    cluster <- readRDS(paste("../data/", i, "_cluster_", j, ".rds", sep = ""))
    out <- run_limma(data1, cluster)
    
    saveRDS(out, paste(i, "_limma_", j, ".rds", sep = ""))
  }
}

rm(list = ls())
