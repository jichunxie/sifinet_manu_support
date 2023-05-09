source("../../../code/run_MAST.R")
  
cluster_setting <- c("true", "seurat", "louvain", "cidr")

for (i in 31:40){
  set.seed(i)
  data1 <- readRDS(paste("../data/", i, "_matrix.rds", sep = ""))
  
  for (j in cluster_setting){
    cluster <- readRDS(paste("../data/", i, "_cluster_", j, ".rds", sep = ""))
    out <- run_MAST(data1, cluster)
    saveRDS(out, paste(i, "_MAST_", j, ".rds", sep = ""))
  }
}

rm(list = ls())
