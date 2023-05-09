source("../../../code/run_DESeq2.R")

cluster_setting <- c("true", "seurat", "louvain", "cidr")

for (i in 51:60){
  set.seed(i)
  data1 <- readRDS(paste("../data/", i, "_matrix.rds", sep = ""))
  for (j in cluster_setting){
    cluster <- readRDS(paste("../data/", i, "_cluster_", j, ".rds", sep = ""))
    out <- run_DESeq2(data1, cluster)

    saveRDS(out, paste(i, "_DESeq2_", j, ".rds", sep = ""))
  }
}

rm(list = ls())
