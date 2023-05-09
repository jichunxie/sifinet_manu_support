library(Matrix)
source("../../../code/run_DESeq2.R")

set.seed(1)
data1 <- readRDS("../matrix.rds")
data1 <- as.matrix(data1)
cluster <- readRDS("../cluster.rds")
out <- run_DESeq2(data1, cluster)
    
saveRDS(out, "DESeq2.rds")

