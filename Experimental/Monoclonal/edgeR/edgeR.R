library(Matrix)
source("../../../code/run_edgeR.R")

set.seed(1)
data1 <- readRDS("../matrix.rds")
data1 <- as.matrix(data1)
cluster <- readRDS("../cluster.rds")
out <- run_edgeR(data1, cluster)

saveRDS(out, "edgeR.rds")
