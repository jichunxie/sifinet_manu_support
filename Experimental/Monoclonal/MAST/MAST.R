library(Matrix)
source("../../../code/run_MAST.R")

set.seed(1)
data1 <- readRDS("../matrix.rds")
data1 <- as.matrix(data1)
cluster <- readRDS("../cluster.rds")
out <- run_MAST(data1, cluster)

saveRDS(out, "MAST.rds")
