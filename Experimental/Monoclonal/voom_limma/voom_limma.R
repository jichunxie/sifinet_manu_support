library(Matrix)
source("../../../code/run_voom_limma.R")

set.seed(1)
data1 <- readRDS("../matrix.rds")
data1 <- as.matrix(data1)
cluster <- readRDS("../cluster.rds")
out <- run_voom_limma(data1, cluster)

saveRDS(out, "voom_limma.rds")
