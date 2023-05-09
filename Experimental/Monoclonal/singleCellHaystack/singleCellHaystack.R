library(Matrix)
source("../../../code/run_singleCellHaystack.R")

set.seed(1)
data1 <- readRDS("../matrix.rds")
data1 <- as.matrix(data1)
out <- run_singleCellHaystack(data1)

saveRDS(out, "singleCellHaystack.rds")
