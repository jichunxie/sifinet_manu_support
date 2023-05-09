library(Matrix)
source("../../../code/run_M3Drop.R")

set.seed(1)
data1 <- readRDS("../matrix.rds")
data1 <- as.matrix(data1)
out <- run_M3Drop(data1)

saveRDS(out, "M3Drop.rds")
