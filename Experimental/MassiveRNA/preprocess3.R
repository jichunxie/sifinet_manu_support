library(Matrix)
#setwd("~/Desktop")
cell_idx <- list()
p <- 27998
fs <- data.frame(matrix(0, p, 10))
totalcell <- 0
for (i in 0:9){
  data <- readMM(paste("matrix", i, "_h.mtx", sep = ""))
  rs <- rowSums(data > 0)
  cell_idx[[i+1]] <- which(rs >= (p * 0.02))
  data <- data[cell_idx[[i+1]], ]
  totalcell <- totalcell + length(cell_idx[[i+1]])
  fs[, i+1] <- colSums(data > 0)
}

feature_idx <- which(rowSums(fs) >= (totalcell * 0.02))
saveRDS(cell_idx, "cell_idx.rds")
saveRDS(feature_idx, "feature_idx.rds")
cell_idx <- readRDS("cell_idx.rds")
feature_idx <- readRDS("feature_idx.rds")

for (i in 0:9){
  data <- readMM(paste("matrix", i, "_h.mtx", sep = ""))
  data <- data[cell_idx[[i+1]], ]
  data <- data[, feature_idx]
  writeMM(data, paste("matrix", i, "_r.mtx", sep = ""))
}

rm(list = ls())
