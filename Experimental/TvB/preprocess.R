#setwd("~/Desktop/SiFINeT/Result/Experimental/TvB/")
library(Matrix)

cell_annotation <- readRDS("../TCell/OriginalData/multi_annotation.rds")

cellidx <- which(cell_annotation %in% c("CD8 Naive", "CD4 Naive", 
                                        "Intermediate B", "Memory B", "Naive B"))
data <- readRDS("../TCell/OriginalData/multi_rna_data.rds")
data <- data[, cellidx]
genenames <- rownames(data)
rpc <- colSums(data)
idx <- which(rpc <= 1.8 * mean(rpc))
data <- data[, intersect(idx, idx)]
idx <- which(colSums(data > 0) > 0.02 * nrow(data))
data <- data[, idx]
idx <- which(rowSums(data > 0) > 0.02 * ncol(data))
data <- data[idx, ]
genenames <- genenames[idx]

saveRDS(data, "matrix.rds")
saveRDS(genenames, "genenames.rds")
rm(list = ls())
