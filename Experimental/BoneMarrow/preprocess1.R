#setwd("~/Desktop/SiFINeT/Result/Experimental/BoneMarrow")

gsvares <- readRDS("gsva_res.rds")
# table(colSums(gsvares > 0))
color <- ifelse(gsvares[1,] > 0, 1, ifelse(gsvares[2,] > 0, 2, 3))
plot(gsvares[1,], gsvares[2,], col = color)

library(Matrix)
data_mat <- readRDS("matrix.rds")
cellinfo <- readRDS("cellinfo.rds")
genename <- readRDS("genename.rds")

id <- which(color == 2)
cellinfo_1 <- cellinfo[id, ]
data_mat_1 <- data_mat[, id]
id3 <- which(rowSums(data_mat_1 > 0) > 0.02 * ncol(data_mat_1))
data_mat_1 <- data_mat_1[id3,]
genenames_1 <- genename[id3]

saveRDS(data_mat_1, "matrix_1.rds")
saveRDS(genenames_1, "genename_1.rds")
saveRDS(cellinfo_1, "cellinfo_1.rds")
rm(list = ls())
