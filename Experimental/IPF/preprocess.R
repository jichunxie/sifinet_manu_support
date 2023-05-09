library(Matrix)

#setwd("~/Desktop/SiFINeT/Result/Experimental/IPF/")
metadata <- read.csv("GSE135893_IPF_metadata.csv")
idx1 <- which((metadata$celltype == "Ciliated") & (metadata$Diagnosis == "IPF"))
idx2 <- which(metadata$celltype %in% c("KRT5-/KRT17+", "Proliferating Epithelial Cells"))

barcode <- read.table("GSE135893_barcodes.tsv")
idx <- match(metadata$X, barcode[,1])
# all.equal(barcode[idx,1], metadata$X)
genename <- read.table("GSE135893_genes.tsv")
data <- readMM("GSE135893_matrix.mtx")
data <- data[, idx]

data1 <- data[, idx1]
data1 <- data1[, colSums(data1 > 0) > 0.05 * nrow(data1)]
i <- rowSums(data1 > 0) > 0.05 * ncol(data1)
data1 <- data1[i,]
saveRDS(data1, "1_matrix.rds")
saveRDS(genename[i,1], "1_genename.rds")

data2 <- data[, idx2]
cellinfo <- metadata$celltype[idx2]
j <- colSums(data2 > 0) > 0.05 * nrow(data2)
data2 <- data2[, j]
cellinfo <- cellinfo[j]
i <- rowSums(data2 > 0) > 0.05 * ncol(data2)
data2 <- data2[i,]
saveRDS(data2, "2_matrix.rds")
saveRDS(genename[i,1], "2_genename.rds")
saveRDS(cellinfo, "2_cellinfo.rds")
rm(list = ls())
