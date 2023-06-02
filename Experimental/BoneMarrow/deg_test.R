library(Seurat)
library(Matrix)

#setwd("~/Desktop/SiFINeT/Result/Experimental/BoneMarrow")
set.seed(1)

data <- readRDS("matrix.rds")
gsva_mat <- readRDS("gsva_res.rds")
SiFINeT_group <- 1 * ((gsva_mat[1,] <= 0) & (gsva_mat[2,] <= 0))

genename <- readRDS("genename.rds")

rownames(data) <- genename
colnames(data) <- 1:ncol(data)

coldata <- data.frame(SiFINeT_group = SiFINeT_group)
rownames(coldata) <- 1:ncol(data)

so <- CreateSeuratObject(data, meta.data = coldata)

Idents(so) <- "SiFINeT_group"
out1 <- FindMarkers(so, ident.1 = 1, return.thresh = 1.1, min.pct = 0, logfc.threshold = 0)
out1$gene <- rownames(out1)

saveRDS(out1, "layer1_test.rds")



set.seed(1)

data <- readRDS("matrix_1.rds")
gsva_mat <- readRDS("gsva_res_1.rds")
SiFINeT_group <- 1 * ((gsva_mat[1,] <= 0) & 
                        (gsva_mat[2,] <= 0) & 
                        (gsva_mat[3,] <= 0))

genename <- readRDS("genename_1.rds")

rownames(data) <- genename
colnames(data) <- 1:ncol(data)

coldata <- data.frame(SiFINeT_group = SiFINeT_group)
rownames(coldata) <- 1:ncol(data)

so <- CreateSeuratObject(data, meta.data = coldata)

Idents(so) <- "SiFINeT_group"
out2 <- FindMarkers(so, ident.1 = 1, return.thresh = 1.1, min.pct = 0, logfc.threshold = 0)
out2$gene <- rownames(out2)

saveRDS(out2, "layer2_test.rds")