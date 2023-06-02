library(Seurat)
library(Matrix)

#setwd("~/Desktop/SiFINeT/Result/Experimental/BoneMarrow")
set.seed(1)

data <- readRDS("PreprocessedData/CD8_rna_matrix.rds")
gsva_mat <- readRDS("gsva_res_rna.rds")
SiFINeT_group <- 1 * ((gsva_mat[1,] <= 0) & (gsva_mat[2,] <= 0) &
                        (gsva_mat[3,] <= 0) & (gsva_mat[4,] <= 0))

genename <- readRDS("PreprocessedData/CD8_genename.rds")

rownames(data) <- genename
colnames(data) <- 1:ncol(data)

coldata <- data.frame(SiFINeT_group = SiFINeT_group)
rownames(coldata) <- 1:ncol(data)

so <- CreateSeuratObject(data, meta.data = coldata)
so <- SCTransform(so, verbose=TRUE, return.only.var.genes = FALSE)

Idents(so) <- "SiFINeT_group"
out1 <- FindMarkers(so, ident.1 = 1, return.thresh = 1.1, min.pct = 0, logfc.threshold = 0)
out1$gene <- rownames(out1)

saveRDS(out1, "rna_test.rds")



set.seed(1)

data <- readRDS("PreprocessedData/CD8_atac_matrix.rds")
gsva_mat <- readRDS("gsva_res_atac.rds")
SiFINeT_group <- 1 * ((gsva_mat[1,] <= 0) & 
                        (gsva_mat[2,] <= 0))

sum(SiFINeT_group) == 1 # 1