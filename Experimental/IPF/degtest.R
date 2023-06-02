library(Seurat)
library(Matrix)

#setwd("~/Desktop/SiFINeT/Result/Experimental/IPF")
set.seed(1)

data <- readRDS("1_matrix.rds")
gsva_mat <- readRDS("gsva_res.rds")
SiFINeT_group <- 1 * ((gsva_mat[2,] >= quantile(gsva_mat[2,], 0.6)) & 
                        (gsva_mat[3,] <= quantile(gsva_mat[3,], 0.4)) & 
                        (gsva_mat[1,] <= quantile(gsva_mat[1,], 0.4)))

gsvares <- readRDS("gsva_res_SenMayo.rds")
group <- (gsvares[1,] >= quantile(gsvares, 1 - sum(SiFINeT_group) / ncol(gsvares))) * 1
genename <- readRDS("1_genename.rds")

library(mclust)
print(adjustedRandIndex(SiFINeT_group, group))


library(stringr)
rownames(data) <- genename
colnames(data) <- 1:ncol(data)

coldata <- data.frame(group = group, SiFINeT_group = SiFINeT_group)
rownames(coldata) <- 1:ncol(data)

so <- CreateSeuratObject(data, meta.data = coldata)
so <- SCTransform(so, verbose=TRUE, return.only.var.genes = FALSE)

Idents(so) <- "group"
out1 <- FindMarkers(so, ident.1 = 1, return.thresh = 1.1, min.pct = 0, logfc.threshold = 0)
out1$gene <- rownames(out1)

Idents(so) <- "SiFINeT_group"
out2 <- FindMarkers(so, ident.1 = 1, return.thresh = 1.1, min.pct = 0, logfc.threshold = 0)
out2$gene <- rownames(out2)


saveRDS(out1, "senmayo_test.rds")
saveRDS(out2, "SiFINeT_test.rds")







data <- readRDS("2_matrix.rds")
group <- readRDS("2_cellinfo.rds")
genename <- readRDS("2_genename.rds")

group <- 1 * (group == "KRT5-/KRT17+")
rownames(data) <- genename
colnames(data) <- 1:ncol(data)
coldata <- data.frame(group = group)
rownames(coldata) <- 1:ncol(data)

so <- CreateSeuratObject(data, meta.data = coldata)
so <- SCTransform(so, verbose=TRUE, return.only.var.genes = FALSE)
Idents(so) <- "group"
out <- FindMarkers(so, ident.1 = 1, return.thresh = 1.1, min.pct = 0, logfc.threshold = 0)

out$gene <- rownames(out)

saveRDS(out, "2_test.rds")
