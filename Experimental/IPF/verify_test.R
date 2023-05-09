library(Seurat)
library(ggplot2)
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
sio <- readRDS("so.rds")
gs <- sio@featureset
genename <- readRDS("1_genename.rds")

library(mclust)
print(adjustedRandIndex(SiFINeT_group, group))




library(Seurat)

#genename <- rownames(data)

library(stringr)
rownames(data) <- 1:nrow(data)
colnames(data) <- 1:ncol(data)

coldata <- data.frame(group = group, SiFINeT_group = SiFINeT_group)
rownames(coldata) <- 1:ncol(data)

so <- CreateSeuratObject(data, meta.data = coldata)
so <- SCTransform(so, #vars.to.regress = c('sample'), 
                  verbose=TRUE, return.only.var.genes = FALSE)


Idents(so) <- "group"
out1 <- FindMarkers(so, ident.1 = 1, return.thresh = 1.1, min.pct = 0, logfc.threshold = 0)
out1$gene <- genename[as.numeric(rownames(out1))]
#saveRDS(out, "2_seurat_sctransform_SenMayo_norm.rds")

Idents(so) <- "SiFINeT_group"
out2 <- FindMarkers(so, ident.1 = 1, return.thresh = 1.1, min.pct = 0, logfc.threshold = 0)
out2$gene <- genename[as.numeric(rownames(out2))]


saveRDS(out1, "senmayo_test.rds")
saveRDS(out2, "SiFINeT_test.rds")
