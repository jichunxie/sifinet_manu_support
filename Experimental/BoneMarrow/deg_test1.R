library(Seurat)
library(Matrix)

#setwd("~/Desktop/SiFINeT/Result/Experimental/BoneMarrow")
set.seed(1)

data <- readRDS("matrix.rds")
gsva_mat <- readRDS("gsva_res.rds")
gsva_mat2 <- rbind(gsva_mat, rep(0, ncol(gsva_mat)))
SiFINeT_group <- apply(gsva_mat2[,], 2, which.max)

genename <- readRDS("genename.rds")

rownames(data) <- genename
colnames(data) <- 1:ncol(data)

coldata <- data.frame(SiFINeT_group = SiFINeT_group)
rownames(coldata) <- 1:ncol(data)

so <- CreateSeuratObject(data, meta.data = coldata)

Idents(so) <- "SiFINeT_group"
test_pair <- rbind(1:2, rep(3, 2))

results <- list()
for(i in 1:ncol(test_pair)) {
  markers <- FindMarkers(so, ident.1 = test_pair[2, i], ident.2 = test_pair[1, i], 
                         return.thresh = 1.1, min.pct = 0, logfc.threshold = 0)
  markers$gene <- rownames(markers)
  results[[i]] <- markers
}

saveRDS(out1, "layer1_test_pair.rds")



set.seed(1)

data <- readRDS("matrix_1.rds")
gsva_mat <- readRDS("gsva_res_1.rds")
gsva_mat2 <- rbind(gsva_mat, rep(0, ncol(gsva_mat)))
SiFINeT_group <- apply(gsva_mat2[,], 2, which.max)

genename <- readRDS("genename_1.rds")

rownames(data) <- genename
colnames(data) <- 1:ncol(data)

coldata <- data.frame(SiFINeT_group = SiFINeT_group)
rownames(coldata) <- 1:ncol(data)

so <- CreateSeuratObject(data, meta.data = coldata)

Idents(so) <- "SiFINeT_group"
test_pair <- rbind(1:3, rep(4, 3))

results <- list()
for(i in 1:ncol(test_pair)) {
  markers <- FindMarkers(so, ident.1 = test_pair[2, i], ident.2 = test_pair[1, i], 
                         return.thresh = 1.1, min.pct = 0, logfc.threshold = 0)
  markers$gene <- rownames(markers)
  results[[i]] <- markers
}

saveRDS(out2, "layer2_test_pair.rds")
