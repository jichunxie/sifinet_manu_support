library(Seurat)
library(Matrix)

set.seed(1)

data <- readRDS("PreprocessedData/CD8_rna_matrix.rds")
gsva_mat <- readRDS("gsva_res_rna.rds")
gsva_mat2 <- rbind(gsva_mat, rep(0, ncol(gsva_mat)))
SiFINeT_group <- apply(gsva_mat2[,], 2, which.max)

genename <- readRDS("PreprocessedData/CD8_genename.rds")

rownames(data) <- genename
colnames(data) <- 1:ncol(data)

coldata <- data.frame(SiFINeT_group = SiFINeT_group)
rownames(coldata) <- 1:ncol(data)

so <- CreateSeuratObject(data, meta.data = coldata)
so <- SCTransform(so, verbose=TRUE, return.only.var.genes = FALSE)
Idents(so) <- "SiFINeT_group"
test_pair <- rbind(1:4, rep(5, 4))

results <- list()
for(i in 1:ncol(test_pair)) {
  markers <- FindMarkers(so, ident.1 = test_pair[2, i], ident.2 = test_pair[1, i], 
                         return.thresh = 1.1, min.pct = 0, logfc.threshold = 0)
  markers$gene <- rownames(markers)
  results[[i]] <- markers
}


saveRDS(results, "rna_test_pair.rds")
