#setwd("~/Desktop/SiFINeT/Result/Experimental/IPF")

# Seurat
library(Seurat)

data <- readRDS("2_matrix.rds")
group <- readRDS("2_cellinfo.rds")
genename <- readRDS("2_genename.rds")

group <- 1 * (group == "KRT5-/KRT17+")
rownames(data) <- 1:nrow(data)
colnames(data) <- 1:ncol(data)
coldata <- data.frame(group = group)
rownames(coldata) <- 1:ncol(data)

so <- CreateSeuratObject(data, meta.data = coldata)
so <- SCTransform(so, #vars.to.regress = c('sample'), 
                  verbose=TRUE, return.only.var.genes = FALSE)
Idents(so) <- "group"
out <- FindMarkers(so, ident.1 = 1, return.thresh = 1.1, min.pct = 0, logfc.threshold = 0)

out$gene <- genename[as.numeric(rownames(out))]

saveRDS(out, "2_test.rds")
