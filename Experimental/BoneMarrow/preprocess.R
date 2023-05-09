#setwd("~/Desktop/SiFINeT/Result/Experimental/BoneMarrow")

designdata <- read.delim("GSE72857_experimental_design.txt", sep = "\t", header = T, skip = 19)

data_mat <- read.delim("GSE72857_umitab.txt")
genenames <- rownames(data_mat)

cells <- designdata$Well_ID[designdata$Batch_desc == "Unsorted myeloid"]

designdata <- designdata[match(cells, designdata$Well_ID), ]
rownames(designdata) <- designdata$Well_ID
data_mat <- data_mat[, match(cells, colnames(data_mat))]


library(Seurat)
library(sctransform)

so <- CreateSeuratObject(data_mat, meta.data = designdata)
so <- SCTransform(so, vars.to.regress = c('Amp_batch_ID'), 
                  verbose=TRUE, return.only.var.genes = FALSE)

data <- so@assays$SCT@data

id2 <- which(colSums(data_mat) > 500)
data <- data[, id2]
designdata <- designdata[id2, ]

id3 <- which(rowSums(data > 0) > 0.02 * ncol(data))
data <- data[id3,]
genenames <- rownames(data)

saveRDS(data, "matrix.rds")
saveRDS(genenames, "genename.rds")
saveRDS(designdata, "cellinfo.rds")
rm(list = ls())
