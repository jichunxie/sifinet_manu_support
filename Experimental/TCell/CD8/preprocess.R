#setwd("~/Desktop/SiFINeT/Result/Experimental/TCell/")
setwd("../")
library(Matrix)

cell_annotation <- readRDS("OriginalData/multi_annotation.rds")
peak_annotation <- read.table("OriginalData/pbmc_granulocyte_sorted_10k_atac_peak_annotation.tsv", 
                              sep = "\t", header = T)

data <- readRDS("OriginalData/multi_atac_data.rds")
data <- as.matrix(data)

peak_annotation$peakname <- gsub("_","-", peak_annotation$peak) 
peak_annotation$peak_id <- 1:nrow(peak_annotation)
peak_annotation$distance <- as.numeric(peak_annotation$distance)

cellidx <- which(cell_annotation %in% c("CD8 Naive", "CD8 TEM_1", "CD8 TEM_2"))
data <- data[, cellidx]


# setting 2: within 5000 bps of TSS, group by genes
peak_annotation2 <- peak_annotation[(peak_annotation$distance <= 0) & 
                                      (peak_annotation$distance >= -5000) & 
                                      !is.na(peak_annotation$distance), ]
peak_annotation2 <- peak_annotation2[order(peak_annotation2$gene), ]

genenames2 <- unique(peak_annotation2$gene)

data2 <- matrix(0, length(genenames2), ncol(data))
peak_annotation2$gene_id <- 1:nrow(peak_annotation2) - 
  cumsum(duplicated(peak_annotation2$gene))


for (i in 1:nrow(data2)){
  if (i %% 100 == 0){
    print(i)
  }
  idx <- peak_annotation2$peak_id[peak_annotation2$gene_id == i]
  datatemp <- data[idx, ]
  if (length(idx) == 1){
    data2[i, ] <- datatemp
  } else {
    data2[i, ] <- colSums(datatemp)
  }
}
# rm(data)

data1 <- readRDS("OriginalData/multi_rna_data.rds")
data1 <- as.matrix(data1)
data1 <- data1[, cellidx]
geneidx <- match(genenames2, rownames(data1))
data1 <- data1[geneidx,]
genenames1 <- rownames(data1)

# atacseq no mito

rpc1 <- colSums(data1)
rpc2 <- colSums(data2)
idx1 <- which(rpc1 <= 1.8 * mean(rpc1))
idx2 <- which(rpc2 <= 1.8 * mean(rpc2))
data1 <- data1[, intersect(idx1, idx2)]
data2 <- data2[, intersect(idx1, idx2)]

idx1 <- which(colSums(data1 > 0) > 0.02 * nrow(data1))
idx2 <- which(colSums(data2 > 0) > 0.02 * nrow(data2))
data1 <- data1[, intersect(idx1, idx2)]
data2 <- data2[, intersect(idx1, idx2)]

idx1 <- which(rowSums(data1 > 0) > 0.02 * ncol(data1))
idx2 <- which(rowSums(data2 > 0) > 0.02 * ncol(data2))
data1 <- data1[intersect(idx1, idx2), ]
data2 <- data2[intersect(idx1, idx2), ]
genenames <- genenames1[intersect(idx1, idx2)]

data <- readRDS("OriginalData/multi_rna_data.rds")
idx <- match(colnames(data1), colnames(data))

rownames(data2) <- rownames(data1)
colnames(data2) <- colnames(data1)

saveRDS(data1, "PreprocessedData/CD8_rna_matrix.rds")
saveRDS(data2, "PreprocessedData/CD8_atac_matrix.rds")
saveRDS(genenames, "PreprocessedData/CD8_genename.rds")
saveRDS(cell_annotation[idx], "PreprocessedData/CD8_celltype.rds")
rm(list = ls())
