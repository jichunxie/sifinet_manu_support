library(Matrix)
#setwd("~/Desktop/SiFINeT/Result/Experimental/Monoclonal/")

bulk <- read.delim("CT26WT1/CT26_bulk_rpm.txt", header = F, sep = " ", stringsAsFactors = F, skip = 1)
colnames(bulk) <- c("gene_id", "WT1", "WT2", "M7", "M16", "gene_name")


sum_bulk <- rowSums(bulk[,2:5] == 0)
selected_gene <- bulk[sum_bulk != 4, c(1,6)]

rm(bulk, sum_bulk)

mito_genes <- unlist(read.table("CT26WT1/mito_genes.txt", stringsAsFactors = F))

feature <- read.delim("CT26WT1/features.tsv", header = F, stringsAsFactors = F)
temp <- match(feature$V1, selected_gene$gene_id)
idx <- which(!is.na(temp))
gene_name <- feature$V2[idx]
saveRDS(gene_name, "all_gene_name.rds")
rm(feature, temp)
data_mat <- readMM("CT26WT1/matrix.mtx")
data_mat <- data_mat[idx,]

mt_idx <- match(tolower(mito_genes), tolower(gene_name))
rpc <- colSums(data_mat)
rpc_mt <- colSums(data_mat[mt_idx,])
mt_ratio <- rpc_mt / rpc
data_mat <- data_mat[, mt_ratio <= 0.05]

rpc <- colSums(data_mat)
data_mat <- data_mat[, rpc <= 1.8 * mean(rpc)]

data_mat <- data_mat[, colSums(data_mat > 0) > 0.1 * nrow(data_mat)]
id2 <- rowSums(data_mat > 0) > 0.1 * ncol(data_mat)
data_mat <- data_mat[id2,]
gene_name <- gene_name[id2]

saveRDS(data_mat, "matrix.rds")
saveRDS(gene_name, "genename.rds")
rm(list = ls())
