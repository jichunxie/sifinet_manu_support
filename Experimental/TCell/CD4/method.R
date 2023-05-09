library(SiFINeT)
library(GSVA)
library(Seurat)
#setwd("~/Desktop/SiFINeT/Result/Experimental/TCell/CD4")

# RNA
set.seed(1)
genename <- readRDS("../PreprocessedData/CD4_genename.rds")
data <- readRDS("../PreprocessedData/CD4_rna_matrix.rds")
so <- create_SiFINeT_object(counts = data, gene.name = genename)
so <- quantile_thres(so)
so <- feature_coexp(so)
so <- create_network(so)
so <- filter_lowexp(so)
so <- cal_connectivity(so)
so <- find_unique_feature(so, t2 = 0.2, t3 = 0.2, t2p = 0.2, t3p = 0.2)
so <- assign_shared_feature(so, 0.4)
so <- enrich_feature_set(so, 0.8)
saveRDS(so, "CD4_rna_so.rds")
featureset <- so@featureset
feature_list <- list()
so <- CreateSeuratObject(data)
so <- NormalizeData(so)
Y <- so@assays$RNA@data
for (i in 1:length(featureset$unique)){
  feature_list[[i]] <- match(c(featureset$unique[[i]], 
                               featureset$shared[[i]],
                               featureset$enriched[[i]]), 
                             genename)
}
outfile <- "gsva_res_rna.rds"
Y <- as.matrix(Y)
rownames(Y) <-  1:nrow(Y)
colnames(Y) <- paste0("s", 1:ncol(Y))
gsva.es <- gsva(Y, feature_list, parallel.sz = 8)
gsva.es.mat <- as.matrix(gsva.es)
saveRDS(gsva.es.mat, outfile)
rm(list = ls())



# ATAC
set.seed(1)
genename <- readRDS("../PreprocessedData/CD4_genename.rds")
data <- readRDS("../PreprocessedData/CD4_atac_matrix.rds")
so <- create_SiFINeT_object(counts = data, gene.name = genename)
so <- quantile_thres(so)
so <- feature_coexp(so)
so <- create_network(so, manual = TRUE, least_edge_prop = 0.005)
so <- filter_lowexp(so)
so <- cal_connectivity(so)
so <- find_unique_feature(so, t2 = 0.08, t3 = 0.08, t2p = 0.1, t3p = 0.1)
so <- assign_shared_feature(so, 0.4)
so <- enrich_feature_set(so, 0.8)
saveRDS(so, "CD4_atac_so.rds")
featureset <- so@featureset
feature_list <- list()
so <- CreateSeuratObject(data)
so <- NormalizeData(so)
Y <- so@assays$RNA@data
for (i in 1:length(featureset$unique)){
  feature_list[[i]] <- match(c(featureset$unique[[i]], 
                               featureset$shared[[i]],
                               featureset$enriched[[i]]), 
                             genename)
}
outfile <- "gsva_res_atac.rds"
Y <- as.matrix(Y)
rownames(Y) <-  1:nrow(Y)
colnames(Y) <- paste0("s", 1:ncol(Y))
gsva.es <- gsva(Y, feature_list, parallel.sz = 8)
gsva.es.mat <- as.matrix(gsva.es)
saveRDS(gsva.es.mat, outfile)
rm(list = ls())
