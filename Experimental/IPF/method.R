library(SiFINeT)
library(GSVA)
library(Seurat)
#setwd("~/Desktop/SiFINeT/Result/Experimental/IPF/")

set.seed(1)
genename <- readRDS("1_genename.rds")
data <- readRDS("1_matrix.rds")
data <- as.matrix(data)
so <- create_SiFINeT_object(counts = data, gene.name = genename)
so <- quantile_thres(so)
so <- feature_coexp(so)
so <- create_network(so)
so <- filter_lowexp(so)
so <- cal_connectivity(so)
so <- find_unique_feature(so, t2 = 0.3, t3 = 0.25, t2p = 0.6, resolution = 0.8)
so <- assign_shared_feature(so)
so <- enrich_feature_set(so)
saveRDS(so, "so.rds")
featureset <- so@featureset

feature_list <- list()
colnames(data) <- 1:ncol(data)
rownames(data) <- genename
so <- CreateSeuratObject(data)
so <- SCTransform(so, return.only.var.genes = FALSE)
Y <- so@assays$SCT@data
for (i in 1:length(featureset$unique)){
  feature_list[[i]] <- match(c(featureset$unique[[i]], 
                               featureset$shared[[i]],
                               featureset$enriched[[i]]), 
                             genename)
}
outfile <- "gsva_res.rds"
Y <- as.matrix(Y)
rownames(Y) <-  1:nrow(Y)
colnames(Y) <- paste0("s", 1:ncol(Y))
gsva.es <- gsva(Y, feature_list, parallel.sz = 8)
gsva.es.mat <- as.matrix(gsva.es)
saveRDS(gsva.es.mat, outfile)
rm(list = ls())
