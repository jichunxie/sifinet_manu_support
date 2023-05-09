library(SiFINeT)
library(GSVA)
library(Seurat)
#setwd("~/Desktop/SiFINeT/Result/Numerical/SD2/")

set.seed(1)
genename <- readRDS("name.rds")
data <- readRDS("matrix.rds")
data <- as.matrix(data)
so <- create_SiFINeT_object(counts = data, gene.name = genename)
so <- quantile_thres(so)
so <- feature_coexp(so)
so <- create_network(so)
so <- filter_lowexp(so)
so <- cal_connectivity(so)
so <- find_unique_feature(so, t3p = 0.4)
so <- assign_shared_feature(so)
so <- enrich_feature_set(so)
saveRDS(so, "so.rds")
featureset <- so@featureset

feature_list <- list()
colnames(data) <- 1:ncol(data)
rownames(data) <- genename
so <- CreateSeuratObject(data)
so <- NormalizeData(so)
Y <- so@assays$RNA@data
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
