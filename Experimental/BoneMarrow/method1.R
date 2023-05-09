library(SiFINeT)
library(GSVA)
library(Seurat)
#setwd("~/Desktop/SiFINeT/Result/Experimental/BoneMarrow")

# run SiFINeT
set.seed(1)

data <- readRDS("matrix_1.rds")
genename <- readRDS("genename_1.rds")

so <- create_SiFINeT_object(counts = data, gene.name = genename)
so <- quantile_thres(so)
so <- feature_coexp(so)
so <- create_network(so)
so <- filter_lowexp(so)
so <- cal_connectivity(so)
so <- find_unique_feature(so, t2 = 0.3, t3 = 0.2, t2p = 0.2, t3p = 0.2, resolution = 0.8)
so <- assign_shared_feature(so, 0.4)
so <- enrich_feature_set(so, 0.4)
saveRDS(so@featureset, "featureset_1.rds")
saveRDS(so, "so_1.rds")
rm(list = ls())

# GSVA
feature_list <- list()
data <- readRDS("matrix_1.rds")
so <- CreateSeuratObject(data)
so <- NormalizeData(so)
Y <- so@assays$RNA@data
genename <- readRDS("genename_1.rds")
featureset <- readRDS("featureset_1.rds")
for (i in 1:length(featureset$unique)){
  feature_list[[i]] <- match(c(featureset$unique[[i]], 
                               featureset$shared[[i]],
                               featureset$enriched[[i]]), 
                             genename)
}
outfile <- "gsva_res_1.rds"
Y <- as.matrix(Y)
rownames(Y) <-  1:nrow(Y)
colnames(Y) <- paste0("s", 1:ncol(Y))
gsva.es <- gsva(Y, feature_list, parallel.sz = 8)
gsva.es.mat <- as.matrix(gsva.es)
saveRDS(gsva.es.mat, outfile)
rm(list = ls())
