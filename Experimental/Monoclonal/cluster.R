library(Seurat)
library(cidr)

set.seed(1)

data <- readRDS("matrix.rds")
colnames(data) <- 1:ncol(data)
rownames(data) <- 1:nrow(data)
so <- CreateSeuratObject(data)
  
so <- NormalizeData(so)
so <- FindVariableFeatures(so, selection.method = "vst", nfeatures = 500)
all.genes <- rownames(so)
so <- ScaleData(so, features = all.genes)


sData <- scDataConstructor(as.matrix(so@assays$RNA@counts), tagType = "raw")
rm(so)
sData <- determineDropoutCandidates(sData)
sData <- wThreshold(sData)
sData <- scDissim(sData)
sData <- scPCA(sData)
sData <- nPC(sData)
maxclu <- max(20, sData@nPC)
sData1 <- scCluster(sData, n = maxclu)

saveRDS(as.numeric(sData1@clusters), "cluster.rds")

sData2 <- scCluster(sData, nCluster = 3)

saveRDS(as.numeric(sData2@clusters), "cluster_k.rds")
