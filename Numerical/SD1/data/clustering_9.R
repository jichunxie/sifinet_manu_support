library(Seurat)
library(igraph)
library(cidr)

knn_self <- function(distmat, k){
  out <- matrix(0, nrow(distmat), ncol(distmat))
  for (i in 1:nrow(distmat)){
    temp <- order(distmat[i, ], decreasing = F)
    out[i, distmat[i, ] <= distmat[i, temp[k]]] <- 1
  }
  return(out)
}

for (i in 81:90){
  set.seed(i)
  data <- readRDS(paste(i, "_matrix.rds", sep = ""))
  
  colnames(data) <- 1:ncol(data)
  rownames(data) <- 1:nrow(data)
  so <- CreateSeuratObject(data)
  
  so <- NormalizeData(so)
  so <- FindVariableFeatures(so, selection.method = "vst", nfeatures = 500)
  all.genes <- rownames(so)
  so <- ScaleData(so, features = all.genes)
  so <- RunPCA(so, features = VariableFeatures(object = so))
  
  so <- FindNeighbors(so, dims = 1:10)
  saveRDS(rep(1:3, c(100, 4900, 1000)), paste(i, "_cluster_true.rds", sep = ""))
  so <- FindClusters(so)
  saveRDS(as.numeric(so$seurat_clusters), paste(i, "_cluster_seurat.rds", sep = ""))
  
  

  
  d <- dist(t(GetAssayData(so, slot = "scale.data")))
  g <- graph_from_adjacency_matrix(knn_self(as.matrix(d), 20), diag = F, mode = "undirected")
  lou_res <- cluster_louvain(g)
  saveRDS(as.numeric(lou_res$membership), paste(i, "_cluster_louvain.rds", sep = ""))
  
  
  sData <- scDataConstructor(as.matrix(so@assays$RNA@counts), tagType = "raw")
  rm(so)
  sData <- determineDropoutCandidates(sData)
  sData <- wThreshold(sData)
  sData <- scDissim(sData)
  sData <- scPCA(sData)
  sData <- nPC(sData)
  maxclu <- max(20, sData@nPC)
  sData <- scCluster(sData, n = maxclu)
  
  saveRDS(as.numeric(sData@clusters), paste(i, "_cluster_cidr.rds", sep = ""))
}
