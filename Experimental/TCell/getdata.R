library(SeuratData)
library(Seurat)
library(Signac)

pbmc.atac <- LoadData("pbmcMultiome", "pbmc.atac")
pbmc.atac@assays$ATAC@key <- "atac_"
saveRDS(pbmc.atac@assays$ATAC@data, "OriginalData/multi_atac_data.rds")

pbmc.rna <- LoadData("pbmcMultiome", "pbmc.rna")
saveRDS(pbmc.rna@assays$RNA@data, "OriginalData/multi_rna_data.rds")
saveRDS(pbmc.rna$seurat_annotations, "OriginalData/multi_annotation.rds")