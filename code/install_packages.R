# R 4.0.3  (4.0.5 better)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("edgeR", force = T)
BiocManager::install("DESeq2")
install.packages('Seurat')
BiocManager::install("MAST")
BiocManager::install("scDD")
BiocManager::install("scde")
BiocManager::install("monocle")

install.packages("Matrix")
#BiocManager::install(version = "3.12")
BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor', 'Matrix.utils'))
BiocManager::install("cole-trapnell-lab/monocle3", force = TRUE) # require Rtools

BiocManager::install("TSCAN", force = TRUE)