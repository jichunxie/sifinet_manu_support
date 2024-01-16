library(tidyverse)
library(monocle3)
library(ggplot2)

set.seed(1)
Y <- readRDS("matrix.rds")
idx <- readRDS("name.rds")

gene_name <- idx
counts <- Y
rownames(counts) <- gene_name
colnames(counts) <- 1:6000
cell_metadata <- matrix(0, 6000, 1) %>% data.frame()
cell_metadata[,1] <- factor(rep(1:5, each = 1200))
rownames(cell_metadata) <- 1:6000
colnames(cell_metadata) <- c("cell_type")
gene_metadata <- matrix(gene_name, nrow(counts), 1)
rownames(gene_metadata) <- gene_name
colnames(gene_metadata) <- c("gene_short_name")

cds <- new_cell_data_set(counts,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_metadata)

cds <- preprocess_cds(cds, num_dim = 100)

cds <- reduce_dimension(cds, reduction_method = "UMAP")

cds <- cluster_cells(cds, k=20)

saveRDS(cds, "SD2_monocle3_cds.rds")
cds <- learn_graph(cds, use_partition = F)
g1 <- plot_cells(cds, color_cells_by = "cell_type", group_label_size = 10, label_cell_groups = F, trajectory_graph_segment_size = 2, 
                 label_leaves = F, label_roots = F, label_branch_points = F) + theme(legend.position = "none")
g1 <- g1 + xlab("UMAP1") + ylab("UMAP2")
saveRDS(g1, "SD2_monocle3_plot.rds")

a <- data.frame(reducedDims(cds)[["UMAP"]])
colnames(a) <- c("UMAP1", "UMAP2")
saveRDS(a, "SD2_monocle3_UMAP.rds")
rm(list = ls())


