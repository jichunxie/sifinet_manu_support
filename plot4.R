library(SiFINeT)
library(tidyverse)
library(ggplot2)
library(ggforce)
library(ggpubr)
library(ggraph)
library(jpeg)
library(grid)

setwd("~/Desktop/SiFINeT/Result_final/")
setwd("Numerical/SD2")
g1 <- readRDS("SD2_monocle3_plot.rds")
a <- readRDS("SD2_monocle3_UMAP.rds")

a$True_cell_cluster <- factor(rep(1:5, each = 1200), levels = 1:5)

middle_pos <- a %>% 
  group_by(True_cell_cluster) %>% 
  summarise(aa = mean(UMAP1), bb = mean(UMAP2))

segment_pos <- data.frame(xstart = middle_pos$aa[1:4], 
                          ystart = middle_pos$bb[1:4], 
                          xend = middle_pos$aa[2:5], 
                          yend = middle_pos$bb[2:5])


g2 <- ggplot() + 
  geom_point(data = a, aes(x = UMAP1, y = UMAP2, color = True_cell_cluster), size = 0.1) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(), 
        axis.line = element_line(colour = "black"), 
        legend.text=element_text(size=10,face="bold"), 
        legend.title=element_text(size=10,face="bold")) +
  guides(colour = guide_legend(override.aes = list(size=5))) + theme(legend.position = "none") + 
  geom_segment(data = segment_pos, aes(x = xstart, y = ystart, xend = xend, yend = yend), size = 2) + 
  guides(size="none")


so <- readRDS("so.rds")
geneset <- so@featureset
gene_list <- list()
for (i in 1:length(geneset$unique)){
  gene_list[[i]] <- match(c(geneset$unique[[i]], geneset$shared[[i]], geneset$enriched[[i]]), so@gene.name)
}

data2 <- readRDS("gsva_res.rds")
group1_a <- apply(data2, 2, which.max)

a$SifiNet_GSVA_cluster <- factor(group1_a, levels = 1:5)

middle_pos <- a %>% 
  group_by(SifiNet_GSVA_cluster) %>% 
  summarise(aa = mean(UMAP1), bb = mean(UMAP2))

group_id <- t(combn(1:5, 2))

segment_pos <- data.frame(xstart = middle_pos$aa[group_id[,1]], 
                          ystart = middle_pos$bb[group_id[,1]], 
                          xend = middle_pos$aa[group_id[,2]], 
                          yend = middle_pos$bb[group_id[,2]],
                          weight = 0.05 * apply(group_id, 1, function(x){length(intersect(gene_list[[x[1]]], gene_list[[x[2]]]))}))
segment_pos <- segment_pos[segment_pos$weight > 0,]

g3 <- ggplot() + 
  geom_point(data = a, aes(x = UMAP1, y = UMAP2, color = SifiNet_GSVA_cluster), size = 0.1) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(), 
        axis.line = element_line(colour = "black"), 
        legend.text=element_text(size=10,face="bold"), 
        legend.title=element_text(size=10,face="bold")) +
  theme(legend.position = "none") + 
  guides(colour = guide_legend(override.aes = list(size=5))) + 
  geom_segment(data = segment_pos, aes(x = xstart, y = ystart, xend = xend, yend = yend, size = weight)) + 
  scale_size(range=c(1,2)) + 
  guides(size="none")



setwd("../../Experimental/CD8")
set.seed(1)

gene_name <- readRDS("PreprocessedData/CD8_genename.rds")

gsva_mat <- readRDS("gsva_res_rna.rds")
gsva_mat_temp <- rbind(gsva_mat, rep(0, ncol(gsva_mat)))
gsva_group_rna <- apply(gsva_mat_temp[,], 2, which.max)
b <- max(gsva_group_rna)
if (b != nrow(gsva_mat)){
  gsva_group_rna[gsva_group_rna == b] <- "NA"
  gsva_group_rna <- factor(gsva_group_rna, levels = c(1:(b-1), "NA"), labels = paste("CD8_R_", c(1:(b-1), "NA"), sep = ""))
} else {
  gsva_group_rna <- factor(gsva_group_rna, levels = 1:b, labels = paste("CD8_R_", 1:b, sep = ""))
}

gsva_mat2 <- readRDS("gsva_res_atac.rds")
gsva_mat_temp2 <- rbind(gsva_mat2, rep(0, ncol(gsva_mat)))
gsva_group_atac <- apply(gsva_mat_temp2[,], 2, which.max)
b <- max(gsva_group_atac)
if (b != nrow(gsva_mat2)){
  gsva_group_atac[gsva_group_atac == b] <- "NA"
  gsva_group_atac <- factor(gsva_group_atac, levels = c(1:(b-1), "NA"), labels = paste("CD8_A_", c(1:(b-1), "NA"), sep = ""))
} else {
  gsva_group_atac <- factor(gsva_group_atac, levels = 1:b, labels = paste("CD8_A_", 1:b, sep = ""))
}

sum(gsva_group_atac == "CD8_A_NA") # 1
nNA_idx <- which(gsva_group_atac != "CD8_A_NA")

ccgroup_ori <- readRDS("PreprocessedData/CD8_celltype.rds")
ccgroup <- factor(ccgroup_ori, levels = c("CD8 Naive", "CD8 TEM_1", "CD8 TEM_2"), 
                  labels = c("Naive", "TEM_1", "TEM_2"))

pcares <- prcomp(t(rbind(gsva_mat[, nNA_idx], gsva_mat2[, nNA_idx])))

pcres <- pcares$x[, 1:2]
a <- data.frame(pcres)
a$SifiNet_rna_cell_cluster <- gsva_group_rna[nNA_idx]
a$SifiNet_rna_cell_cluster <- factor(a$SifiNet_rna_cell_cluster, 
                                     levels = c("CD8_R_1", "CD8_R_2", "CD8_R_3",
                                                "CD8_R_4", "CD8_R_NA"),
                                     labels = c("TEM_1", "TCM_1", "TEM_2", 
                                                "TNaive", "TCM_2"))
a$SifiNet_atac_cell_cluster <- gsva_group_atac[nNA_idx]
a$SifiNet_atac_cell_cluster <- factor(a$SifiNet_atac_cell_cluster,
                                      levels = c("CD8_A_1", "CD8_A_2"),
                                      labels = c("TCM", "TEM"))
a$celltype <- ccgroup[nNA_idx]


g4 <- ggplot() + 
  geom_point(data = a, aes(x = PC1, y = PC2, color = celltype), size = 0.5) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(), 
        axis.line = element_line(colour = "black"), 
        legend.text=element_text(face="bold"), 
        legend.title=element_text(face="bold"), 
        plot.title = element_text(hjust = 0.5, face="bold")) +
  scale_color_manual(values=c("#88CCEE", "#CC6677", "#DDCC77")) +
  guides(colour = guide_legend(override.aes = list(size=2))) + 
  guides(size="none") + labs(color=guide_legend(title="Seurat"))


g5 <- ggplot() + 
  geom_point(data = a, aes(x = PC1, y = PC2, color = SifiNet_rna_cell_cluster), size = 0.5) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(), 
        axis.line = element_line(colour = "black"), 
        legend.text=element_text(face="bold"), 
        legend.position = "right", 
        legend.justification = "bottom",
        legend.title=element_text(face="bold"), 
        plot.title = element_text(hjust = 0.5, face="bold")) +
  guides(colour = guide_legend(override.aes = list(size=2))) + 
  scale_color_manual(values=c("#117733", "#332288", "#AA4499", 
                              "#44AA99", "#888888")) +
  scale_size(range=c(0.5,1)) + 
  guides(size="none") + coord_cartesian(clip = "off") + 
  labs(color=guide_legend(title="SifiNet RNA"))


g6 <- ggplot() + 
  geom_point(data = a, aes(x = PC1, y = PC2, color = SifiNet_atac_cell_cluster), size = 0.5) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(), 
        axis.line = element_line(colour = "black"), 
        legend.text=element_text(face="bold"), 
        legend.position = "right", 
        legend.justification = "bottom", 
        legend.title=element_text(face="bold"), 
        plot.title = element_text(hjust = 0.5, face="bold")) +
  guides(colour = guide_legend(override.aes = list(size=2))) + 
  scale_color_manual(values=c("#999933", "#882255", "#888888")) +
  scale_size(range=c(0.5,1)) + 
  guides(size="none") + coord_cartesian(clip = "off") + 
  labs(color=guide_legend(title="SifiNet ATAC"))


so <- readRDS("CD8_rna_so.rds")
id_list <- list()
for (i in 1:length(so@featureset$unique)){
  id_list[[i]] <- match(c(so@featureset$unique[[i]], 
                          so@featureset$shared[[i]], 
                          so@featureset$enriched[[i]]), so@gene.name)
}
edge <- data.frame(t(combn(1:length(id_list), 2)))
edge_mat <- 1 * (abs(so@coexp - so@est_ms$mean) >= so@thres)
edge$edge_weight <- apply(edge, 1, function(x){mean(edge_mat[id_list[[x[1]]], id_list[[x[2]]]])})
edge <- edge[edge$edge_weight > 0.2, ]
edge_width_range <- range(edge$edge_weight)

gt1 <- geneset_topology(so, weightthres = 0.2, 
                        set_name = c("TEM_1", "TCM", "TEM_2", "TNaive"), 
                        node_color = c("#117733", "#332288", "#AA4499", "#44AA99"))
gt1
ggsave("gt1.jpeg", device='jpeg', dpi=1200, width = 3, height = 2)
image1 <- jpeg::readJPEG("gt1.jpeg")

so <- readRDS("CD8_atac_so.rds")
id_list <- list()
for (i in 1:length(so@featureset$unique)){
  id_list[[i]] <- match(c(so@featureset$unique[[i]], 
                          so@featureset$shared[[i]], 
                          so@featureset$enriched[[i]]), so@gene.name)
}
edge <- data.frame(t(combn(1:length(id_list), 2)))
edge_mat <- 1 * (abs(so@coexp - so@est_ms$mean) >= so@thres)
edge$edge_weight <- apply(edge, 1, function(x){mean(edge_mat[id_list[[x[1]]], id_list[[x[2]]]])})
edge <- edge[edge$edge_weight > 0.2, ]
new_width <- 0.5 + 1.5 * 
  (edge$edge_weight[1] - edge_width_range[1]) / 
  (edge_width_range[2] - edge_width_range[1])
gt2 <- geneset_topology(so, weightthres = 0.2, set_name = c("TCM", "TEM"), 
                        node_color = c("#117733", "#332288"))
gt2 + scale_edge_width(range=c(new_width, new_width))
ggsave("gt2.jpeg", device='jpeg', dpi=1200, width = 3, height = 2)
image2 <- jpeg::readJPEG("gt2.jpeg")


g7 <- g5 + annotation_custom(grob = rasterGrob(image1, interpolate = TRUE, 
                                               width=unit(0.8,'npc'),
                                               x = unit(0.85,"npc"), y = unit(0.6,"npc"),
                                               hjust = 0, vjust = 0)) +
  coord_cartesian(clip = "off")
g7

g8 <- g6 + annotation_custom(grob = rasterGrob(image2, interpolate = TRUE, 
                                               width=unit(0.8,'npc'),
                                               x = unit(0.93,"npc"), y = unit(0.6,"npc"),
                                               hjust = 0, vjust = 0)) +
  coord_cartesian(clip = "off")
g8


setwd("../../")
ggarrange(g3, g2, g1, g4, g7, g8,
          nrow = 2, ncol = 3,
          labels = "auto")
ggsave("Fig4.jpeg", width = 12, height = 6, 
       units = "in", device='jpeg', dpi=1200)
