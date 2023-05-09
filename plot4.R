library(SiFINeT)
library(tidyverse)
library(ggplot2)
library(ggforce)
library(ggpubr)
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

segment_pos <- data.frame(xstart = middle_pos$aa[1:4], 
                          ystart = middle_pos$bb[1:4], 
                          xend = middle_pos$aa[2:5], 
                          yend = middle_pos$bb[2:5],
                          weight = 0.05 * sapply(1:4, function(x){length(intersect(gene_list[[x]], gene_list[[x+1]]))}))


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



setwd("../../Experimental/TCell/CD8")
set.seed(1)

gene_name <- readRDS("../PreprocessedData/CD8_genename.rds")

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

ccgroup <- readRDS("../PreprocessedData/CD8_celltype.rds")

pcares <- prcomp(t(rbind(gsva_mat, gsva_mat2)))

pcres <- pcares$x[, 1:2]
a <- data.frame(pcres)
a$SifiNet_rna_cell_cluster <- gsva_group_rna
a$SifiNet_rna_cell_cluster <- factor(a$SifiNet_rna_cell_cluster, 
                                     labels = c("CD8 TEM_1",
                                                "CD8 TCM",
                                                "CD8 TEM_2",
                                                "CD8 TNaive",
                                                "Unassigned"))
a$SifiNet_atac_cell_cluster <- gsva_group_atac
a$SifiNet_atac_cell_cluster <- factor(a$SifiNet_atac_cell_cluster,
                                      labels = c("CD8 TCM",
                                                 "CD8 TEM",
                                                 "Unassigned"))
a$celltype <- ccgroup


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
  labs(color=guide_legend(title="SifiNet-RNA"))


g6 <- ggplot() + 
  geom_point(data = a, aes(x = PC1, y = PC2, color = SifiNet_atac_cell_cluster), size = 0.5) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(), 
        axis.line = element_line(colour = "black"), 
        legend.text=element_text(face="bold"), 
        legend.title=element_text(face="bold"), 
        plot.title = element_text(hjust = 0.5, face="bold")) +
  guides(colour = guide_legend(override.aes = list(size=2))) + 
  scale_color_manual(values=c("#999933", "#882255", "#888888")) +
  scale_size(range=c(0.5,1)) + 
  guides(size="none") + 
  labs(color=guide_legend(title="SifiNet-ATAC"))






setwd("../CD4")
set.seed(1)

gene_name <- readRDS("../PreprocessedData/CD4_genename.rds")

gsva_mat <- readRDS("gsva_res_rna.rds")
gsva_mat_temp <- rbind(gsva_mat, rep(0, ncol(gsva_mat)))
gsva_group_rna <- apply(gsva_mat_temp[,], 2, which.max)
b <- max(gsva_group_rna)
if (b != nrow(gsva_mat)){
  gsva_group_rna[gsva_group_rna == b] <- "NA"
  gsva_group_rna <- factor(gsva_group_rna, levels = c(1:(b-1), "NA"), labels = paste("CD4_R_", c(1:(b-1), "NA"), sep = ""))
} else {
  gsva_group_rna <- factor(gsva_group_rna, levels = 1:b, labels = paste("CD4_R_", 1:b, sep = ""))
}

gsva_mat2 <- readRDS("gsva_res_atac.rds")
gsva_mat_temp2 <- rbind(gsva_mat2, rep(0, ncol(gsva_mat)))
gsva_group_atac <- apply(gsva_mat_temp2[,], 2, which.max)
b <- max(gsva_group_atac)
if (b != nrow(gsva_mat2)){
  gsva_group_atac[gsva_group_atac == b] <- "NA"
  gsva_group_atac <- factor(gsva_group_atac, levels = c(1:(b-1), "NA"), labels = paste("CD4_A_", c(1:(b-1), "NA"), sep = ""))
} else {
  gsva_group_atac <- factor(gsva_group_atac, levels = 1:b, labels = paste("CD4_A_", 1:b, sep = ""))
}

ccgroup <- readRDS("../PreprocessedData/CD4_celltype.rds")

pcares <- prcomp(t(rbind(gsva_mat, gsva_mat2)))

pcres <- pcares$x[, 1:2]
a <- data.frame(pcres)
a$SifiNet_rna_cell_cluster <- gsva_group_rna
a$SifiNet_rna_cell_cluster <- factor(a$SifiNet_rna_cell_cluster,
                                     labels = c("CD4 TCM",
                                                "CD4 Treg",
                                                "CD4 TNaive",
                                                "Unassigned"))
a$SifiNet_atac_cell_cluster <- gsva_group_atac
a$SifiNet_atac_cell_cluster <- factor(a$SifiNet_atac_cell_cluster,
                                      labels = c("CD4_A_C1",
                                                 "CD4_A_C2",
                                                 "CD4_A_C3"))
a$celltype <- ccgroup


g7 <- ggplot() + 
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
  guides(size="none") + 
  labs(color=guide_legend(title="Seurat"))


g8 <- ggplot() + 
  geom_point(data = a, aes(x = PC1, y = PC2, color = SifiNet_rna_cell_cluster), size = 0.5) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(), 
        axis.line = element_line(colour = "black"), 
        legend.position = "right", 
        legend.justification = "bottom",
        legend.text=element_text(face="bold"), 
        legend.title=element_text(face="bold"), 
        plot.title = element_text(hjust = 0.5, face="bold")) +
  guides(colour = guide_legend(override.aes = list(size=2))) + 
  scale_color_manual(values=c("#117733", "#332288", "#AA4499", "#888888")) +
  scale_size(range=c(0.5,1)) + 
  guides(size="none") + 
  labs(color=guide_legend(title="SifiNet-RNA"))


g9 <- ggplot() + 
  geom_point(data = a, aes(x = PC1, y = PC2, color = SifiNet_atac_cell_cluster), size = 0.5) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(), 
        axis.line = element_line(colour = "black"), 
        legend.text=element_text(face="bold"), 
        legend.title=element_text(face="bold"), 
        plot.title = element_text(hjust = 0.5, face="bold")) +
  guides(colour = guide_legend(override.aes = list(size=2))) + 
  scale_color_manual(values=c("#44AA99", "#999933", "#882255", "#888888")) +
  scale_size(range=c(0.5,1)) + 
  guides(size="none") + 
  labs(color=guide_legend(title="SifiNet-ATAC"))


setwd("../CD8")
gene_name <- readRDS("../PreprocessedData/CD8_genename.rds")
so <- readRDS("CD8_rna_so.rds")
gt1 <- geneset_topology(so, weightthres = 0.2, set_name = c("CD8 TEM_1",
                                                            "CD8 TCM",
                                                            "CD8 TEM_2",
                                                            "CD8 TNaive"), 
                        node_color = c("#117733", "#332288", "#AA4499", "#44AA99"))
gt1
ggsave("gt1.jpeg", device='jpeg', dpi=1200, width = 3, height = 2)
image1 <- jpeg::readJPEG("gt1.jpeg")

setwd("../CD4")
gene_name <- readRDS("../PreprocessedData/CD4_genename.rds")
so <- readRDS("CD4_rna_so.rds")
gt2 <- geneset_topology(so, set_name = c("CD4 TCM",
                                         "CD4 Treg",
                                         "CD4 TNaive"), 
                        node_color = c("#117733", "#332288", "#AA4499"),
                        shiftsize = 0.2)
gt2
ggsave("gt2.jpeg", device='jpeg', dpi=1200, width = 3, height = 1)
image2 <- jpeg::readJPEG("gt2.jpeg")


g10 <- g5 + annotation_custom(grob = rasterGrob(image1, interpolate = TRUE, 
                                               width=unit(0.8,'npc'),
                                               x = unit(0.93,"npc"), y = unit(0.6,"npc"),
                                               hjust = 0, vjust = 0)) +
  coord_cartesian(clip = "off")
g10

g11 <- g8 + annotation_custom(grob = rasterGrob(image2, interpolate = TRUE, 
                                                width=unit(0.8,'npc'),
                                                x = unit(0.93,"npc"), y = unit(0.6,"npc"),
                                                hjust = 0, vjust = 0)) +
  coord_cartesian(clip = "off")
g11

setwd("../../../")
ggarrange(g3, g2, g1, g4, g10, g6, g7, g11, g9, 
          nrow = 3, ncol = 3,
          labels = "auto")
ggsave("Fig4.jpeg", width = 12, height = 9, 
       units = "in", device='jpeg', dpi=1200)
