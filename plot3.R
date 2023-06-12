library(Seurat)
library(mclust)
library(grid)
library(ggpubr)
library(ggplot2)
library(M3C)
library(Matrix)
library(tidyverse)
library(ggforce)

setwd("~/Desktop/SiFINeT/Result_final")
setwd("Experimental/Monoclonal/")
set.seed(1)

sio <- readRDS("so.rds")
gs <- sio@featureset
rm(sio)

data <- readRDS("matrix.rds")
genename <- readRDS("genename.rds")
colnames(data) <- 1:ncol(data)
rownames(data) <- tolower(genename)
so <- CreateSeuratObject(data)

so <- NormalizeData(so)
so <- FindVariableFeatures(so, selection.method = "vst", nfeatures = length(unique(unlist(gs))))
all.genes <- rownames(so)
so <- ScaleData(so, features = all.genes)

gsva_mat <- readRDS("gsva_res.rds")
gene_set <- tolower(unique(unlist(gs)))

so <- RunPCA(so, features = unique(c(gene_set, VariableFeatures(object = so))))
so <- RunUMAP(so, dims = 1:50, n.components = 5, seed.use = 42)

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
so <- CellCycleScoring(so, s.features = tolower(s.genes), 
                       g2m.features = tolower(g2m.genes), set.ident = TRUE)

ccgroup <- so[[]]$Phase
gsva_mat2 <- rbind(gsva_mat, rep(0, ncol(gsva_mat)))
gsva_group <- apply(gsva_mat2[,], 2, which.max)
cluster_d <- readRDS("cluster.rds")
cluster_k <- readRDS("cluster_k.rds")
so <- AddMetaData(so, 3 - gsva_group, col.name = "gsva_cell_type")
so <- AddMetaData(so, ccgroup, col.name = "ccgroup")
so <- AddMetaData(so, cluster_d, col.name = "CIDRcluster_default")
so <- AddMetaData(so, cluster_k, col.name = "CIDRcluster")

p1 <- DimPlot(so, reduction = "umap", group.by = "ccgroup", dims = c(1,2), 
              cols = c("orange", "blue", "green")) + labs(title = NULL)
p2 <- DimPlot(so, reduction = "umap", group.by = "gsva_cell_type", dims = c(1,2), 
              cols = c("orange", "blue", "green")) + labs(title = NULL)
p3 <- DimPlot(so, reduction = "umap", group.by = "CIDRcluster", dims = c(1, 2), 
              cols = c("blue", "green", "orange")) + labs(title = NULL) + 
  geom_circle(aes(x0=1, y0=-0.8, r=0.8), col = "red", inherit.aes = FALSE) + 
  geom_circle(aes(x0=-1, y0=-2, r=0.8), col = "red", inherit.aes = FALSE)

adjustedRandIndex(gsva_group, ccgroup)
adjustedRandIndex(cluster_k, ccgroup)

p4 <- ggarrange(p2, p1, p3, nrow = 1, common.legend = TRUE, 
                legend = "right", legend.grob = get_legend(p1), 
                labels = c("a", "b", "c"))

psupp <- DimPlot(so, reduction = "umap", group.by = "CIDRcluster_default", dims = c(1, 2)) + labs(title = NULL)
psupp
ggsave("../../Supp_mono_clu.jpeg", width = 2.5, height = 2.5, 
       units = "in", device='jpeg', dpi=600)

setwd("../BoneMarrow")
set.seed(1)

# SiFINeT+GSVA
gsvares <- readRDS("gsva_res.rds")
gsvares <- rbind(gsvares[2,], gsvares[1,])
color1 <- ifelse(gsvares[1,] > 0, 1, ifelse(gsvares[2,] > 0, 2, 3))

gsvares1 <- readRDS("gsva_res_1.rds")
gsvares1_temp <- rbind(gsvares1, rep(0, ncol(gsvares1)))
color2 <- apply(gsvares1_temp, 2, which.max)
pcares <- prcomp(t(gsvares1))

plotdata1 <- data.frame(cbind(t(gsvares), color1))
colnames(plotdata1) <- c("GSVA_layer1_GS1", "GSVA_layer1_GS2", "CellGroup")
plotdata1$CellGroup <- factor(plotdata1$CellGroup, levels = 1:7, 
                              labels = c("Granulocyte/macrophage",
                                         "Megakaryocyte/erythrocyte",
                                         "Early progenitors",
                                         "Granulocyte",
                                         "Macrophage",
                                         "Monocyte/dentritic",
                                         "T/B"))
g1 <- ggplot() + 
  geom_point(data = plotdata1, aes(x = GSVA_layer1_GS1, y = GSVA_layer1_GS2, color = CellGroup), size = 1) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(), 
        axis.line = element_line(colour = "black"), 
        legend.text=element_text(size=10,face="bold"), 
        legend.title=element_text(size=10,face="bold")) +
  scale_color_manual(name = "Progenitor Cell Subpopulation",
                     values=c("#E69F00", "#56B4E9", "#009E73", 
                              "#F0E442", "#0072B2", "#D55E00", "#999999"),
                     drop = FALSE) + 
  guides(size="none") + 
  guides(color = guide_legend(override.aes = list(size=5)))#+ theme(legend.position = "none")

plotdata2 <- data.frame(cbind(pcares$x[,1:2], color2 + 3))
colnames(plotdata2) <- c("GSVA_layer2_PC1", "GSVA_layer2_PC2", "CellGroup")
plotdata2$CellGroup <- factor(plotdata2$CellGroup, levels = 1:7, 
                              labels = c("Granulocyte/macrophage",
                                         "Megakaryocyte/erythrocyte ",
                                         "Early progenitors",
                                         "Granulocyte",
                                         "Macrophage",
                                         "Monocyte/dentritic",
                                         "T/B"))
g2 <- ggplot() + 
  geom_point(data = plotdata2, aes(x = GSVA_layer2_PC1, y = GSVA_layer2_PC2, color = CellGroup), size = 1) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(), 
        axis.line = element_line(colour = "black"), 
        legend.text=element_text(size=10,face="bold"), 
        legend.title=element_text(size=10,face="bold")) +
  scale_color_manual(name = "Progenitor Cell Subpopulation",
                     values = c("#E69F00", "#56B4E9", "#009E73", 
                                "#F0E442", "#0072B2", "#D55E00", "#999999"),
                     drop = FALSE) + 
  guides(size="none") + 
  guides(color = guide_legend(override.aes = list(size=5))) + 
  theme(legend.position = "none")

g3 <- g1  + 
  annotate(geom = "segment", x = 0, xend = 0.65, y = -0.65, yend = -0.65) + 
  annotate(geom = "segment", x = 0, xend = 0.65, y = 0, yend = 0) + 
  annotate(geom = "segment", x = 0, xend = 0, y = 0, yend = -0.65) + 
  annotate(geom = "segment", x = 0.65, xend = 0.65, y = 0, yend = -0.65) + 
  annotate(geom = "segment", x = 0.65, xend = 0.75, y = 0, yend = 0.7) + 
  annotate(geom = "segment", x = 0.65, xend = 0.75, y = -0.65, yend = -0.7)

g5 <- ggarrange(g3, g2, nrow = 1, labels = c("d", "e"), common.legend = TRUE, legend="right") 

# Gene expression
table(color2)
table(color1)
cell_group <- c(unname(table(color2)), unname(table(color1))[3:2])

selected_gene <- c("Cebpa", "Elane", "Mpo", "Cebpe", "Ndufa4", "Sub1", "Irf8", "Ctss", 
                   "Cst3", "Cd34", "Cd27", "Gata2", "Apoe", "Gata1", "Gfi1b", "Klf1")
a <- order(color2)

data <- readRDS("matrix.rds")
data_temp <- data[match(selected_gene, rownames(data)), ]
data1 <- data_temp[, color1 == 1]
data1 <- data1[, a]
data2 <- data_temp[, color1 == 2]
data3 <- data_temp[, color1 == 3]

data_final <- cbind(data1, data3, data2)
a1 <- rownames(data_final)
a2 <- colnames(data_final)

rm(data1, data2, data3, data, data_temp)

data_final <- as.matrix(data_final)
data_final2 <- data_final
for (i in 1:16){
  data_final2[i, ] <- (data_final[i, ] - min(data_final[i, ]))/diff(range(data_final[i, ]))
}
selected_gene_label <- selected_gene
rownames(data_final2) <- selected_gene_label

data_final1 <- data_final2 %>% 
  as.data.frame() %>%
  rownames_to_column("Gene_name") %>%
  pivot_longer(-c(Gene_name), names_to = "Sample_name", values_to = "counts")
data_final1$Gene_name <- factor(data_final1$Gene_name, levels = selected_gene_label)
data_final1$Sample_name <- factor(data_final1$Sample_name, levels = a2)

g1 <- ggplot(data_final1, aes(Sample_name, Gene_name, fill= counts)) + 
  geom_tile() + 
  scale_fill_gradient2(low="navy", mid="white", high="red", midpoint = 0) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme(legend.text=element_text(size=10,face="bold"), 
        legend.title=element_text(size=12,face="bold"),
        axis.text.x=element_blank(), #remove x axis labels
        axis.title.x=element_blank(), 
        axis.title.y=element_blank(), 
        axis.ticks.x=element_blank(),
        plot.margin = unit(c(0.01, 0.01, 0.1, 0.03), "npc")) + theme(legend.position = "none") + 
  geom_vline(xintercept = cumsum(cell_group))
g1

g2 <- g1 + 
  annotate("text",x = c(0, cumsum(cell_group[1:(length(cell_group) - 1)])) + cell_group / 2, 
           y=-0.5, label = c("Granulocyte",
                             "Macrophage",
                             "Monocyte/dentritic",
                             "T/B",
                             "Early progenitors",
                             "Megakaryocyte/erythrocyte"), size = 3) + 
  coord_cartesian(ylim=c(1, 16), clip="off")


g7 <- ggarrange(p4, g5, g2, ncol = 1, labels = c("", "", "f"))
ggsave("../../Fig3.jpeg", width = 8, height = 9, 
       units = "in", device='jpeg', dpi=1200)

rm(list = ls())
