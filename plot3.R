library(Seurat)
library(mclust)
library(grid)
library(ggpubr)
library(ggplot2)
library(M3C)
library(Matrix)
library(tidyverse)
library(ggforce)
library(caret)

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
gsva_group_o <- ifelse(gsva_group == 3, 0, ifelse(gsva_group == 1, 1, 2))
cluster_d <- readRDS("cluster.rds")
cluster_k <- readRDS("cluster_k.rds")
so <- AddMetaData(so, gsva_group_o, col.name = "gsva_cell_type")
so <- AddMetaData(so, ccgroup, col.name = "ccgroup")
so <- AddMetaData(so, cluster_d, col.name = "CIDRcluster_default")
so <- AddMetaData(so, cluster_k, col.name = "CIDRcluster")

p1 <- DimPlot(so, reduction = "umap", group.by = "ccgroup", dims = c(1,2), 
              cols = c("orange", "blue", "green")) + labs(title = "Seurat score")
p2 <- DimPlot(so, reduction = "umap", group.by = "gsva_cell_type", dims = c(1,2), 
              cols = c("orange", "green", "blue")) + labs(title = "SifiNet")
p2n <- DimPlot(so, reduction = "umap", group.by = "gsva_cell_type", dims = c(1,2)) +
  scale_color_manual(labels=factor(c("Early G1","Late G1/S","G2/M"), levels = c("Early G1","Late G1/S","G2/M")),
                     values=c("orange", "green", "blue")) + labs(title = "SifiNet")
p3 <- DimPlot(so, reduction = "umap", group.by = "CIDRcluster", dims = c(1, 2), 
              cols = c("blue", "green", "orange")) + labs(title = "CIDR clustering") + 
  geom_circle(aes(x0=1, y0=-0.8, r=0.8), col = "red", inherit.aes = FALSE) + 
  geom_circle(aes(x0=-1, y0=-2, r=0.8), col = "red", inherit.aes = FALSE)

adjustedRandIndex(gsva_group, ccgroup)
adjustedRandIndex(cluster_k, ccgroup)

p4 <- ggarrange(p2n, p1, p3, nrow = 1, common.legend = TRUE, 
                legend = "right", legend.grob = get_legend(p2n), 
                labels = c("a", "b", "c"))

psupp <- DimPlot(so, reduction = "umap", group.by = "CIDRcluster_default", dims = c(1, 2)) + labs(title = NULL)
ggsave("../../Supp_mono_clu.pdf", width = 2.5, height = 2.5, 
       units = "in", device='pdf', dpi=1200)

ccgroupn <- ifelse(ccgroup == "G1", 0, ifelse(ccgroup == "G2M", 2, 1))
conf1 <- confusionMatrix(factor(gsva_group_o, levels = 2:0), factor(ccgroupn, levels = 2:0), 
                         dnn = c("Prediction", "Reference"))
plt <- as.data.frame(conf1$table)
plt$Prediction <- factor(plt$Prediction, levels=rev(levels(plt$Prediction)))
ari1 <- adjustedRandIndex(gsva_group, ccgroupn)
pconf1 <- ggplot(plt, aes(Prediction,Reference, fill= Freq)) +
  geom_tile() + geom_text(aes(label=Freq)) +
  scale_fill_gradient(low="white", high="#009194") +
  labs(x = "SifiNet",y = "Seurat score") +
  ggtitle(paste("ARI = ", round(ari1, 3), sep = "")) +
  scale_x_discrete(labels=c("Early G1","Late G1/S","G2/M")) +
  scale_y_discrete(labels=c("G2/M","Late G1/S","Early G1")) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5))

cidrn <- ifelse(cluster_k == 1, 2, ifelse(cluster_k == 2, 1, 0))
conf2 <- confusionMatrix(factor(cidrn, levels = 2:0), factor(ccgroupn, levels = 2:0), 
                         dnn = c("Prediction", "Reference"))
plt <- as.data.frame(conf2$table)
plt$Prediction <- factor(plt$Prediction, levels=rev(levels(plt$Prediction)))
ari2 <- adjustedRandIndex(cidrn, ccgroupn)
pconf2 <- ggplot(plt, aes(Prediction,Reference, fill= Freq)) +
  geom_tile() + geom_text(aes(label=Freq)) +
  scale_fill_gradient(low="white", high="#009194") +
  labs(x = "CIDR",y = "Seurat score") +
  ggtitle(paste("ARI = ", round(ari2, 3), sep = "")) +
  scale_x_discrete(labels=c("Early G1","Late G1/S","G2/M")) +
  scale_y_discrete(labels=c("G2/M","Late G1/S","Early G1")) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5))

conf3 <- confusionMatrix(factor(gsva_group_o, levels = 2:0), factor(cidrn, levels = 2:0), 
                         dnn = c("Prediction", "Reference"))
plt <- as.data.frame(conf3$table)
plt$Prediction <- factor(plt$Prediction, levels=rev(levels(plt$Prediction)))
ari3 <- adjustedRandIndex(gsva_group, cidrn)
pconf3 <- ggplot(plt, aes(Prediction,Reference, fill= Freq)) +
  geom_tile() + geom_text(aes(label=Freq)) +
  scale_fill_gradient(low="white", high="#009194") +
  labs(x = "SifiNet",y = "CIDR") +
  ggtitle(paste("ARI = ", round(ari3, 3), sep = "")) +
  scale_x_discrete(labels=c("Early G1","Late G1/S","G2/M")) +
  scale_y_discrete(labels=c("G2/M","Late G1/S","Early G1")) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5))

pconf <- ggarrange(pconf1, pconf2, pconf3, nrow = 1, 
                   labels = c("a", "b", "c"))
pconf
ggsave("../../Supp_mono_confusion.pdf", width = 8, height = 2.5, 
       units = "in", device='pdf', dpi=1200)
setwd("../BoneMarrow")
set.seed(1)

# SiFINeT+GSVA
gsvares <- readRDS("gsva_res.rds")
gsvares <- rbind(gsvares[2,], gsvares[1,])
color1 <- ifelse(gsvares[1,] > 0, 1, ifelse(gsvares[2,] > 0, 2, 3))
pcares <- prcomp(t(gsvares))

gsvares1 <- readRDS("gsva_res_1.rds")
gsvares1_temp <- rbind(gsvares1, rep(0, ncol(gsvares1)))
color2 <- apply(gsvares1_temp, 2, which.max)
pcares1 <- prcomp(t(gsvares1))

boundary <- matrix(c(0, 0, 0, -0.65, 0.65, -0.65, 0.65, -0.4), ncol = 2, byrow = T)
b <- (boundary - matrix(pcares$center, nrow = nrow(boundary), ncol = 2, byrow = T)) %*% pcares$rotation


plotdata1 <- data.frame(cbind(pcares$x[,1:2], color1))
colnames(plotdata1) <- c("GSVA_layer1_PC1", "GSVA_layer1_PC2", "CellGroup")
plotdata1$CellGroup <- factor(plotdata1$CellGroup, levels = 1:7, 
                              labels = c("Granulocyte/macrophage",
                                         "Megakaryocyte/erythrocyte",
                                         "Early myeloid progenitors",
                                         "Granulocyte",
                                         "Macrophage",
                                         "Monocyte/dendritic",
                                         "Lymphoid progenitors"))
g1 <- ggplot() + 
  geom_point(data = plotdata1, aes(x = GSVA_layer1_PC1, y = GSVA_layer1_PC2, color = CellGroup), size = 1) +
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

plotdata2 <- data.frame(cbind(pcares1$x[,1:2], color2 + 3))
colnames(plotdata2) <- c("GSVA_layer2_PC1", "GSVA_layer2_PC2", "CellGroup")
plotdata2$CellGroup <- factor(plotdata2$CellGroup, levels = 1:7, 
                              labels = c("Granulocyte/macrophage",
                                         "Megakaryocyte/erythrocyte ",
                                         "Early myeloid progenitors",
                                         "Granulocyte",
                                         "Macrophage",
                                         "Monocyte/dendritic",
                                         "Lymphoid progenitors"))
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
  annotate(geom = "segment", x = b[1,1], xend = b[2,1], y = b[1,2], yend = b[2,2]) + 
  annotate(geom = "segment", x = b[2,1], xend = b[3,1], y = b[2,2], yend = b[3,2]) + 
  annotate(geom = "segment", x = b[3,1], xend = b[4,1], y = b[3,2], yend = b[4,2]) + 
  annotate(geom = "segment", x = b[4,1], xend = b[1,1], y = b[4,2], yend = b[1,2]) + 
  annotate(geom = "segment", x = b[2,1], xend = 1.1, y = b[2,2], yend = 0.4) + 
  annotate(geom = "segment", x = b[1,1], xend = 1.1, y = b[1,2], yend = -0.3)

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
        plot.margin = unit(c(0.1, 0.01, 0.15, 0.03), "npc")) + theme(legend.position = "none") + 
  geom_vline(xintercept = cumsum(cell_group))
g1

g2 <- g1 + 
  annotate("text",x = c(0, cumsum(cell_group[1:(length(cell_group) - 1)])) + cell_group / 2, 
           y=-0.9, label = c("Granulocyte",
                             "Macrophage",
                             "Monocyte/dendritic",
                             "Lymphoid progenitors",
                             "Early myeloid progenitors",
                             "Megakaryocyte/erythrocyte"), size = 3, angle = 10) + 
  coord_cartesian(ylim=c(1, 16), clip="off")


g7 <- ggarrange(p4, g5, g2, ncol = 1, labels = c("", "", "f"))
ggsave("../../Fig3.pdf", width = 8, height = 9, 
       units = "in", device='pdf', dpi=1200)

rm(list = ls())
