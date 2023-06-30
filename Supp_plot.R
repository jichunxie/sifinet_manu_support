library(Seurat)
library(ggplot2)
library(ggpubr)
library(mclust)
library(tidyverse)
library(cogena)
library(Matrix)
library(ggnewscale)
library(stringr)

set.seed(1)


# Supp_sd1_res
# In plot2.R

# Supp_sd1_umap
setwd("~/Desktop/SiFINeT/Result_final/")
setwd("Numerical/SD1/")
data <- readRDS("data/18_matrix.rds")
group1 <- readRDS("data/18_cluster_true.rds")
group2 <- readRDS("data/18_cluster_seurat.rds")
colnames(data) <- 1:ncol(data)
rownames(data) <- 1:nrow(data)
so <- CreateSeuratObject(data)

so <- NormalizeData(so)
so <- FindVariableFeatures(so, selection.method = "vst", nfeatures = 500)
all.genes <- rownames(so)
so <- ScaleData(so, features = all.genes)
so <- RunPCA(so, features = VariableFeatures(object = so))

so <- FindNeighbors(so, dims = 1:10)
so <- FindClusters(so)
so <- RunUMAP(so, dims = 1:10)
so <- AddMetaData(so, data.frame(group1, group2))

g1 <- DimPlot(so, reduction = "umap", group.by = "group1") + 
  labs(title = NULL) + 
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()
  )
g2 <- DimPlot(so, reduction = "umap", group.by = "group2") + 
  labs(title = NULL) + 
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()
  )

ggarrange(g1, g2, nrow = 1, labels = "auto")
ggsave("../../Supp_sd1_umap.jpeg", width = 6, height = 2.5, 
       units = "in", device='jpeg', dpi=600)

# Supp_mono_clu
# In plot3

# Supp_cd8_trans
setwd("../../Experimental/CD8")
set.seed(1)

# SiFINeT+GSVA
gsva_mat <- readRDS("gsva_res_rna.rds")
gsva_mat_temp <- rbind(gsva_mat, rep(0, ncol(gsva_mat)))
cell_group <- apply(gsva_mat_temp[,], 2, which.max)

selected_gene <- c("Ccr7", "Foxo1", "Foxp1", "Tcf7", "Gzma", "Klrd1", "Klrg1")

data <- readRDS("PreprocessedData/CD8_rna_matrix.rds")
data_temp <- data[match(toupper(selected_gene), rownames(data)), ]
data1 <- list()
for (i in 1:5){
  data1[[i]] <- data_temp[, cell_group == i]
}

data_final <- cbind(data1[[1]], data1[[3]], data1[[5]], data1[[2]])
data_final <- as.matrix(data_final)
data_final2 <- data_final
for (i in 1:7){
  data_final2[i, ] <- (data_final[i, ] - min(data_final[i, ]))/diff(range(data_final[i, ]))
}
selected_gene_label <- selected_gene
rownames(data_final2) <- selected_gene_label

data_final1 <- data_final2 %>% 
  as.data.frame() %>%
  rownames_to_column("Gene_name") %>%
  pivot_longer(-c(Gene_name), names_to = "Sample_name", values_to = "counts")
data_final1$Gene_name <- factor(data_final1$Gene_name, levels = selected_gene_label)
data_final1$Sample_name <- factor(data_final1$Sample_name, levels = colnames(data_final2))

group_count <- sapply(data1, ncol)
group_count <- group_count[c(1,3,5,2)]

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
  geom_vline(xintercept = cumsum(group_count))
g1

g2 <- g1 + 
  annotate("text",x = c(0, cumsum(group_count[1:(length(group_count) - 1)])) + group_count / 2, 
           y=-0.05, label = c("TEM_1",
                              "TEM_2",
                              "TCM_2",
                              "TCM_1"), size = 3) + 
  coord_cartesian(ylim=c(1, 7), clip="off")
g2
ggsave("../../Supp_cd8_trans.jpeg", width = 6, height = 2.5, 
       units = "in", device='jpeg', dpi=600)


# Supp_mdata_time_mem
setwd("../../Experimental/MassiveRNA/")

timeres <- read.table("time_res.txt")
timeres <- timeres[, 3]
memres <- read.table("mem_res.txt")
memres <- memres[,2]

time_mem_data <- data.frame(ncell = c(10000, 20000, 50000, 100000, 200000, 400000, 
                                      600000, 800000, 1000000, 1287071),
                            total_time = timeres, peak_memory = memres)

time_mem_data$log_ncell <- log10(time_mem_data$ncell)
time_mem_data$log_total_time <- log10(time_mem_data$total_time)
time_mem_data$log_peak_memory <- log10(time_mem_data$peak_memory)

p1 <- ggplot(time_mem_data, aes(log_ncell, log_total_time)) + 
  geom_point(size = 3) + 
  geom_smooth(method="lm") + 
  labs(x = "log10 number of cells", y = "log10 time cost (s)") + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(), 
        axis.line = element_line(colour = "black"), 
        legend.text=element_text(face="bold"), 
        legend.title=element_blank()) + 
  annotate("text",x=4.5,y=4.5,label=(paste("Slope =",round(coef(lm(time_mem_data$log_total_time~time_mem_data$log_ncell))[2], 3))))

p2 <- ggplot(time_mem_data, aes(log_ncell, log_peak_memory)) + 
  geom_point(size = 3) + 
  geom_smooth(data=subset(time_mem_data, log_ncell >= 5), method="lm") + 
  geom_smooth(data=subset(time_mem_data, log_ncell < 5), method="lm") + 
  labs(x = "log10 number of cells", y = "log10 peak memory usage (MB)") + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(), 
        axis.line = element_line(colour = "black"), 
        legend.text=element_text(face="bold"), 
        legend.title=element_blank()) + 
  annotate("text",x=5.7,y=4,label=(paste("Slope =",round(coef(lm(time_mem_data$log_peak_memory[4:10]~time_mem_data$log_ncell[4:10]))[2], 3)))) + 
  annotate("text",x=4.3,y=4,label=(paste("Slope =",round(coef(lm(time_mem_data$log_peak_memory[1:3]~time_mem_data$log_ncell[1:3]))[2], 3))))

p3 <- ggarrange(p1, p2, nrow = 1,
                labels = "auto")
ggsave("../../Supp_mdata_time_mem.jpeg", width = 8, height = 3, 
       units = "in", device='jpeg', dpi=600)



# Supp_sd1_cluster
setwd("../../Numerical/SD1/data")
cluster_setting <- c("_cluster_seurat.rds", 
                     "_cluster_louvain.rds",
                     "_cluster_cidr.rds")

res <- matrix(0, 100, 3)

for (i in 1:100){
  for (j in 1:3){
    group <- readRDS(paste(i, cluster_setting[j], sep = ""))
    res[i, j] <- adjustedRandIndex(group, rep(1:3, c(100, 4900, 1000)))
  }
}

colnames(res) <- c("Seurat", "Louvain", "CIDR")
res <- data.frame(res)

plotdata <- pivot_longer(res, c("Seurat", "Louvain", "CIDR"))
colnames(plotdata) <- c("Setting", "ARI")
plotdata$Method <- rep(c("Seurat", "Simple Louvain", "CIDR"), 100)
plotdata$Method <- factor(plotdata$Method, levels = c("Seurat", "Simple Louvain", "CIDR"))
plotdata$Setting <- factor(plotdata$Setting, 
                           levels = c("Seurat", "Louvain", "CIDR"), 
                           labels = c("Seurat", "Louvain", "CIDR"))

library(ggplot2)
ggplot(plotdata, aes(x=Setting, y=ARI, fill = Method)) + 
  geom_boxplot() + theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + ylab("ARI") + xlab("") + 
  theme(axis.text=element_text(size=12,face="bold"), 
        axis.title=element_text(size=14,face="bold"),
        legend.text=element_text(size=12,face="bold"), 
        legend.title=element_text(size=14,face="bold"))

ggsave("../../../Supp_sd1_cluster.jpeg", width = 4, height = 3, 
       units = "in", device='jpeg', dpi=600)

# Supp_imm_coex

setwd("../../../Experimental/TvB/")
so <- readRDS("so.rds")

data <- gmt2list("c7.all.v2023.1.Hs.symbols.gmt")
a <- names(data)
b1 <- grepl("TCELL_VS_BCELL_UP", a, fixed = TRUE)
d1 <- unlist(data[b1])
b2 <- grepl("TCELL_VS_BCELL_DN", a, fixed = TRUE)
d2 <- unlist(data[b2])
e1 <- table(d1)
upgenes <- names(e1)[e1 >= 2]
e2 <- table(d2)
dngenes <- names(e2)[e2 >= 2]
rm(a,b1,b2,d1,d2,e1,e2)


group1id <- match(tolower(upgenes), tolower(so@gene.name))
group1id <- group1id[!is.na(group1id)]
group2id <- match(tolower(dngenes), tolower(so@gene.name))
group2id <- group2id[!is.na(group2id)]
group3id <- c(group1id, group2id)
group4id <- setdiff(1:length(so@gene.name), group3id)



coexp <- so@coexp
rm(so)
a <- c(coexp[group1id, group1id][upper.tri(coexp[group1id, group1id])],
       coexp[group2id, group2id][upper.tri(coexp[group2id, group2id])])
b <- c(coexp[group1id, group2id])
d <- c(coexp[group3id, group4id])
e <- c(coexp[group4id, group4id][upper.tri(coexp[group4id, group4id])])

da <- density(a)
db <- density(b)
dd <- density(d)
de <- density(e)


refx <- seq(-5, 5, length.out = length(da$x))
refy <- dnorm(refx)



plotx <- c(da$x, db$x, dd$x, de$x, refx)
ploty <- c(da$y, db$y, dd$y, de$y, refy)
mvec <- c(mean(a), mean(b), mean(d), mean(e), 0)
sdvec <- c(sd(a), sd(b), sd(d), sd(e), 1)
coextype <- c("Same", "Different", "One", "None", "Reference")
lege <- paste(coextype, " (mean: ", round(mvec, 2), ", sd: ", round(sdvec, 2), ")", sep = "")

label <- rep(coextype, 
             c(length(da$x), length(db$x), length(dd$x), length(de$x), length(de$x)))
label <- factor(label, levels = coextype, labels = lege)
plotdata <- data.frame(plotx, ploty, label)
plotdata <- plotdata[abs(plotdata$plotx) <= 10, ]


ggplot(plotdata, aes(plotx, ploty, color = label)) +
  geom_line(size = 2, alpha = 0.5) + theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        legend.title=element_blank()) + 
  labs(x = "Coexpression", y = "Density") + 
  scale_color_manual(values = c("orange", "green", "red", "blue", "gray")) + 
  guides(fill=guide_legend(title="Coexpression type"))

ggsave("../../Supp_imm_coex.jpeg", width = 5, height = 3, 
       units = "in", device='jpeg', dpi=600)




# Supp_thres
setwd("../../Experimental/Monoclonal/")

so <- readRDS("so.rds")
con1 <- so@conn

pathgene <- readRDS("cellcycle/pathgene.rds")
con_path <- data.frame(con1)[(con1$name %in% str_to_title(pathgene)) & (con1$C1 != 0), ]
ccgene <- readRDS("cellcycle/ccgene.rds")
con_ccgene <- data.frame(con1)[(con1$name %in% str_to_title(ccgene)) & (con1$C1 != 0), ]
hubgene <- readRDS("cellcycle/hubgene.rds")
con_hub <- data.frame(con1)[(con1$name %in% str_to_title(hubgene)) & (con1$C1 != 0), ]
con_ccgene$C1 <- log2(con_ccgene$C1 / (nrow(so@coexp) - 1))
con_path$C1 <- log2(con_path$C1 / (nrow(so@coexp) - 1))
con_hub$C1 <- log2(con_hub$C1 / (nrow(so@coexp) - 1))


g1 <- ggplot(data = con_ccgene, aes(x = C2, y = C3)) + theme_bw() + 
  stat_density2d(aes(colour = ..level.., fill = ..level..), geom = "polygon", show.legend = FALSE) + 
  scale_colour_gradient(low = "lightpink", high = "red") +  
  scale_fill_gradient(low = "lightpink", high = "red") + 
  new_scale_color() +
  stat_density2d(data = con_path, aes(x = C2, y = C3, colour = ..level..), show.legend = FALSE) + 
  scale_colour_gradient(low = "lightgreen", high = "darkgreen") + 
  new_scale_color() +
  stat_density2d(data = con_hub, aes(x = C2, y = C3, colour = ..level..), show.legend = FALSE) + 
  scale_colour_gradient(low = "lightblue", high = "blue") + 
  geom_vline(xintercept = 0.4, color = "black") + 
  geom_hline(yintercept = 0.3, color = "black") + 
  labs(x = "C2i", y = "C3i") + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        legend.title=element_blank())



g2 <- ggplot(data = con_ccgene, aes(x = C1, y = C3)) + theme_bw() + 
  stat_density2d(aes(colour = ..level.., fill = ..level..), geom = "polygon", show.legend = FALSE) + 
  scale_colour_gradient(low = "lightpink", high = "red") +  
  scale_fill_gradient(low = "lightpink", high = "red") + 
  new_scale_color() +
  stat_density2d(data = con_path, aes(x = C1, y = C3, colour = ..level..), show.legend = FALSE) + 
  scale_colour_gradient(low = "lightgreen", high = "darkgreen") + 
  new_scale_color() +
  stat_density2d(data = con_hub, aes(x = C1, y = C3, colour = ..level..), show.legend = FALSE) + 
  scale_colour_gradient(low = "lightblue", high = "blue") + 
  geom_vline(xintercept = log2(5 / (nrow(so@coexp) - 1)), color = "black") + 
  geom_hline(yintercept = 0.3, color = "black") + 
  labs(x = "C1i", y = "C3i") + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        legend.title=element_blank())



g3 <- ggplot(data = con_ccgene, aes(x = C1, y = C2)) + theme_bw() + 
  stat_density2d(aes(colour = ..level.., fill = ..level..), geom = "polygon", show.legend = FALSE) + 
  scale_colour_gradient(low = "lightpink", high = "red") +  
  scale_fill_gradient(low = "lightpink", high = "red") + 
  new_scale_color() +
  stat_density2d(data = con_path, aes(x = C1, y = C2, colour = ..level..), show.legend = FALSE) + 
  scale_colour_gradient(low = "lightgreen", high = "darkgreen") + 
  new_scale_color() +
  stat_density2d(data = con_hub, aes(x = C1, y = C2, colour = ..level..), show.legend = FALSE) + 
  scale_colour_gradient(low = "lightblue", high = "blue") + 
  geom_vline(xintercept = log2(5 / (nrow(so@coexp) - 1)), color = "black") + 
  geom_hline(yintercept = 0.4, color = "black") + 
  labs(x = "C1i", y = "C2i") + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        legend.title=element_blank())

g4 <- ggarrange(g1, g2, g3, ncol = 1,
                labels = "auto")


ggsave("../../Supp_thres.jpeg", width = 5, height = 9, 
       units = "in", device='jpeg', dpi=600)
