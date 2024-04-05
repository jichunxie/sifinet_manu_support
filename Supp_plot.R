library(Seurat)
library(ggplot2)
library(ggpubr)
library(mclust)
library(tidyverse)
library(cogena)
library(Matrix)
library(ggnewscale)
library(stringr)
library(SiFINeT)
library(ggvenn)
library(igraph)
library(ggraph)
library(reshape2)
library(caret)

set.seed(1)

setwd("~/Desktop/SiFINeT/Result_final/")


#Supp fig 2
cal_conn(matrix(c(0,1,1,1,1,0,0,0,
                  1,0,1,1,1,0,0,0,
                  1,1,0,1,1,0,0,0,
                  1,1,1,0,1,1,1,1,
                  1,1,1,1,0,1,1,1,
                  0,0,0,1,1,0,1,1,
                  0,0,0,1,1,1,0,1,
                  0,0,0,1,1,1,1,0), nrow = 8, ncol = 8), thres = 0.5)

#Supp fig 3
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

ggsave("../../../Supp_sd1_cluster.pdf", width = 4, height = 3, 
       units = "in", device='pdf', dpi=1200)

#Supp fig 4
# In plot2.R

#Supp fig 5
setwd("~/Desktop/SiFINeT/Result_final/")
setwd("Numerical/SD1/")
set.seed(1)
data <- readRDS("data/18_matrix.rds")
group1 <- readRDS("data/18_cluster_true.rds")
group2 <- readRDS("data/18_cluster_cidr.rds")
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
ggsave("../../Supp_sd1_umap.pdf", width = 6, height = 2.5, 
       units = "in", device='pdf', dpi=1200)

#Supp fig 6
setwd("../../Experimental/Monoclonal/")
#https://stackoverflow.com/questions/37378744/igraph-grouped-layout-based-on-attribute
set.seed(1)

sio <- readRDS("so.rds")
genename <- sio@gene.name
out <- sio@featureset
edgemat <- abs(sio@coexp - sio@est_ms$mean) >= sio@thres
rm(sio)
set1 <- c(out$unique[[1]], out$shared[[1]], out$enriched[[1]])
set2 <- c(out$unique[[2]], out$shared[[2]], out$enriched[[2]])
set_inter <- intersect(set1, set2)
set1_only <- setdiff(set1, set_inter)
set2_only <- setdiff(set2, set_inter)
idx1 <- match(set1_only, genename)
idx2 <- match(set2_only, genename)
idx3 <- match(set_inter, genename)
idx <- c(idx1, idx2, idx3)

group <- rep(1:3, c(length(idx1), length(idx2), length(idx3)))
edgemat1 <- edgemat[idx, idx]
colnames(edgemat1) <- genename[idx]
rownames(edgemat1) <- genename[idx]
idxlist <- list(idx1, idx2, idx3)
plt <- data.frame(matrix(0, 0, 3))
colnames(plt) <- c("Prediction", "Reference", "Prop")
count <- 1
for (i in 1:3){
  for (j in 1:3){
    if (i != j){
      temp_prop <- sum(edgemat[idxlist[[i]], idxlist[[j]]]) / length(idxlist[[i]]) / length(idxlist[[j]])
    } else {
      temp_prop <- sum(edgemat[idxlist[[i]], idxlist[[j]]]) / length(idxlist[[i]]) / (length(idxlist[[j]]) - 1)
    }
    tempmat <- edgemat1[group == i, group == j]
    tempmat[tempmat != 0] <- temp_prop^5
    edgemat1[group == i, group == j] <- tempmat
    plt[count,] <- c(i, j, temp_prop)
    count <- count + 1
  }
}

plt$Prediction <- factor(plt$Prediction, levels=c(1,3,2))
plt$Reference <- factor(plt$Reference, levels=c(2,3,1))
p5 <- ggplot(plt, aes(Prediction,Reference, fill= Prop)) +
  geom_tile() + geom_text(aes(label=round(Prop, 3))) +
  scale_fill_gradient(low="white", high="#009194") +
  labs(x = "SifiNet feature gene set",y = "SifiNet feature gene set") +
  scale_x_discrete(labels=c("Set 1 unique","Shared","Set 2 unique")) +
  scale_y_discrete(labels=c("Set 2 unique","Shared","Set 1 unique")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5))





G <- graph_from_adjacency_matrix(edgemat1, "undirected", weighted = TRUE)
G_Grouped <- G

for(i in unique(group)) {
  GroupV = which(group == i)
  G_Grouped = add_edges(G_Grouped, combn(GroupV, 2), attr=list(weight=2))
} 
LO <- layout_with_fr(G_Grouped)

pos <- data.frame(LO, genename[idx])
colnames(pos) <- c("x", "y", "name")
pos$group <- factor(group)
G_edge <- get.data.frame(G)
G_edge$start.x <- pos[match(G_edge$from, pos[,3]), 1] 
G_edge$start.y <- pos[match(G_edge$from, pos[,3]), 2]
G_edge$end.x <- pos[match(G_edge$to, pos[,3]), 1]
G_edge$end.y <- pos[match(G_edge$to, pos[,3]), 2]

p1 <- ggplot() +
  geom_segment(data=G_edge,aes(x=start.x,xend = end.x, y=start.y,yend = end.y,size=weight),colour="lightgray", size = 0.00002) +
  geom_point(data=pos,aes(x=x,y=y,col=group),size=2) +  # adds a black border around the nodes
  scale_color_manual(name="Gene set", values = c("#1b9e77", "#d95f02", "#e6ab02"), labels = c("Set 1 unique", "Set 2 unique", "Shared")) + 
  theme_bw()+  # use the ggplot black and white theme
  theme(legend.position = "none",
        axis.text.x = element_blank(),  # remove x-axis text
        axis.text.y = element_blank(), # remove y-axis text
        axis.ticks = element_blank(),  # remove axis ticks
        axis.title.x = element_blank(), # remove x-axis labels
        axis.title.y = element_blank(), # remove y-axis labels
        panel.background = element_blank(), 
        panel.border =element_blank(), 
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank())


ref1 <- sample(1:2, length(idx1), replace = T)
ref2 <- sample(1:2, length(idx2), replace = T)
ref3 <- sample(1:2, length(idx3), replace = T)

geneset_topology2 <- function(id_list, edge_mat, weightthres = 0.3, edge_method = 2, node_color = "black",
                              shiftsize = 0.05, boundsize = 0.1, prefix = "", set_name = NULL){
  if ((edge_method == 1) & (weightthres < 1)){
    weightthres <- 5
  }
  
  nodes <- data.frame(1:length(id_list), sapply(id_list, length))
  if (is.null(set_name)){
    nodes[,1] <- paste(prefix, nodes[,1], sep = "")
  } else{
    nodes[,1] <- set_name
  }
  edge <- data.frame(t(combn(1:length(id_list), 2)))
  edge$edge_weight <- apply(edge, 1, function(x){mean(edge_mat[id_list[[x[1]]], id_list[[x[2]]]])})
  if (is.null(set_name)){
    edge[,1] <- paste(prefix, edge[,1], sep = "")
    edge[,2] <- paste(prefix, edge[,2], sep = "")
  } else{
    edge[,1] <- set_name[edge[,1]]
    edge[,2] <- set_name[edge[,2]]
  }
  edge <- edge[edge$edge_weight > weightthres, ]
  me <- max(edge$edge_weight)
  edge$scaled_edge_weight <- edge$edge_weight / me
  colnames(nodes) <- c("node_id", "gene_count")
  
  g <- graph_from_data_frame(edge, directed = FALSE, nodes)
  
  lay_temp <- create_layout(g, layout = "circle")
  shift <- shiftsize * (max(lay_temp$y) - min(lay_temp$y) + max(lay_temp$x) - min(lay_temp$x))
  bound <- boundsize * (max(lay_temp$y) - min(lay_temp$y) + max(lay_temp$x) - min(lay_temp$x))
  gt <- ggraph(lay_temp) + 
    geom_edge_link(aes(width = scaled_edge_weight))  + 
    scale_edge_width(range=c(0.5,2)) +
    geom_node_point(aes(col = node_color), size = 5, shape = 16) +
    scale_color_manual(name="Gene set", values = node_color, labels = set_name) + 
    theme_void() + 
    theme(legend.text=element_text(size = 12), legend.title=element_text(size = 14)) + 
    xlim(c(min(lay_temp$x) - bound, max(lay_temp$x) + bound)) + 
    ylim(c(min(lay_temp$y) - bound, max(lay_temp$y) + bound)) +
    guides(size = "none", edge_width = "none",
           colour = guide_legend(override.aes = list(size=5)))
  return(gt)
}

p2 <- geneset_topology2(list(idx1[ref1 == 1], idx1[ref1 == 2], 
                             idx2[ref2 == 1], idx2[ref2 == 2],
                             idx3[ref3 == 1], idx3[ref3 == 2]), 
                        edgemat, 
                        node_color = c("#1b9e77", "#1b9e78", 
                                       "#d95f02", "#d95f03", 
                                       "#e6ab02", "#e6ab03"),
                        set_name = c("Set 1 unique 1", "Set 1 unique 2", 
                                     "Set 2 unique 1", "Set 2 unique 2", 
                                     "Shared 1", "Shared 2"))

p3 <- ggarrange(p1, p2, nrow = 1, common.legend = TRUE, 
                legend = "right", legend.grob = get_legend(p1), 
                labels = c("a", "b"))




s.gene <- cc.genes$s.genes
g2m.gene <- cc.genes$g2m.genes

idx1 <- match(s.gene, toupper(genename))
idx2 <- match(g2m.gene, toupper(genename))
idx1 <- idx1[!is.na(idx1)]
idx2 <- idx2[!is.na(idx2)]
idx <- c(idx1, idx2)

group <- rep(1:2, c(length(idx1), length(idx2)))
edgemat1 <- edgemat[idx, idx]
colnames(edgemat1) <- genename[idx]
rownames(edgemat1) <- genename[idx]

G <- graph_from_adjacency_matrix(edgemat1, "undirected", weighted = TRUE)
LO <- layout_with_fr(G)
pos <- data.frame(LO, genename[idx])
colnames(pos) <- c("x", "y", "name")
pos$group <- factor(group)
G_edge <- get.data.frame(G)
G_edge$start.x <- pos[match(G_edge$from, pos[,3]), 1] 
G_edge$start.y <- pos[match(G_edge$from, pos[,3]), 2]
G_edge$end.x <- pos[match(G_edge$to, pos[,3]), 1]
G_edge$end.y <- pos[match(G_edge$to, pos[,3]), 2]

p4 <- ggplot() +
  geom_segment(data=G_edge,aes(x=start.x,xend = end.x, y=start.y,yend = end.y,size=weight),colour="lightgray", size = 0.02) +
  geom_point(data=pos,aes(x=x,y=y,col=group),size=2) +  # adds a black border around the nodes
  scale_color_manual(name="Seurat cell\ncycle genes", values = c("#1b9e77", "#d95f02"), labels = c("G1/S", "G2/M")) + 
  theme_bw() + 
  theme(axis.text.x = element_blank(),  # remove x-axis text
        axis.text.y = element_blank(), # remove y-axis text
        axis.ticks = element_blank(),  # remove axis ticks
        axis.title.x = element_blank(), # remove x-axis labels
        axis.title.y = element_blank(), # remove y-axis labels
        panel.background = element_blank(), 
        panel.border =element_blank(), 
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank())


plt <- data.frame(matrix(0, 0, 3))
colnames(plt) <- c("Prediction", "Reference", "Prop")
count <- 1
for (i in 1:2){
  for (j in 1:2){
    if (i != j){
      temp_prop <- sum(edgemat1[(group == i), (group == j)]) / sum((group == i)) / sum((group == j))
    } else {
      temp_prop <- sum(edgemat1[(group == i), (group == j)]) / sum((group == i)) / (sum((group == j)) - 1)
    }
    plt[count,] <- c(i, j, temp_prop)
    count <- count + 1
  }
}

plt$Prediction <- factor(plt$Prediction, levels=c(1,2))
plt$Reference <- factor(plt$Reference, levels=c(2,1))
p6 <- ggplot(plt, aes(Prediction,Reference, fill= Prop)) +
  geom_tile() + geom_text(aes(label=round(Prop, 3))) +
  scale_fill_gradient(low="white", high="#009194") +
  labs(x = "Seurat cell cycle genes",y = "Seurat cell cycle genes") +
  scale_x_discrete(labels=c("G1/S", "G2/M")) +
  scale_y_discrete(labels=c("G2/M", "G1/S")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5))

p7 <- ggarrange(p4, p5, p6, nrow = 1, 
                labels = c("c", "d", "e"))

p8 <- ggarrange(p3, p7, nrow = 2)

ggsave("../../Supp_mono_net.pdf", width = 10, height = 6, 
       units = "in", device='pdf', dpi=1200)

#Supp fig 7
setwd("../../Experimental/Monoclonal/")

set.seed(1)
genenames <- readRDS("genename.rds")
ccgene <- readRDS("cellcycle/ccgene.rds")
ccid <- match(ccgene, tolower(genenames))
ccgene <- ccgene[!is.na(ccid)]
ccid <- ccid[!is.na(ccid)]


sio <- readRDS("so.rds")
out <- sio@featureset
rm(sio)
sifi_res <- unique(unlist(out))
sifi_id <- match(sifi_res, genenames)
out <- readRDS("DESeq2/DESeq2.rds")
deseq_id <- which(out[,2] <= 0.05)
deseq_res <- genenames[deseq_id]

plotdata <- list(ccgene, tolower(sifi_res), tolower(deseq_res))
names(plotdata) <- c("Benchmark", "SifiNet", "CIDR-DESeq2-TS")
ggvenn(plotdata)
ggsave("../../Supp_mono_venn.pdf", width = 4, height = 5, 
       units = "in", device='pdf', dpi=1200)

#Supp fig 8
# In plot3.R

#Supp fig 9
# In plot3.R

#Supp fig 10
setwd("Experimental/CD8")
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
ggsave("../../Supp_cd8_trans.pdf", width = 6, height = 2.5, 
       units = "in", device='pdf', dpi=1200)

#Supp fig 11
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
ggsave("../../Supp_mdata_time_mem.pdf", width = 8, height = 3, 
       units = "in", device='pdf', dpi=1200)


#Note fig 1
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

ggsave("../../Supp_imm_coex.pdf", width = 5, height = 3, 
       units = "in", device='pdf', dpi=1200)

#Note fig 2
set.seed(1)
gene.means <- rgamma(1000, shape = 81, rate = 9)
gene.means.fc <- matrix(gene.means, nrow = 1000, ncol = 1000, byrow = FALSE)
Y <- matrix(rpois(1000000, lambda = gene.means.fc), 
            nrow = 1000, ncol = 1000)
ncell <- c(100, 1000)
out <- list()
for (i in 1:2){
  so <- create_SiFINeT_object(counts = Y[, 1:ncell[i]], gene.name = 1:nrow(Y))
  so <- quantile_thres(so)
  so <- feature_coexp(so)
  coexp <- so@coexp
  coexp_value <- coexp[upper.tri(coexp)]
  out[[i]] <- coexp_value
}

x <- seq(-5, 5, length.out = 1000)
yn <- dnorm(x)

plotdata <- data.frame(value = c(out[[1]], out[[2]]), 
                       N = factor(rep(c(100, 1000), 
                                      c(length(out[[1]]), length(out[[2]]))), 
                                  levels = c(100, 1000)))
plotdata2 <- data.frame(density = yn, value = x)

ggplot() +
  stat_density(data = plotdata, aes(x=value, group=N, col = N), geom="line", 
               position = "identity", linewidth = 1.2) + 
  geom_line(data = plotdata2, aes(x=value, y=density), alpha = 0.5, linewidth = 0.5) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + 
  guides(colour = guide_legend(override.aes=list(size=1)))

ggsave("../../Supp_thm1.pdf", width = 5, height = 4, 
       units = "in", device='pdf', dpi=1200)

#Note fig 3
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


ggsave("../../Supp_thres.pdf", width = 5, height = 9, 
       units = "in", device='pdf', dpi=1200)

#Note fig 4
a <- matrix(0.0, nrow = 9, ncol = 9)
a[1, ] <- c(0.8, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
a[2, ] <- c(0.6, 0.6, 0.5, 0.4, 0.0, 0.0, 0.0, 0.0, 0.0)
a[3, ] <- c(0.0, 0.8, 0.8, 0.8, 0.0, 0.0, 0.0, 0.0, 0.0)
a[4, ] <- c(0.0, 0.4, 0.5, 0.6, 0.6, 0.0, 0.0, 0.0, 0.0)
a[5, ] <- c(0.0, 0.0, 0.0, 0.0, 0.8, 0.0, 0.0, 0.0, 0.0)
a[6, ] <- c(0.0, 0.0, 0.0, 0.0, 0.6, 0.6, 0.5, 0.4, 0.0)
a[7, ] <- c(0.0, 0.0, 0.0, 0.0, 0.0, 0.8, 0.8, 0.8, 0.0)
a[8, ] <- c(0.0, 0.0, 0.0, 0.0, 0.0, 0.4, 0.5, 0.6, 0.6)
a[9, ] <- c(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.8)

b <- a
b[b != 0] <- exp(b[b != 0] + 0.02)
b[b == 0] <- 1.0
b <- cbind(c("G1", "G2", "G3", "G4", 
             "G5", "G6", "G7", "G8", "G9"), b)
b <- data.frame(b)
colnames(b) <- c("Geneset", "1", "2a", "2b", "2c", 
                 "3", "4a", "4b", "4c", "5")
plotdata <- melt(b, id.vars = c("Geneset"), variable.name = "Cellpop")
plotdata$value <- as.numeric(plotdata$value)
plotdata$Geneset <- as.factor(plotdata$Geneset)

ggplot(plotdata, aes(Cellpop, Geneset)) +
  geom_tile(aes(fill = value)) +
  scale_fill_gradient(low = "white", high = "red") + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme(axis.ticks.x=element_blank()) + 
  geom_vline(xintercept = c(1.5, 4.5, 5.5, 8.5)) + 
  labs(x = "Cell type", y = "Feature gene set", fill = "Mean fold change")

ggsave("../../Supp_sd2_setting.pdf", width = 6, height = 2.5, 
       units = "in", device='pdf', dpi=1200)




