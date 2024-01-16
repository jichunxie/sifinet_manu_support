library(Seurat)
library(mclust)
library(grid)
library(ggpubr)
library(ggplot2)
library(magick)
library(tiff)

setwd("~/Desktop/SiFINeT/Result_final")
setwd("Numerical/SD1/")

set.seed(1)

res1_sen <- matrix(0, 4, 100)
res1_spe <- matrix(0, 4, 100)
res1_f1s <- matrix(0, 4, 100)
res2_sen <- matrix(0, 5, 400)
res2_spe <- matrix(0, 5, 400)
res2_f1s <- matrix(0, 5, 400)
meth1 <- c("SifiNet", "M3Drop", "scry", "singleCellHaystack")
meth2 <- c("DESeq2", "edgeR", "limma", "voom_limma", "MAST")
setting_list <- c("true", "seurat", "louvain", "cidr")
for (j in 1:100){
  gene_name <- readRDS(paste("data/", j, "_name.rds", sep = ""))
  for (i in 1:4){
    if(i == 1){
      fset <- readRDS(paste("data/featureset/", j, "_featureset.rds", sep = ""))
      detection <- unlist(fset)
    } else if (i == 3){
      res <- readRDS(paste(meth1[i], "/", j, "_", meth1[i], ".rds", sep = ""))
      t_gene <- sum(gene_name <= 120)
      ord <- order(res, decreasing = T)
      detection <- gene_name[ord[1:t_gene]]
    } else {
      res <- readRDS(paste(meth1[i], "/", j, "_", meth1[i], ".rds", sep = ""))
      detection <- gene_name[res[res[,2] <= 0.05, 1]]
    }
    tp <- sum(detection <= 120)
    fp <- sum(detection > 120)
    fn <- sum(gene_name <= 120) - tp
    tn <- sum(gene_name > 120) - fp
    res1_sen[i, j] <- tp / (tp + fn)
    res1_spe[i, j] <- tn / (tn + fp)
    res1_f1s[i, j] <- 2 * tp / (2 * tp + fp + fn)
  }
}

for (k in 1:4){
  for (j in 1:100){
    clu <- readRDS(paste("data/", j, "_cluster_", setting_list[k], ".rds", sep = ""))
    gene_name <- readRDS(paste("data/", j, "_name.rds", sep = ""))
    for (i in 1:5){
      res <- readRDS(paste(meth2[i], "/", j, "_", meth2[i], "_", setting_list[k], ".rds", sep = ""))
      detection <- gene_name[res[res[,2] <= 0.05, 1]]
      detection <- detection[!is.na(detection)]
      tp <- sum(detection <= 120)
      fp <- sum(detection > 120)
      fn <- sum(gene_name <= 120) - tp
      tn <- sum(gene_name > 120) - fp
      res2_sen[i, (k - 1) * 100 + j] <- tp / (tp + fn)
      res2_spe[i, (k - 1) * 100 + j] <- tn / (tn + fp)
      res2_f1s[i, (k - 1) * 100 + j] <- 2 * tp / (2 * tp + fp + fn)
    }
  }
}
apply(res1_f1s, 1, median)

a <- c()
b <- c()
for (j in 1:100){
  gene_name <- readRDS(paste("data/", j, "_name.rds", sep = ""))
  a[j] <- sum(gene_name <= 120)
  b[j] <- length(gene_name)
}
mean(a)
mean(b) - mean(a)

Method <- rep(c("SifiNet", "M3Drop", "scry", "singleCellHaystack",
                "DESeq2-T", "edgeR-T", "limma-T", "limma-voom-T", "MAST-T",
                "DESeq2-Seurat-D", "edgeR-Seurat-D", "limma-Seurat-D", "limma-voom-Seurat-D", "MAST-Seurat-D",
                "DESeq2-Louvain-D", "edgeR-Louvain-D", "limma-Louvain-D", "limma-voom-Louvain-D", "MAST-Louvain-D",
                "DESeq2-CIDR-D", "edgeR-CIDR-D", "limma-CIDR-D", "limma-voom-CIDR-D", "MAST-CIDR-D"), each = 100)
sen_value <- c(res1_sen[1, ], res1_sen[2, ], res1_sen[3, ], res1_sen[4, ], 
               res2_sen[1, 1:100], res2_sen[2, 1:100], res2_sen[3, 1:100], res2_sen[4, 1:100], res2_sen[5, 1:100],
               res2_sen[1, 101:200], res2_sen[2, 101:200], res2_sen[3, 101:200], res2_sen[4, 101:200], res2_sen[5, 101:200],
               res2_sen[1, 201:300], res2_sen[2, 201:300], res2_sen[3, 201:300], res2_sen[4, 201:300], res2_sen[5, 201:300],
               res2_sen[1, 301:400], res2_sen[2, 301:400], res2_sen[3, 301:400], res2_sen[4, 301:400], res2_sen[5, 301:400])

spe_value <- c(res1_spe[1, ], res1_spe[2, ], res1_spe[3, ], res1_spe[4, ], 
               res2_spe[1, 1:100], res2_spe[2, 1:100], res2_spe[3, 1:100], res2_spe[4, 1:100], res2_spe[5, 1:100],
               res2_spe[1, 101:200], res2_spe[2, 101:200], res2_spe[3, 101:200], res2_spe[4, 101:200], res2_spe[5, 101:200],
               res2_spe[1, 201:300], res2_spe[2, 201:300], res2_spe[3, 201:300], res2_spe[4, 201:300], res2_spe[5, 201:300],
               res2_spe[1, 301:400], res2_spe[2, 301:400], res2_spe[3, 301:400], res2_spe[4, 301:400], res2_spe[5, 301:400])

f1s_value <- c(res1_f1s[1, ], res1_f1s[2, ], res1_f1s[3, ], res1_f1s[4, ], 
               res2_f1s[1, 1:100], res2_f1s[2, 1:100], res2_f1s[3, 1:100], res2_f1s[4, 1:100], res2_f1s[5, 1:100], 
               res2_f1s[1, 101:200], res2_f1s[2, 101:200], res2_f1s[3, 101:200], res2_f1s[4, 101:200], res2_f1s[5, 101:200],
               res2_f1s[1, 201:300], res2_f1s[2, 201:300], res2_f1s[3, 201:300], res2_f1s[4, 201:300], res2_f1s[5, 201:300], 
               res2_f1s[1, 301:400], res2_f1s[2, 301:400], res2_f1s[3, 301:400], res2_f1s[4, 301:400], res2_f1s[5, 301:400])

Method_type <- rep(c("SifiNet", "Cluster-free methods", "True cell clusters", 
                     "Default resolution"), c(100, 300, 500, 1500))

plotdata <- data.frame(Method, sen_value, spe_value, f1s_value, Method_type)
plotdata$Method <- factor(plotdata$Method, 
                          levels = c("DESeq2-T", "edgeR-T", "limma-T", "limma-voom-T", "MAST-T", 
                                     "SifiNet", "M3Drop", "scry", "singleCellHaystack",
                                     "DESeq2-Seurat-D", "edgeR-Seurat-D", "limma-Seurat-D", "limma-voom-Seurat-D", "MAST-Seurat-D",
                                     "DESeq2-Louvain-D", "edgeR-Louvain-D", "limma-Louvain-D", "limma-voom-Louvain-D", "MAST-Louvain-D",
                                     "DESeq2-CIDR-D", "edgeR-CIDR-D", "limma-CIDR-D", "limma-voom-CIDR-D", "MAST-CIDR-D"), 
                          labels = c("DESeq2-OC", "edgeR-OC", "limma-OC", "limma-voom-OC", "MAST-OC", 
                                     "SiftNet", "M3Drop", "scry", "singleCellHaystack",
                                     "DESeq2-Seurat-TS", "edgeR-Seurat-TS", "limma-Seurat-TS", "limma-voom-Seurat-TS", "MAST-Seurat-TS",
                                     "DESeq2-Louvain-TS", "edgeR-Louvain-TS", "limma-Louvain-TS", "limma-voom-Louvain-TS", "MAST-Louvain-TS",
                                     "DESeq2-CIDR-TS", "edgeR-CIDR-TS", "limma-CIDR-TS", "limma-voom-CIDR-TS", "MAST-CIDR-TS"))

plotdata$Method_type <- factor(plotdata$Method_type, levels = c("SifiNet", "Cluster-free methods", "True cell clusters", 
                                                                "Default resolution"))


wide <- plotdata
colnames(wide) <- c("Method", "Sensitivity", "Specificity", "F1 score", "Method_type")
long <- tidyr::pivot_longer(
  wide, c("Sensitivity", "Specificity", "F1 score"), 
  values_to = "value", names_to = "variable"
)
long$variable <- factor(long$variable, levels = c("Sensitivity", "Specificity", "F1 score"))
a <- ifelse(c("DESeq2-OC", "edgeR-OC", "limma-OC", "limma-voom-OC", "MAST-OC", 
              "SifiNet", "M3Drop", "scry", "singleCellHaystack",
              "DESeq2-Seurat-TS", "edgeR-Seurat-TS", "limma-Seurat-TS", "limma-voom-Seurat-TS", "MAST-Seurat-TS", 
              "DESeq2-Louvain-TS", "edgeR-Louvain-TS", "limma-Louvain-TS", "limma-voom-Louvain-TS", "MAST-Louvain-TS",
              "DESeq2-CIDR-TS", "edgeR-CIDR-TS", "limma-CIDR-TS", "limma-voom-CIDR-TS", "MAST-CIDR-TS") == "SifiNet", "red", "black")
p1 <- ggplot(long, aes(x=Method, y=value, fill=Method_type)) + 
  facet_wrap(~ variable) + 
  geom_boxplot() + # geom_dotplot(binaxis='y', stackdir='center', dotsize=0.2) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.title.y=element_blank(), 
        axis.title.x=element_blank()) +
  theme(axis.text=element_text(size=10,face="bold"), 
        axis.text.y=element_text(color=rev(a)),
        axis.title=element_text(size=12,face="bold")) + 
  geom_vline(aes(xintercept=15.5), colour="#990000")+ 
  geom_vline(aes(xintercept=18.5), colour="#990000") + 
  geom_vline(aes(xintercept=19.5), colour="#990000") + 
  theme(legend.position = "none") + coord_flip() +
  scale_x_discrete(limits=rev)
p1



tiff("../../Supp_sd1_res.tiff", width = 8, height = 8, units = "in", res = 1200)
p1 + theme(plot.margin = unit(c(0, 0.15, 0, 0), "npc")) 
grid.text(unit(0.86, "npc"), unit(0.75, "npc"), 
          label=expression(paste('SifiNet'),'type=4'), 
          hjust = 0, vjust=0, gp=gpar(fontsize=7, fontface = "bold"))
grid.text(unit(0.86, "npc"), unit(0.67, "npc"), 
          label=expression(paste('Cluster-independent\napproach'),'type=4'), 
          hjust = 0, vjust=0, gp=gpar(fontsize=7, fontface = "bold"))
grid.text(unit(0.86, "npc"), unit(0.86, "npc"), 
          label=expression(paste('Cluster-based methods\nOracle cluster'),'type=4'), 
          hjust = 0, vjust=0, gp=gpar(fontsize=7, fontface = "bold"))
grid.text(unit(0.86, "npc"), unit(0.32, "npc"), 
          label=expression(paste('Cluster-based methods\nTwo step'),'type=4'), 
          hjust = 0, vjust=0, gp=gpar(fontsize=7, fontface = "bold"))
dev.off()


pdf("../../Supp_sd1_res.pdf")
plot(image_read("../../Supp_sd1_res.tiff"))
dev.off()



Method <- rep(c("SifiNet", "M3Drop", "scry", "singleCellHaystack",
                "DESeq2", "edgeR", "limma", "limma-voom", "MAST"), each = 100)

sen_value <- c(res1_sen[1, ], res1_sen[2, ], res1_sen[3, ], res1_sen[4, ], 
               res2_sen[1, 301:400], res2_sen[2, 301:400], res2_sen[3, 301:400], res2_sen[4, 301:400], res2_sen[5, 301:400])

spe_value <- c(res1_spe[1, ], res1_spe[2, ], res1_spe[3, ], res1_spe[4, ], 
               res2_spe[1, 301:400], res2_spe[2, 301:400], res2_spe[3, 301:400], res2_spe[4, 301:400], res2_spe[5, 301:400])

f1s_value <- c(res1_f1s[1, ], res1_f1s[2, ], res1_f1s[3, ], res1_f1s[4, ], 
               res2_f1s[1, 301:400], res2_f1s[2, 301:400], res2_f1s[3, 301:400], res2_f1s[4, 301:400], res2_f1s[5, 301:400])


Method_type <- rep(c("SifiNet", "Cluster-free methods", "Cluster-based methods"), 
                   c(100, 300, 500))

plotdata <- data.frame(Method, sen_value, spe_value, f1s_value, Method_type)
plotdata$Method <- factor(plotdata$Method, 
                          levels = c("SifiNet", "M3Drop", "scry", "singleCellHaystack",
                                     "DESeq2", "edgeR", "limma", "limma-voom", "MAST"), 
                          labels = c("SifiNet", "M3Drop", "scry", "singleCellHaystack",
                                     "DESeq2-CIDR-TS", "edgeR-CIDR-TS", "limma-CIDR-TS", 
                                     "limma-voom-CIDR-TS", "MAST-CIDR-TS"))

plotdata$Method_type <- factor(plotdata$Method_type, levels = c("SifiNet", "Cluster-free methods", "Cluster-based methods"))
wide <- plotdata
colnames(wide) <- c("Method", "Sensitivity", "Specificity", "F1 score", "Method_type")
long <- tidyr::pivot_longer(
  wide, c("Sensitivity", "Specificity", "F1 score"), 
  values_to = "value", names_to = "variable"
)
long$variable <- factor(long$variable, levels = c("Sensitivity", "Specificity", "F1 score"))

p2 <- ggplot(long, aes(x=Method, y=value, fill=Method_type)) + 
  facet_wrap(~ variable) + 
  geom_boxplot() + # geom_dotplot(binaxis='y', stackdir='center', dotsize=0.2) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.title.y=element_blank(), 
        axis.title.x=element_blank()) + #ylab("") + xlab("") + 
  theme(axis.text=element_text(size=10,face="bold"), 
        axis.title=element_text(size=12,face="bold"),
        legend.text=element_text(size=10,face="bold"), 
        legend.title=element_text(size=12,face="bold")) + 
  geom_vline(aes(xintercept=5.5), colour="#990000")+ 
  geom_vline(aes(xintercept=8.5), colour="#990000") + 
  theme(legend.position = "none") + coord_flip() +
  scale_x_discrete(limits=rev)
p2



setwd("../../Experimental/Monoclonal/")

set.seed(1)
genenames <- readRDS("genename.rds")
ccgene <- readRDS("cellcycle/ccgene.rds")
ccid <- match(ccgene, tolower(genenames))
ccgene <- ccgene[!is.na(ccid)]
ccid <- ccid[!is.na(ccid)]

res <- data.frame(matrix(0, 8, 3))

methods <- c("SifiNet", "M3Drop", "singleCellHaystack", "DESeq2", 
             "edgeR", "limma", "voom_limma", "MAST")
clu <- readRDS("cluster.rds")
for (i in 1:8){
  if(i == 1){
    sio <- readRDS("so.rds")
    out <- sio@featureset
    rm(sio)
    detection <- match(unique(unlist(out)), genenames)
  } else {
    out <- readRDS(paste(methods[i], "/", methods[i], ".rds", sep = ""))
    detection <- which(out[,2] <= 0.05)
  }
  tp <- length(intersect(tolower(detection), ccid))
  fp <- length(setdiff(tolower(detection), ccid))
  fn <- length(setdiff(ccid, tolower(detection)))
  tn <- length(genenames) - tp - fp - fn
  res[i, 1] <- tp / (tp + fn)
  res[i, 2] <- tn / (tn + fp)
  res[i, 3] <- 2 * tp / (2 * tp + fp + fn)
}


colnames(res) <- c("Sensitivity", "Specificity", "F1 score")
res$Method <- c("SifiNet", "M3Drop", "singleCellHaystack", "DESeq2", "edgeR", 
                "limma", "limma-voom", "MAST")
res$Method_type <- rep(c("SifiNet", "Cluster-free methods", "Cluster-based methods"), 
                       c(1, 2, 5))
res$Method <- factor(res$Method, 
                     levels = c("SifiNet", "M3Drop", "singleCellHaystack", "DESeq2", 
                                "edgeR", "limma", "limma-voom", "MAST"),
                     labels = c("SifiNet", "M3Drop", "singleCellHaystack", "DESeq2-CIDR-TS", 
                                "edgeR-CIDR-TS", "limma-CIDR-TS", "limma-voom-CIDR-TS", "MAST-CIDR-TS"))
res$Method_type <- factor(res$Method_type, levels = c("SifiNet", "Cluster-free methods", "Cluster-based methods"))


wide <- res
long <- tidyr::pivot_longer(
  wide, c("Sensitivity", "Specificity", "F1 score"), 
  values_to = "value", names_to = "variable"
)
long$variable <- factor(long$variable, levels = c("Sensitivity", "Specificity", "F1 score"))


p3 <- ggplot(long, aes(x=Method, y=value, fill=Method_type)) + 
  facet_wrap(~ variable) + 
  geom_bar(stat="identity") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.title.y=element_blank(), 
        axis.title.x=element_blank()) + 
  theme(axis.text=element_text(size=10,face="bold"), 
        axis.title=element_text(size=12,face="bold"),
        legend.text=element_text(size=10,face="bold"), 
        legend.title=element_text(size=12,face="bold")) + 
  geom_vline(aes(xintercept=5.5), colour="#990000")+ 
  geom_vline(aes(xintercept=7.5), colour="#990000") + 
  theme(legend.position = "none") + coord_flip() +
  scale_x_discrete(limits=rev)
p3

setwd("../../Experimental/IPF/")

genename <- readRDS("2_genename.rds")
testres <- readRDS("2_test.rds")
so <- readRDS("so.rds")
gs <- so@featureset
rm(so)

gene_list <- list()
for (i in 1:length(gs$unique)){
  gene_list[[i]] <- c(gs$unique[[i]], gs$shared[[i]], gs$enriched[[i]])
}

plotdata <- data.frame(gene = unlist(gene_list), 
                       Geneset = rep(paste("Set", 1:length(gene_list), sep = ""), 
                                     sapply(gene_list, length)))
plotdata <- plotdata[!duplicated(plotdata$gene), ]

plotdata$lfc <- testres$avg_log2FC[match(plotdata$gene, testres$gene)]
plotdata$Geneset <- factor(plotdata$Geneset, levels = paste("Set", 1:length(gene_list), sep = ""))
plotdata$gene <- factor(plotdata$gene, levels = plotdata$gene)

plotdata$pvalue <- testres$p_val_adj[match(plotdata$gene, testres$gene)]
plotdata$nlpvalue <- -log10(plotdata$pvalue)
plotdata$sig <- ifelse(floor(plotdata$nlpvalue) >= 6, 3, floor(plotdata$nlpvalue / 2))
plotdata$sig[is.na(plotdata$sig)] <- 0
#https://stackoverflow.com/questions/69936965/how-to-add-an-asterix-significance-above-specific-bars-in-faceted-bar-graph-g
p4 <- ggplot(plotdata, aes(fill=Geneset, y=lfc, x=gene)) + 
  geom_bar(position="dodge", stat="identity") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + 
  scale_fill_manual(values=c("#E69F00", "#56B4E9", "#009E73")) + 
  theme(axis.text=element_text(size=6,face="bold"), 
        axis.title=element_text(size=12,face="bold"),
        legend.text=element_text(size=10,face="bold"), 
        legend.title=element_text(size=12,face="bold"),
        axis.title.x=element_blank(),
        legend.position="left") + 
  geom_vline(xintercept = 50.5, linewidth = .2) + 
  geom_vline(xintercept = 74.5, linewidth = .2) +
  geom_text(aes(label = ifelse(sig == 3, "***", ifelse(sig == 2, "**", ifelse(sig == 1, "*", "")))), 
            position = position_dodge(width = .9), vjust = 0.5 - sign(plotdata$lfc) * 0.5, size = 8 / .pt) + 
  ylab("Log2 fold change")
p4



p6 <- ggarrange(p2, p3, p4, nrow = 3,
                labels = "auto")
p6
ggsave("../../Fig2.pdf", width = 8, height = 8, 
       units = "in", device='pdf', dpi=1200)
