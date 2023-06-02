library(GSVA)
library(Seurat)
library(readxl)

#setwd("~/Desktop/SiFINeT/Result/Experimental/IPF")
set.seed(1)

gene_list <- list()
data <- readRDS("1_matrix.rds")
genename <- readRDS("1_genename.rds")
rownames(data) <- genename
colnames(data) <- 1:ncol(data)
so <- CreateSeuratObject(data)
so <- SCTransform(so, return.only.var.genes = FALSE)
Y <- so@assays$SCT@data

SenMayo <- read_excel("41467_2022_32552_MOESM4_ESM.xlsx")
geneset <- SenMayo$`Gene(human)`

so <- readRDS("so.rds")
sifiset <- tolower(unique(unlist(so@featureset)))
geneset1 <- setdiff(tolower(geneset), sifiset)

temp <- match(tolower(geneset1), tolower(genename))
gene_list[[1]] <- temp[!is.na(temp)]

outfile <- "gsva_res_SenMayo.rds"

Y <- as.matrix(Y)
rownames(Y) <-  1:nrow(Y)
colnames(Y) <- paste0("s", 1:ncol(Y))


gsva.es <- gsva(Y, gene_list, parallel.sz = 8)

gsva.es.mat <- as.matrix(gsva.es)

saveRDS(gsva.es.mat, outfile)
