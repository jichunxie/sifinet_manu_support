# DESeq2
# https://github.com/csoneson/conquer_comparison/blob/master/scripts/apply_DESeq2nofilt.R
# https://github.com/csoneson/conquer_comparison/blob/master/scripts/apply_DESeq2devel.R
library(DESeq2)

run_DESeq2 <- function(data, group){
  rownames(data) <- 1:nrow(data)
  colnames(data) <- 1:ncol(data)
  data[data == 0] <- 1
  group <- factor(group)
  coldata <- data.frame(group = group)
  
  dds <- DESeqDataSetFromMatrix(countData=data, colData=coldata, design=~group)
  
  dds <- DESeq(dds)
  
  res <- results(dds, pAdjustMethod = "BH", alpha = 0.05)
  
  outdata <- data.frame(as.numeric(rownames(res)), res$padj)
  return(outdata)
}