# limma
# https://github.com/csoneson/conquer_comparison/blob/master/scripts/apply_limmatrend.R
library(edgeR)

run_limma <- function(data, group){
  rownames(data) <- 1:nrow(data)
  colnames(data) <- 1:ncol(data)
  d <- DGEList(data)
  d <- calcNormFactors(d)
  group <- factor(group)
  y <- edgeR::cpm(d, log = TRUE, prior.count = 1)
  mm <- model.matrix(~group)
  
  
  fit <- lmFit(y, mm)
  fit <- eBayes(fit)
  top.table <- topTable(fit,sort="none",n=Inf, adjust.method = "BH")
  outdata <- data.frame(as.numeric(rownames(top.table)), top.table$adj.P.Val)
  return(outdata)
}

