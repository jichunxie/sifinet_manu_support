# edgeR (QL)
# https://github.com/csoneson/conquer_comparison/blob/master/scripts/apply_edgeRQLF.R
library(edgeR)

run_edgeR <- function(data, group){
  rownames(data) <- 1:nrow(data)
  colnames(data) <- 1:ncol(data)
  d <- DGEList(data)
  d <- calcNormFactors(d)
  group <- factor(group)
  mm <- model.matrix(~group)
  
  d <- estimateDisp(d, design = mm)
  fit <- glmQLFit(d, design = mm)
  ql_t <- glmQLFTest(fit)
  top.table <- topTags(ql_t,sort="none",n=Inf, adjust.method = "BH")
  top.table <- top.table$table
  
  outdata <- data.frame(as.numeric(rownames(top.table)), top.table$FDR)
  return(outdata)
}
