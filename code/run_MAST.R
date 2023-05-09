# MAST
# https://github.com/csoneson/conquer_comparison/blob/master/scripts/apply_MASTcpm.R
library(MAST)
library(edgeR)

run_MAST <- function(data, group){
  rownames(data) <- 1:nrow(data)
  colnames(data) <- 1:ncol(data)
  group <- factor(group)
  dge <- DGEList(counts = data)
  dge <- edgeR::calcNormFactors(dge)
  cpms <- edgeR::cpm(dge)
  sca <- FromMatrix(exprsArray = log2(cpms + 1), 
                    cData = data.frame(group = group))
  
  zlmdata <- zlm(~group, sca)
  mast <- lrTest(zlmdata, "group")
  temp <- mast[, "hurdle", "Pr(>Chisq)"]
  temp <- p.adjust(temp, method = "BH")
  outdata <- data.frame(as.numeric(names(temp)), temp)
  return(outdata)
}