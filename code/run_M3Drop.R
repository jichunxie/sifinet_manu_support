# M3Drop
# https://www.bioconductor.org/packages/release/bioc/vignettes/M3Drop/inst/doc/M3Drop_Vignette.pdf
library(M3Drop)

run_M3Drop <- function(data){
  rownames(data) <- 1:nrow(data)
  norm <- M3DropConvertData(data, is.counts=TRUE)
  #norm <- M3DropConvertData(log2(norm+1), is.log=TRUE, pseudocount=1)
  M3Drop_genes <- M3DropFeatureSelection(norm, mt_method="BH", mt_threshold=0.05)
  out <- M3Drop_genes[, c(1,4)]
  out[,1] <- as.numeric(out[,1])
  return(out)
}
