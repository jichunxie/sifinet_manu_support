# singleCellHaystack
# https://alexisvdb.github.io/singleCellHaystack/articles/examples/a03_example_highD_advanced.html
# https://rdrr.io/cran/singleCellHaystack/src/R/haystack.R
library(singleCellHaystack)

run_singleCellHaystack <- function(data.expression){
  median.per.gene <- apply(data.expression,1,median)
  data.detection <- data.expression > median.per.gene
  general.detection = apply(data.detection, 2, sum)
  
  data.pca <- prcomp(data.expression)
  data.pc <- data.pca$rotation[,1:50]
  
  #set.seed(123)
  res.pc50.adv <- haystack(x = data.pc, detection = data.detection, 
                           method = "highD", use.advanced.sampling = general.detection)
  
  pvalue <- 10^res.pc50.adv$results$log.p.vals
  qvalue <- p.adjust(pvalue, method = "BH")

  out <- data.frame(id = 1:nrow(data.expression), qvalue = qvalue)
  
  return(out)
}

