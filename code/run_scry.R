# scry
# https://www.bioconductor.org/packages/release/bioc/vignettes/scry/inst/doc/scry.html
library(scry)

run_scry <- function(data){
  sce<-devianceFeatureSelection(data)
  return(sce)
}
