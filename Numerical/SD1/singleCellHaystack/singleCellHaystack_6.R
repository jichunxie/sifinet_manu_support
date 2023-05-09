source("../../../code/run_singleCellHaystack.R")

for (i in 51:60){
  set.seed(i)
	data1 <- readRDS(paste("../data/", i, "_matrix.rds", sep = ""))
  
  out <- run_singleCellHaystack(data1)
  
  saveRDS(out, paste(i, "_singleCellHaystack.rds", sep = ""))
}

rm(list = ls())
