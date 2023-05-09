source("../../../code/run_scry.R")

for (i in 1:100){
	set.seed(i)
  data1 <- readRDS(paste("../data/", i, "_matrix.rds", sep = ""))
  
  out <- run_scry(data1)
  
  saveRDS(out, paste(i, "_scry.rds", sep = ""))
}
