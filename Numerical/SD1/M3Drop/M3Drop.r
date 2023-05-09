source("../../../code/run_M3Drop.R")

for (i in 1:100){
	set.seed(i)
  data1 <- readRDS(paste("../data/", i, "_matrix.rds", sep = ""))

  out <- run_M3Drop(data1)
    
  saveRDS(out, paste(i, "_M3Drop.rds", sep = ""))
}
