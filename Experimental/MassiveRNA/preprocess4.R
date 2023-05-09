library(Matrix)

Z <- c()
datalist <- list()
for (i in 0:9){
  print(i)
  data <- readMM(paste("matrix", i, "_r.mtx", sep = ""))
  datalist[[i+1]] <- data
  Z <- c(Z, rowMeans(data))
}
saveRDS(Z, "Z.rds")
p <- ncol(data)
npd <- ceiling(p / 10)


for (i in 0:9){
  print(i)
  data <- rbind(datalist[[1]][, (i*npd + 1):min(p, (i+1) * npd)],
                datalist[[2]][, (i*npd + 1):min(p, (i+1) * npd)],
                datalist[[3]][, (i*npd + 1):min(p, (i+1) * npd)],
                datalist[[4]][, (i*npd + 1):min(p, (i+1) * npd)],
                datalist[[5]][, (i*npd + 1):min(p, (i+1) * npd)],
                datalist[[6]][, (i*npd + 1):min(p, (i+1) * npd)],
                datalist[[7]][, (i*npd + 1):min(p, (i+1) * npd)],
                datalist[[8]][, (i*npd + 1):min(p, (i+1) * npd)],
                datalist[[9]][, (i*npd + 1):min(p, (i+1) * npd)],
                datalist[[10]][, (i*npd + 1):min(p, (i+1) * npd)])
  writeMM(data, paste("matrix", i, "_v.mtx", sep = ""))
}

rm(list = ls())
