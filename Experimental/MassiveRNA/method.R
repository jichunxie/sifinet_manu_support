library(Matrix)
library(SiFINeT)
set.seed(1)
tt <- sum(.Internal(gc(FALSE, TRUE, TRUE))[13:14])

# step 1
Z <- readRDS("Z.rds")
Z <- matrix(Z, ncol = 1)

q5 <- c()

for (i in 0:9){
  data <- readMM(paste("matrix", i, "_v.mtx", sep = ""))
  #data <- as(data, "dgCMatrix")
  q5 <- c(q5, apply(data, 2, median))
  so <- create_SiFINeT_object(counts = data, sparse = TRUE, meta.data = Z, rowfeature = F)
  so <- quantile_thres(so)
  writeMM(so@data.thres$data1, paste("rqres", i, ".mtx", sep = ""))
  rm(so)
}
rm(data)
rm(Z)
saveRDS(q5, "q5.rds")

# step 1p
datalist <- list()
for (i in 0:9){
  print(i)
  data <- readMM(paste("rqres", i, ".mtx", sep = ""))
  data <- as(data * 1, "dgCMatrix")
  datalist[[i+1]] <- data
}
data <- cbind(datalist[[1]], datalist[[2]],
              datalist[[3]], datalist[[4]],
              datalist[[5]], datalist[[6]],
              datalist[[7]], datalist[[8]],
              datalist[[9]], datalist[[10]])
writeMM(data, "rqres_all.mtx")
rm(data)
rm(datalist)

# step 2
genenames <- readRDS("genenames.rds")
feature_id <- readRDS("feature_idx.rds")
genenames <- genenames[feature_id]
datathres <- readMM("rqres_all.mtx")
datathres <- as(datathres * 1, "dgCMatrix")
q5 <- readRDS("q5.rds")

so <- create_SiFINeT_object(counts = matrix(1:4,2,2), sparse = TRUE, rowfeature = F)
so@data.thres$data1 <- datathres
rm(datathres)
so@q5 <- q5
so <- feature_coexp(so)
so@data.thres <- list()
so <- create_network(so)
so <- filter_lowexp(so)
so <- cal_connectivity(so)
so <- find_unique_feature(so)
so <- assign_shared_feature(so)
so <- enrich_feature_set(so)
saveRDS(so, "so.rds")

temp <- sum(.Internal(gc(FALSE, FALSE, TRUE))[13:14]) - tt
print(temp)
rm(list = ls())
