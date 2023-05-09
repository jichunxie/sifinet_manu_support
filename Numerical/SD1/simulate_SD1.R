library(GeneScape)

for (seed in 1:100){
  set.seed(seed)
  sim <- GeneScape(nCells = 6000, nGenes = 5000, lib.size.loc = 9.3, lib.size.scale = 0.25, 
  					gene.mean.shape = 0.3, gene.mean.rate = 0.15, 
                    nGroups = 3, groups = rep(1:3, c(100, 4900, 1000)),
                    de.n = 50,
                    de.fc.loc = list(c(rep(0.5, 20), rep(0.3, 30)), 
                                     c(rep(0.3, 30), rep(0.5, 20)),
                                     rep(0.5, 50)), 
                    de.fc.scale = rep(0.2, 3), de.id = list(1:50, 21:70, 71:120),
                    add.sub = FALSE, 
                    add.cor = TRUE, cor.n = 3, cor.cor = 5:7/10, cor.id = list(121:140, 141:160, 161:180), band.width = 10,
                    add.hub = TRUE, hub.n = 20, hub.cor = rep(rep(c(0.3, 0.6), each = 5), 2), hub.id = 201:220, 
                    hub.size = rep(c(10, 20), each = 10),
                    drop = TRUE, dropout.location = -2, dropout.slope = -1)

  Y <- as.matrix(sim[[1]])
  Y <- Y[,colSums(Y > 0) > 0.1 * nrow(Y)]
  idx <- which(rowSums(Y > 0) > 0.1 * ncol(Y))
  Y <- Y[idx,]
  
  saveRDS(Y, paste("data/", seed,"_matrix.rds", sep = ""))
  saveRDS(idx, paste("data/", seed,"_name.rds", sep = ""))
}
