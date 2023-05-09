library(GeneScape)

set.seed(161334)
sim <- GeneScape(nCells = 6000, nGenes = 5000, lib.size.loc = 9.3, lib.size.scale = 0.25, 
  				  gene.mean.shape = 0.3, gene.mean.rate = 0.15,
  				  nGroups = 9, groups = rep(1:9, c(1200, 400, 400, 400, 1200, 400, 400, 400, 1200)),
  				  de.n = 60,
                  de.fc.loc = list(c(rep(0.8, 40), rep(0.6, 20)),
                                   c(rep(0.6, 20), rep(0.8, 20), rep(0.4, 20)),
                                   c(rep(0.5, 20), rep(0.8, 20), rep(0.5, 20)),
                                   c(rep(0.4, 20), rep(0.8, 20), rep(0.6, 20)),
                                   c(rep(0.6, 20), rep(0.8, 20), rep(0.6, 20)),
                                   c(rep(0.6, 20), rep(0.8, 20), rep(0.4, 20)),
                                   c(rep(0.5, 20), rep(0.8, 20), rep(0.5, 20)),
                                   c(rep(0.4, 20), rep(0.8, 20), rep(0.6, 20)),
                                   c(rep(0.6, 20), rep(0.8, 40))), 
                  de.fc.scale = rep(0.2, 9), 
                  de.id = list(1:60, 
                               41:100, 41:100, 41:100, 
                               81:140,
                               121:180, 121:180, 121:180,
                               161:220),
                  add.sub = FALSE, 
                  add.cor = TRUE, cor.n = 3, cor.cor = 5:7/10, cor.id = list(401:420, 421:440, 441:460), band.width = 10,
                  add.hub = TRUE, hub.n = 20, hub.cor = rep(rep(c(0.3, 0.6), each = 5), 2), 
                  hub.id = 481:500, hub.size = rep(c(10, 20), each = 10), 
                  drop = TRUE, dropout.location = -2, dropout.slope = -1)

Y <- as.matrix(sim[[1]])
Y <- Y[,colSums(Y > 0) > 0.1 * nrow(Y)]
idx <- which(rowSums(Y > 0) > 0.1 * ncol(Y))
Y <- Y[idx,]
  
saveRDS(Y, "matrix.rds")
saveRDS(idx, "name.rds")
