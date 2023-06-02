library(MASS)
set.seed(1)
sim1 <- mvrnorm(500, c(-1,1,1), matrix(c(1, 0, 0, 0, 1, 0, 0, 0, 1), 3, 3))
sim2 <- mvrnorm(500, c(-1,2,3), matrix(c(1, 0.5, 0, 0.5, 1, 0.5, 0, 0.5, 1), 3, 3))
sim3 <- mvrnorm(500, c(-2,3,1), matrix(c(1, 0.5, 0, 0.5, 1, 0, 0, 0, 1), 3, 3))
sim4 <- mvrnorm(500, c(-4,2,1), matrix(c(1, 0, 0, 0, 1, 0.5, 0, 0.5, 1), 3, 3))
sim <- rbind(sim1, sim2, sim3, sim4)
res <- prcomp(sim)
textx <- c(mean(res$x[1:500,1]), mean(res$x[501:1000,1]), 
           mean(res$x[1001:1500,1]), mean(res$x[1501:2000,1]))
texty <- c(mean(res$x[1:500,2]), mean(res$x[501:1000,2]), 
           mean(res$x[1001:1500,2]), mean(res$x[1501:2000,2]))
textc <- c("CS1", "CS3", "CS2", "CS4")
jpeg("Fig1_sim.jpeg", width = 4, height = 3, units = "in", res = 300)
par(mar = c(0.1, 0.1, 0.1, 0.1))
plot(res$x[,1], res$x[,2], 
     col = rep(c("green", "light blue", "hotpink", "yellow"), each = 500), 
     cex = 0.5, pch = 16, xlab = "", ylab = "", xaxt = "n", yaxt = "n")
text(textx, texty, textc, cex=1.6, col="black")
dev.off()
