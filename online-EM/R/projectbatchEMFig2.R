nRuns <- 500
nSamples <- 100
sim <- simulate.data(nSamples)
beta1.mat <- matrix(nrow = nRuns, ncol = 3)
beta2.mat <- matrix(nrow = nRuns, ncol = 3)


library(mixtools)
for(i in 1:nRuns){
  r <- sim$data
  U <- sim$U
  z <- matrix(nrow = nSamples, ncol = 3)
  z[,1] <- rep(1, nSamples)
  z[,2] <- U
  z[,3] <- U^2/10 
  n <- length(r)
  beta <- matrix(nrow = 3, ncol = 2, data = c(0,5, 0, 15, 10, -10))
  result <- regmixEM(r, z[,2:3], lambda = c(.5,.5), beta = beta, addintercept = T)
  beta1.mat[i,] <- result$beta[,1]
  beta2.mat[i,] <- result$beta[,2]
}

save(beta1.mat, file = "beta1EMFig2.rda")
save(beta2.mat, file = "beta2EMFig2.rda")
