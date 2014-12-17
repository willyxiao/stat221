args <- as.numeric(commandArgs(trailingOnly = TRUE))
taskID <- args[1]

simulate.data = function(nsamples){
  data = rep(0, nsamples)
  class = rep(0, nsamples)
  Us = rep(0, nsamples)
  
 for(i in 1:nsamples){
    U = runif(1, 0, 10)
    Us[i] = U
    V = rnorm(1, 0, 9)
    
    if(runif(1, 0, 1) > .5){
      class[i] = 1
      data[i] = 5*U + V
    } else {
      class[i] = 2
      data[i] = 15 + 10*U - U**2 + V
    }
  }
  list(data=data, class=class, U=Us)
}


nRuns <- 1
nSamples <- 10000
beta1.mat <- matrix(nrow = nRuns, ncol = 3)
beta2.mat <- matrix(nrow = nRuns, ncol = 3)

library(segmented, lib.loc = "~/R/library")
library(mixtools, lib.loc = "~/R/library")
for(i in 1:nRuns){
  sim <- simulate.data(nSamples)
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

save(beta1.mat, file = sprintf("beta1_%i.rda", taskID))
save(beta2.mat, file = sprintf("beta2_%i.rda", taskID))
