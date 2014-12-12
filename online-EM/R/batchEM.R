rm(list = ls())

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

comp.ll <- function(par, z, r, w, ind){
  beta1 <- as.matrix(par[1:3])
  beta2 <- as.matrix(par[4:6])
  sigma2 <- par[7:8]
  stopifnot(sigma2[1] > 0, sigma2[2] > 0)

  ll <- sum((log(w) - (1/2))*(log(sigma2[1]) + 1/(sigma2[1])*(r - t(beta1)%*%t(z))^2)*ind +
        (log(1 - w) - (1/2)*(log(sigma2[2]) + 1/(sigma2[2])*(r - t(beta2)%*%t(z))^2))*(1-ind))
  if(is.na(ll)){
    browser()
  }
  return(ll)
}


batch.EM <- function(r, z, w, sigma2, beta1, beta2){
  nSteps <- length(r)
  beta1.mat <- matrix(nrow = 3, ncol = nSteps)
  beta2.mat <- matrix(nrow = 3, ncol = nSteps)
  ind <- rep(NA, length(r))
  sep <- quantile(r, 1 - w)
  ind[r <= sep] <- 0
  ind[r > sep] <- 1
  
  for(n in 1:nSteps){
  #E-step
  w <- sum((w/sigma2[1])*exp((-w/2*sigma2[1])*(r - (beta1%*%t(z)))^2))/
       sum(((w/sigma2[1])*exp((-w/2*sigma2[1])*(r - (beta1%*%t(z)))^2) + ((1-w)/sigma2[1])*exp((-(1-w)/2*sigma2[2])*(r - (beta2%*%t(z)))^2)))   
  if(w == 1){
    w <- 0.99
  } else if(w == 0){
    w <- 0.01
  }
  #M-step
  opt <- optim(par = c(beta1, beta2, sigma2), fn = comp.ll, control = list(fnscale = -1), method = "L-BFGS-B",
                      z = z,  r = r, w = w, ind = ind, lower = c(rep(-1e8,6), 1e-10, 1e-10),upper = c(rep(1e8,6), 1e8, 1e8))
  beta1 <- opt$par[1:3]
  beta2 <- opt$par[4:6]
  sigma2 <- opt$par[7:8]
  sep <- quantile(r, 1 - w)
  ind[r <= sep] <- 0
  ind[r > sep] <- 1
  #store results
  beta1.mat[,n] <- beta1
  beta2.mat[,n] <- beta2
  
  print(n)
  }
  return(list(beta1.mat = beta1.mat, beta2.mat = beta2.mat))

}


#Run EM for fig 2
system.time()
nRuns <- 1
beta2 <- matrix(nrow = 3, ncol = nRuns)
beta1 <- matrix(nrow = 3, ncol = nRuns)
for(i in 1:nRuns){
  nSamples <- 10000
  sim <- simulate.data(nSamples)
  r <- sim$data
  U <- sim$U
  z <- matrix(nrow = nSamples, ncol = 3)
  z[,1] <- rep(1, nSamples)
  z[,2] <- U
  z[,3] <- U^2/10
  w.0 <- .5
  sigma2.0 <- rep(1,2)
  beta1.0 <- rep(1,3)
  beta2.0 <- rep(1,3)
  result <- batch.EM(r, z, w.0, sigma2.0, beta1.0, beta2.0)
  beta1[,i] <- result$beta1.mat[,100]
  beta2[,i] <- result$beta2.mat[,100]
  print(c("looped",i))
}
system.time()

setwd("C:/Users/Hannah/Dropbox/Stat 221/Project")
save(beta1, file = "beta1EM.rda")
save(beta2, file = "beta2EM.rda")


