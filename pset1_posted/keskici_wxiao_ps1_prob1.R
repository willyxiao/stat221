data = read.table("dataLogisticNorm3D.txt", header = TRUE)

covariance_matrix = function(d, alpha, beta){
  entries = rep(-beta, d*d)
  m = matrix(entries,nrow = d, ncol = d)
  for(i in 1:d) {
    m[i,i] = alpha
  }
  return(m)
}
  
dlogisticnorm = function(u, mu, alpha, beta){
  cov_mtrx = covariance_matrix(length(u),alpha, beta)
  
  if(sum(u) >= 1){
    return (0) #do we want an error instead?
  }
  
  u_d_plus_one = 1 - sum(u)
  
  a = abs(det(2*pi*cov_mtrx))^-.5
  b = 1/(prod(u)*(u_d_plus_one))
  c = as.matrix(log(u/u_d_plus_one) - mu) %*% (1/cov_mtrx) %*% t(log(u/u_d_plus_one) - mu)
  return (a*b*exp(-.5*c))
}

logisticnorm.mle = function(U){
  d = length(U)
  n = length(U[,1])
  mu.hat = colSums(U) / n
    
  temp = 0
  for (i in 1:d){
    temp = temp + var(U[,i])
  }
  alpha.hat = temp / d
  
  temp = 0
  for (i in 1:d){
    j = i +1
    while(j <= d){
      temp = temp + cov(U[,i],U[,j])
      j = j + 1
    }
  }
  beta.hat = -temp / choose(d,2)

  MLEs = list("mu.hat" = mu.hat, "alpha.hat" = alpha.hat, "beta.hat" = beta.hat)
  return(MLEs)  
}

###################
#Since we took a numerical approach to finding our estimates, let's compare 
#them with results from optim

#This part is currently Out of Order

run.mle = function(data){
  y = data
  n = length(y[,1])
  
  log.lik = function(par) {
    l = 0
    for (i in 1:n){
      l = l + log(dlogisticnorm(y[i,], mu=par[1], alpha=par[2], beta=par[3]))
    }
    return(l)
  }
  
  out = optim(par=c(c(1,1,1), 1,1), fn = log.lik, control=list(fnscale=-1), 
              method="L-BFGS-B")
  return(out$par)
  
}


