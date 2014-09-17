data = read.table("dataLogisticNorm3D.txt", header = TRUE)

covariance_matrix = function(d, alpha, beta){
  return(matrix(-beta, d, d) + diag(alpha + beta, d))
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
  mu.hat = colSums(U) / n #note: this is the mle for log points to get mu so we still need
                          #to tranform back but b/c of invariance property of mle's 
                          #this is ok
    
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

tranform_data = function(log_points){
  ### u_dplus_one for each observed point i is 1 - sum(data[i,])
  #then each coordinate in the row becomes log(point/u_dplus_one)
  
}

###################
#Since we took a numerical approach to finding our estimates, let's compare 
#them with results from optim

#This part is currently Out of Order

run.mle = function(data){
  y = data
  n = length(y[,1])
  
  contrained_dlogisticnorm = function(u, mu, alpha, transformed_beta){
    beta = (alpha - transformed_beta)/(length(u) -1)
    return(dlogisticnorm(u, mu, alpha, beta))
  }
  
  
  log.lik = function(par) {
    l = 0
    for (i in 1:n){
      l = l + log(contrained_dlogisticnorm(y[i,], mu=par[1], alpha=par[2], transformed_beta=par[3]))
    }
    return(l)
  }
  
  out = optim(par=c(c(.1,.1,1), 1,1), fn = log.lik, control=list(fnscale=-1),
              method="L-BFGS-B", lower=c(c(1e-6,1e-6,1e-6), 1e-6, 1e-6))
  return(out$par)
  
}
run.mle(data)

