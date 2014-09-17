data = read.table("dataLogisticNorm3D.txt", header = TRUE)

covariance_matrix = function(d, alpha, beta){
  return(matrix(-beta, d, d) + diag(alpha + beta, d))
}

tranform_data = function(U){
  ### u_dplus_one for each observed point i is 1 - sum(data[i,])
  #then each coordinate in the row becomes log(point/u_dplus_one)
  u_d_plus_one = 1 - rowSums(U) 
  X = U / u_d_plus_one
  return (log(X))
}

#1.3
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
  X = tranform_data(U)
  d = length(X[1,])
  n = length(X[,1])
  mu.hat = colSums(X) / n 
    
  temp = 0
  for (i in 1:d){
    temp = temp + var(X[,i])
  }
  alpha.hat = temp / d
  
  temp = 0
  for (i in 1:d){
    j = i +1
    while(j <= d){
      temp = temp + cov(X[,i],X[,j])
      j = j + 1
    }
  }
  beta.hat = -temp / choose(d,2)

  MLEs = list("mu.hat" = mu.hat, "alpha.hat" = alpha.hat, "beta.hat" = beta.hat)
  return(MLEs)  
}



#1.4
logisticnorm.mle(data)

###################
#Since we took a numerical approach to finding our estimates in 1.2, let's 
#simulate data and compare MLE's with Parameters

library(mvtnorm)
run.mle_check = function(n, mu, sigma){
  sim_data = rmvnorm(n, mean = mu, sigma = sigma)
  sim_data = exp(sim_data)
  sim_data = sim_data/(1 + rowSums(sim_data))
  
  result = logisticnorm.mle(sim_data)
  return(result)
  
}

run.mle_check(250, c(.1,.2,.3), covariance_matrix(3,.8,.3))
