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
  
  if(sum(u) >= 1 | length(u) != length(mu)){
    return (0) #do we want an error instead?
  }
  
  u_d_plus_one = 1 - sum(u)
  
  a = abs(det(2*pi*cov_mtrx))^-.5
  b = 1/(prod(u)*(u_d_plus_one))
  c = t(log(u/u_d_plus_one) - mu) %*% (1/cov_mtrx) %*% (log(u/u_d_plus_one) - mu)
  return (a*b*exp(-.5*c))
}

