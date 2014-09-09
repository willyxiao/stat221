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
  2*pi*covariance_matrix
}