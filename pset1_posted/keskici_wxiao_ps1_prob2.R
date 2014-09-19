data_area2_path = 'dat/data1985_area2.csv'
theta0List_path = 'dat/theta0list.Rdata'

#2.2
# G should be a vector of G_low
# theta is a matrix and X is a matrix
ll = function(G, theta, X){
  n = length(X) - 1
  total_ll = 0

  for(i in 1:n) { 
    for(j in 1:length(theta)){
      g_Li = G[i]
      g_Hi = 1 - g_Li
      
      theta_Hj = theta[[j]]$high
      theta_Lj = theta[[j]]$low
      
      k = X[i+1,j + 1] + 1
      
      # print (c(as.integer(i), as.integer(j), g_Hi, g_Li, theta_Hj[k], theta_Lj[k]))
      p = g_Hi * theta_Hj[k] + g_Li * theta_Lj[k]
      total_ll = total_ll + log(p)
    }
  }
  return (total_ll)
}

data_area2 = read.table(data_area2_path, header=T)
data_area2[,2] = data_area2[,2] - 1
load(theta0List_path) #Note: variable name is theta0List

g = matrix(.5, 2, length(data_area2) - 1)

llOptim = function(g, theta_low, theta_high, X){
  big_negative = -Inf  

  theta = combine(theta_low, theta_high)
  likelihood = ll(g,theta, X)
  
  return (ifelse(likelihood == -Inf, big_negative, likelihood))
}

#2.3
gomMLE = function(X, G0, theta0){
  
  return (optim(par=list(G=G0, theta=theta0), fn=llOptim, X=X, control=list(fnscale=-1))$par)
}