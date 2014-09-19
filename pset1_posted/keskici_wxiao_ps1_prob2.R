data_area2_path = 'dat/data1985_area2.csv'
theta0List_path = 'dat/theta0list.Rdata'

#2.2
ll = function(G, theta, X){
  n = length(X) - 1
  # FIXME should g_h be G$high and g_l be G$low?
  g_h = G[1,]
  g_l = G[2,]
#  g_h = rep(.5, n)
#  g_l = rep(.5, n)
#  theta_h = theta[1]
#  theta_l = theta[2]

  total_ll = 0
  for(i in 1:n) { 
    for(j in 1:length(theta)){
      g_Hi = g_h[i]
      g_Li = g_l[i]
      
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

llOptim = function(par, X){
  G = par$G
  theta = par$theta
  return (ll(G, theta, X))
}

#2.3
gomMLE = function(X, G0, theta0){
  return (optim(par=list(G=G0, theta=theta0), fn=llOptim, X=X, control=list(fnscale=-1))$par)
}