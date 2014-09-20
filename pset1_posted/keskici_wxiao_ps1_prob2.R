data_area2_path = 'dat/data1985_area2.csv'
theta0List_path = 'dat/theta0list.Rdata'

data_area2 = read.table(data_area2_path, header=T)
data_area2[,2] = data_area2[,2] - 1
load(theta0List_path) #Note: variable name is theta0List

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

llNoNegInf = function(G, theta, X){
  likelihood = ll(G, theta, X)
  return (ifelse(likelihood == -Inf, -.Machine$double.xmax, likelihood))
}


llOptimG = function(g, theta, X){
  return (llNoNegInf(g,theta,X))
}

llOptimTheta = function(theta_L, theta_H, j, g, theta, X){
  theta[[j]]$high = theta_H
  theta[[j]]$low = theta_L
  return (llNoNegInf(g, theta, X))
}

#2.3
gomMLE = function(X, G0, theta0){
  lik = -Inf
  lik1 = 0

  G = G0
  theta = theta0
  
  N = length(X) - 1
  J = length(theta)
  
  while(lik != lik1) {
    
    # g_L,n for n = 1,...,N
    res = optim(par=c(g=G), fn=llOptimG, X=X, theta=theta, control=list(fnscale=-1))
    G = res$par
    
    # theta_l,j for j = 1,...,J
    for(j in 1:J){
      theta_j = theta[[j]]
      theta_Lj = theta_j$low
      theta_Hj = theta_j$high
      res = optim(par=c(theta_L=theta_Lj), 
                  fn=llOptimTheta,
                  X=X, theta=theta, theta_H=theta_Hj, g=G, j=j,
                  control=list(fnscale=-1))
      theta[[j]]$low = res$par
    }
    
    for(j in 1:J){
      theta_Hj = theta[[j]]$high
      res = optim(par=c(theta_H=theta_Hj), 
                  fn=llOptimTheta,
                  X=X, theta=theta, theta_H=theta_Hj, j=j, g=G,
                  control=list(fnscale=-1))
      theta[[j]]$high = res$par
    
      lik1 = lik
      lik = res$value
    }
    

  }

  return (list(G.hat=G, theta.hat=theta, maxlik=lik))
}

run.gomMLE = function(X=data_area2,G0=FALSE,theta0=theta0List) {
  if(!G0){
    G0 = rep(.5,length(X) - 1)
  }
  return (gomMLE(X,G0,theta0))  
}