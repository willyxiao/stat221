data_area2_path = 'dat/data1985_area2.csv'
theta0List_path = 'dat/theta0list.Rdata'

data_area2 = read.table(data_area2_path, header=T)

# fix one of the rows in data
data_area2[,2] = data_area2[,2] - 1

# variable name is theta0List
load(theta0List_path)

data = dataToMatrix(data_area2)
G = rep(.5, dim(data)[1])

getLowProbability = function(x, theta){
  theta$low[x]
}
getHighProbability = function(x, theta){
  theta$high[x]
}

#2.2
# G is a N x 1 vector where G[1] is G_Li
# X is a N x J matrix where X[n,j] is the category of plot 1 feature j
# theta is a J x 2 matrix where theta[j][i] j is feature, i -> (1=high, 2=low) 
llV = function(G, theta, X){
  if(all(G < 0 | G > 1)){
    return (-Inf)
  }
  XT = t(X) + 1
  lowProbs = matrix(mapply(getLowProbability,x=XT,theta=theta), nrow=dim(XT)[1], ncol=dim(XT)[2])
  highProbs = matrix(mapply(getHighProbability,x=XT, theta=theta), nrow=dim(XT)[1], ncol=dim(XT)[2])
  retval = sum(log(G*lowProbs + (1 - G)*highProbs))
  return (retval)
}

dataToMatrix = function(data){
  # throw out first column of ids
  return (data[,2:dim(data_area2)[2]])
}

data = dataToMatrix(data_area2)

# G should be an 
ll = function(G, theta, X){
  # penalize for not meeting constraint
  if(all(G < 0 | G > 1)){
    return (-Inf)
  }

  n = length(X) - 1
  total_ll = 0

  for(i in 1:n) { 
    for(j in 1:length(theta)){
      g_Li = G[i]
      g_Hi = 1 - g_Li
      
      theta_Hj = theta[[j]]$high
      theta_Lj = theta[[j]]$low
      
      # penalize for not meeting constraint
      if(!all(theta_Hj > 0) | !all(theta_Lj > 0)){
        return (-Inf)
      }
      
      sum_H = sum(theta_Hj)
      sum_L = sum(theta_Lj)
      
      # TODO how to penalize?
      if(sum_H > 1.1 || sum_H < 0.9){
        return (-Inf)
      } else if(sum_L > 1.1 || sum_L < 0.9){
        return (-Inf)
      }

      k = X[i+1,j + 1] + 1
      
      # print (c(as.integer(i), as.integer(j), g_Hi, g_Li, theta_Hj[k], theta_Lj[k]))
      p = g_Hi * theta_Hj[k] + g_Li * theta_Lj[k]
      total_ll = total_ll + log(p)
    }
  }
  return (total_ll)
}

transform = function(v){
  return(exp(v)/sum(exp(v)))
}

big_negative = -1e200

llNoNegInf = function(G, theta, X){
  #   likelihood = llV(G, theta, X)
  likelihood = ll(G, theta, X)
  retval = ifelse(likelihood == -Inf, big_negative, likelihood)
  return (retval)
}


llOptimG = function(G, theta, X){
  return (llNoNegInf(G,theta,X))
}

llOptimTheta = function(theta_L, theta_H, j, G, theta, X, high_flag){
  if (high_flag){
    theta[[j]]$high = transform(theta_H)
    theta[[j]]$low = theta_L
  }
  else{
    theta[[j]]$high = theta_H
    theta[[j]]$low = transform(theta_L)
  }
  
  theta_Hj = theta[[j]]$high
  theta_Lj = theta[[j]]$low
  
  if(!all(theta_Hj > 0) | !all(theta_Lj > 0)){
    return (big_negative)
  }
  
  sum_H = sum(theta_Hj)
  sum_L = sum(theta_Lj)
  
  # TODO how to penalize?
  if(sum_H > 1.1 || sum_H < 0.9){
    return (big_negative)
  } else if(sum_L > 1.1 || sum_L < 0.9){
    return (big_negative)
  }

  return (llNoNegInf(G, theta, X))
}

#2.3
# OPTIMIZE (gomMLE takes like wayyy too long to run)
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
                  method="L-BFGS-B",
                  X=X, theta=theta, theta_H=theta_Hj, G=G, j=j, high_flag = FALSE, lower = rep(1e-6, length(theta_Lj)),
                  upper = rep(1, length(theta)),
                  control=list(fnscale=-1))
      dummy = transform(res$par)
      theta[[j]]$low = dummy
      print(dummy)
    }
    print("finished lows")
    
    for(j in 1:J){
      theta_j = theta[[j]]
      theta_Lj = theta_j$low
      theta_Hj = theta_j$high  
      res = optim(par=c(theta_H=theta_Hj), 
                  fn=llOptimTheta,
                  method="L-BFGS-B",
                  X=X, theta=theta, theta_L=theta_Lj, j=j, G=G, high_flag = TRUE, lower = rep(1e-6, length(theta_Lj)),
                  upper = rep(1, length(theta)),
                  control=list(fnscale=-1))
      dummy = transform(res$par)
      theta[[j]]$high = dummy
      print(dummy) 
      lik_holder = res$value
    }
    lik1 = lik
    lik = lik_holder    
    print(c(lik=lik,lik1=lik1))
  }

  return (list(G.hat=G, theta.hat=theta, maxlik=lik))
}

# run.gomMLE = function(X=data,G0=FALSE,theta0=theta0List) {
run.gomMLE = function(X=data_area2,G0=FALSE,theta0=theta0List) {
  if(!G0){
    G0 = rep(.5,length(X) - 1)
  }
  return (gomMLE(X,G0,theta0))  
}
#-2120.93     -Inf 
#-2090.045 -2120.930 
#-2086.919 -2090.045 
#-2083.276 -2086.919
#-2078.717 -2083.276 
#-2073.462 -2078.717
#-2067.277 -2073.462 
#-2060.340 -2067.277 
#-2051.935 -2060.340 
#-2043.091 -2051.935 
#-2033.181 -2043.091
#-2021.950 -2033.181 
#-2011.095 -2021.950 
#
#Error in optim(par = c(theta_L = theta_Lj), fn = llOptimTheta, method = "L-BFGS-B",  : 
#                 L-BFGS-B needs finite values of 'fn'
#               In addition: Warning message:
#                 In log(p) : NaNs produced
