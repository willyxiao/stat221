library(logging)
#basicConfig()
#addHandler(writeToFile, file="gomMLE.log", level="INFO")

sourceMe = function() {
  source("./keskici_wxiao_ps1_prob2.R")
}


data_area2_path = 'dat/data1985_area2.csv'
theta0List_path = 'dat/theta0list.Rdata'

data_area2 = read.table(data_area2_path, header=T)

# fix one of the rows in data
data_area2[,2] = data_area2[,2] - 1

# variable name is theta0List
load(theta0List_path)
load("./dat/g0.Rdata")

dataToMatrix = function(data){
  noIds = data[,2:dim(data_area2)[2]]
  t(t(noIds)) # to make it a matrix
}

makeNonZero = function(theta){
  for(i in 1:length(theta)){
    theta[[i]]$high = theta[[i]]$high + .01
    theta[[i]]$low = theta[[i]]$low + .01
  }
  theta
}
theta0List = makeNonZero(theta0List)

data = dataToMatrix(data_area2)

getLowProbability = function(x, theta){
  theta$low[x]
}
getHighProbability = function(x, theta){
  theta$high[x]
}

#2.2
# G is a N x 1 vector where G[1] is G_Li
# X is a N x J matrix where X[n,j] is the category of plot 1 feature j
llV = function(G, theta, X){
  if(!all(0 <= G & G <= 1)){
    return (-Inf)
  }
  XT = t(X) + 1
  lowProbs = matrix(mapply(getLowProbability,x=XT,theta=theta), nrow=dim(XT)[1], ncol=dim(XT)[2])
  highProbs = matrix(mapply(getHighProbability,x=XT, theta=theta), nrow=dim(XT)[1], ncol=dim(XT)[2])
  liks = G*lowProbs + (1-G)*highProbs
  retval = sum(log(liks))
  if(retval == -Inf){
    browser()
  }
  return (retval)
}

# G should be an 
ll = function(G, theta, X){
  # penalize for not meeting constraint
  if(all(G < 0 | G > 1)){
    print("ll: G out of bounds")
    return (-Inf)
  }

  n = dim(X)[1]
  total_ll = 0

  # m = matrix(0,n,length(theta))
  
  for(i in 1:n) { 
    # res = c()
    for(j in 1:length(theta)){
      g_Li = G[i]
      g_Hi = 1 - g_Li
      
      theta_Hj = theta[[j]]$high
      theta_Lj = theta[[j]]$low
      
      # penalize for not meeting constraint
      if(!all(theta_Hj > 0) | !all(theta_Lj > 0)){
        browser()
        print("ll broken range")
        return (-Inf)
      }
      
      sum_H = sum(theta_Hj)
      sum_L = sum(theta_Lj)
      
      # TODO how to penalize?
      if(sum_H > 1.1 || sum_H < 0.9){
        print("ll broken range2")
        return (-Inf)
      } else if(sum_L > 1.1 || sum_L < 0.9){
        print("ll broken range3")
        return (-Inf)
      }

      k = X[i,j + 1] + 1
      
      # print (c(as.integer(i), as.integer(j), g_Hi, g_Li, theta_Hj[k], theta_Lj[k]))
      p = g_Hi * theta_Hj[k] + g_Li * theta_Lj[k]
      # res[length(res) + 1] = p
      total_ll = total_ll + log(p)
    }
    # m[i,] = res
  }
  # browser()
  return (total_ll)
}

transform = function(v){
  return(exp(v)/sum(exp(v)))
}

big_negative = -1e200

llNoNegInf = function(G, theta, X){
  lik = llV(G, theta, data)
  retval = ifelse(lik == -Inf || is.na(lik), big_negative, lik)
  return (retval)
}


llOptimG = function(G, theta, X){
  return (llNoNegInf(G,theta,X))
}

# G is a N x 1 vector where G[1] is G_Li
# X is a N x J matrix where X[n,j] is the category of plot i feature j
llOptimTheta2 = function(theta_L, theta_H, j, G, theta, X, high_flag){
  k_j = X[,j] + 1
  
  if(high_flag){
    theta_H = transform(theta_H)
  } else {
    theta_L = transform(theta_L)
  }
  
  retval = sum(log(G*theta_L[k_j] + (1-G)*theta_H[k_j]))
  retval
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

llOptimG2 = function(g, theta, X, i){
  if(!(0 <= g & g <= 1)){
    return (big_negative)
  }
  lowProbs = mapply(getLowProbability, x=X[i,] + 1,theta=theta)
  highProbs = mapply(getHighProbability, x=X[i,] + 1,theta=theta)
  browser()
  sum(log(g*lowProbs + (1-g)*highProbs))
}

#2.3
# OPTIMIZE (gomMLE takes like wayyy too long to run)
gomMLE = function(X, G0, theta0){
  lik = -Inf
  lik1 = 0

  G = G0
  theta = theta0
  
  #N = length(X) - 1
  N = length(G)
  J = length(theta)
  
  while(lik != lik1) {
    
    # g_L,n for n = 1,...,N
    loginfo("Running G...")
    G_old = G
    for(i in 1:N){
      res = optim(par=c(g=G[i]), 
                  fn=llOptimG2, 
                  method="L-BFGS-B", 
                  X=X, theta=theta, i=i,
                  lower=0, 
                  upper=1, 
                  control=list(fnscale=-1))
      G[i] = res$par
#      loginfo("G_L%d is %f", i, res$par)
    }
    
#    res = optim(par=c(g=G), fn=llOptimG, method="L-BFGS-B", X=X, theta=theta, control=list(fnscale=-1))
    loginfo("Distance: %f", mean(abs(G - G_old)))
#    G = res$par
#    G <<- res$par
    loginfo("Result: %f", llV(G,theta,X))

    loginfo("Running lows...")
    # theta_l,j for j = 1,...,J
    for(j in 1:J){
      theta_j = theta[[j]]
      theta_Lj = theta_j$low
      theta_Hj = theta_j$high
      res = optim(par=c(theta_L=theta_Lj), 
                  fn=llOptimTheta2,
                  method="L-BFGS-B",
                  X=X, theta=theta, theta_H=theta_Hj, G=G, j=j, high_flag = FALSE,
                  control=list(fnscale=-1))
      dummy = transform(res$par)
      theta[[j]]$low = dummy
#      loginfo("On feature: %2d, loglik: %f", j, llV(G, theta, X))
    }

    loginfo("Running highs")
    for(j in 1:J){
      theta_j = theta[[j]]
      theta_Lj = theta_j$low
      theta_Hj = theta_j$high  
      res = optim(par=c(theta_H=theta_Hj), 
                  fn=llOptimTheta2,
                  method="L-BFGS-B",
                  X=X, theta=theta, theta_L=theta_Lj, j=j, G=G, high_flag = TRUE,
                  control=list(fnscale=-1))
      dummy = transform(res$par)
      theta[[j]]$high = dummy
#      loginfo("On feature: %2d, loglik: %f", j, llV(G, theta, X))
      lik_holder = llV(G, theta, X)
    }
    lik1 = lik
    lik = lik_holder    
    loginfo("Old: %f, New: %f", lik1,lik)
    
    MLES <<- list(G.hat=G, theta.hat=theta, maxlik=lik)
    browser()
  }

  return (MLES)
}

run.gomMLE = function(X=data,G0=G,theta0=theta0List) {
# run.gomMLE = function(X=data_area2,G0=FALSE,theta0=theta0List) {
  if(!G0){
    G0 = rep(.5,dim(X)[1])
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
