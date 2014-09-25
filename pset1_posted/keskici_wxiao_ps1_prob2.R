library(RUnit)
library(logging)
basicConfig()
addHandler(writeToFile, file="gomMLE.log", level="INFO")

BIG.NEGATIVE = -1e200

data = NA
theta0List = NA

data.to.matrix = function(data){
  noIds = data[,2:dim(data)[2]]
  t(t(noIds)) # to make it a matrix
}

theta.non.zero = function(theta){
  for(i in 1:length(theta)){
    theta[[i]]$high = theta[[i]]$high + .0001
    theta[[i]]$low = theta[[i]]$low + .0001
    theta[[i]]$high = theta[[i]]$high / sum(theta[[i]]$high)
    theta[[i]]$low = theta[[i]]$low / sum(theta[[i]]$low)
  }
  theta
}

load.defaults = function(){
  data_area2_path = 'dat/data1985_area2.csv'
  theta0List_path = 'dat/theta0list.Rdata'
  g0_path = 'dat/g0.Rdata'
  
  data_area2 = read.table(data_area2_path, header=T)
  # fix one of the rows in data
  data_area2[,2] = data_area2[,2] - 1
  
  load(theta0List_path) # theta0List
  load(g0_path) # g0
  
  theta0List <<- theta.non.zero(theta0List)
  data <<- data.to.matrix(data_area2)  
}

load.defaults()

#2.2
ll = function(G, theta, X){
  if(!all(0 < G & G < 1)){
    return (-Inf)
  }
  log.density = c()
  n = length(G)
  if(n != nrow(X)){
    stop("Not correct data size.")
  }
  for(j in 1:ncol(X)){
    theta.Lj = theta[[j]]$low
    theta.Hj = theta[[j]]$high
    log.density = c(log.density, log(G*theta.Lj[X[,j] + 1] + (1 - G) * theta.Hj[X[,j] + 1]))
  }
  sum(log.density)
}

simplex <- function(v){
  u.vector = exp(v) / (1 + sum(exp(v)))
  c(u.vector, 1 - sum(u.vector))
}

logistic <- function(v){
  return(exp(v)/sum(exp(v)))
}

# computes the log likelihood but replaces an NA or -Inf with BIG.NEGATIVE
optim.G <- function(G, Theta, X){
  lik = ll(G, theta, data)
  ifelse(lik == -Inf || is.na(lik), BIG.NEGATIVE, lik)
}

# computes the log-likelihood over a single feature j
ll.theta <- function(G, theta_vec, Xj){
  k.j <- Xj + 1
  sum(log(G*theta_vec[k.j]))
}


llTheta.check <- function(){
  x <- c(0,1,2)
  theta_Lj <- c(.2,.7,.1)
  theta_Hj <- c(.1,.5,.4)
  G <- c(.1,.5,.2)

  #.1 * (.2) + .9 * (.1) = .02 + .09 = .11
  #.5 * (.7) + .5 * (.5) = .35 + .25 = .60
  #.2 * (.1) + .8 * (.4) = .02 + .32 = .34
  checkTrue(round(llTheta(G, theta_Lj, theta_Hj, x), 4) == -3.7969)
}
# G is a N x 1 vector where G[1] is G_Li
# X is a N x J matrix where X[n,j] is the category of plot i feature j
optim.theta <- function(theta_L, theta_H, j, G, theta, X, high_flag){
  #ld_thetas <- theta_L
  if(high_flag){
    #theta_H <- logistic(theta_H)
    return(ll.theta((1-G), logistic(theta_H), X[,j]))
  } else {
    #theta_L <- logistic(theta_L)
    return(ll.theta(G, logistic(theta_L), X[,j]))
  }

#  stopifnot(!all(theta_H > 0 & theta_L > 0))
  #if(!all(theta_H > 0) || !all(theta_L > 0)){
  #  browser()
  #  return (BIG.NEGATIVE)
  #} 
  
  #llTheta(G, theta_L, theta_H, X[,j])
}

#2.3
# OPTIMIZE (gomMLE takes like wayyy too long to run)
gomMLE <- function(X, G0, theta0){
  lik <- -Inf
  lik1 <- 0
  
  G = G0
  theta <- theta0
  
  N <- length(X)
  J <- length(theta)
  
  while(lik != lik1) {
    
    # g_L,n for n = 1,...,N
    loginfo("Running G...")
    print(G)
    res = optim(par=G, fn=optim.G, method="L-BFGS-B", X=X, Theta=theta, control=list(fnscale=-1))
    loginfo("Distance: %f", sum((G - res$par)**2))
    G <- res$par    
    loginfo("Result: %f", res$value)
    print(G)
    
    loginfo("Running lows...")
    # theta_l,j for j = 1,...,J
    for(j in 1:J){
      theta_j <- theta[[j]]
      theta_Lj <- theta_j$low
      theta_Hj <- theta_j$high

      res = optim(par=rep(.001, length(theta_Lj)), 
                  fn=optim.theta,
                  method="L-BFGS-B",
                  X=X, theta_H=theta_Hj, G=G, j=j, high_flag = FALSE, 
                  lower = rep(-10,length(theta_Lj)), upper = rep(10, length(theta_Lj)),
                  control=list(fnscale=-1))
      dummy = logistic(res$par)
      theta[[j]]$low <- dummy

      loginfo("On feature: %2d, loglik: %f", j, ll(G, theta, X))
    }

    loginfo("Running highs")
    for(j in 1:J){
      theta_j <- theta[[j]]
      theta_Lj <- theta_j$low
      theta_Hj <- theta_j$high 
      old_theta_Hj <- theta_Hj
      old_val <- ll(G, theta, X)
#      old_val <- lloptim.theta2check(theta_L=theta_Lj,theta_H=theta_Hj,j=j,G=G,X=X)
      #old_theta <- unlist(theta)
      res <- optim(par=rep(.001, length(theta_Hj)), 
                  fn=optim.theta,
                  method="L-BFGS-B",
                  X=X, theta_L=theta_Lj, G=G, j=j, high_flag = TRUE, 
                  lower = rep(-10,length(theta_Hj)), upper = rep(10, length(theta_Hj)),
                  control=list(fnscale=-1))
      dummy <- logistic(res$par)
      
      theta[[j]]$high <- dummy
      #if(ll(G, theta, X) < old_val){
#        print(old_theta - unlist(theta))
#        print(theta[[j]]$high)
#        print(old_theta_Hj)
        #print(res$par)
      #  browser()
      #}
                  
      loginfo("On feature: %2d, loglik: %f", j, ll(G, theta, X))
      lik_holder <- ll(G, theta, X)
    }
    lik1 <- lik
    lik <- lik_holder    
    loginfo("Old: %f, New: %f", lik1,lik)
    
    MLES <<- list(G.hat=G, theta.hat=theta, maxlik=lik)
  }

  return (MLES)
}

run.gomMLE <- function(X=data,G0=FALSE,theta0=theta0List) {
  if(!G0){
    
    G0 <- rbeta(dim(X)[1],100,100)
    #G0 <- rep(.5,dim(X)[1])
  }
  return (gomMLE(X,G0,theta0))  
}