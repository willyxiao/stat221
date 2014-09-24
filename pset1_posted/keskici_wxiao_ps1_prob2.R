library(RUnit)
library(logging)
basicConfig()
addHandler(writeToFile, file="gomMLE.log", level="INFO")

sourceMe <- function() {
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

makeNonZero = function(theta){
  for(i in 1:length(theta)){
    theta[[i]]$high = theta[[i]]$high + .0001
    theta[[i]]$low = theta[[i]]$low + .0001
    theta[[i]]$high = theta[[i]]$high / sum(theta[[i]]$high)
    theta[[i]]$low = theta[[i]]$low / sum(theta[[i]]$low)
  }
  theta
}
theta0List = makeNonZero(theta0List)

getLowProbability = function(x, theta){
  theta$low[x]
}
getHighProbability = function(x, theta){
  theta$high[x]
}

#2.2
# G is a N x 1 vector where G[1] is G_Li
# X is a N x J matrix where X[n,j] is the category of plot 1 feature j
llV2 = function(G, theta, X){
  if(!all(0 < G & G < 1)){
    return (-Inf)
  }
  XT = t(X) + 1
  lowProbs = matrix(mapply(getLowProbability,x=XT,theta=theta), nrow=dim(XT)[1], ncol=dim(XT)[2])
  highProbs = matrix(mapply(getHighProbability,x=XT, theta=theta), nrow=dim(XT)[1], ncol=dim(XT)[2])
  liks = G*lowProbs + (1-G)*highProbs
  retval = sum(log(liks))
  return (retval)
}

llV = function(G, theta, X){
  log.density = c()
  n = length(G)
  if(n != nrow(X)){
    stop("Not correct data size.")
  }
  for(j in 1:ncol(X)){
    theta.Lj = theta[[j]]$low
    theta.Hj = theta[[j]]$high
    log.density = c(log.density, log(G * theta.Lj[X[,j] + 1] + (1 - G) * theta.Hj[X[,j] + 1]))
  }
  sum(log.density)
}
dataToMatrix = function(data){
  noIds = data[,2:dim(data_area2)[2]]
  t(t(noIds)) # to make it a matrix
}

data = dataToMatrix(data_area2)

logistic <- function(v){
  return(exp(v)/sum(exp(v)))
}

big_negative = -1e200

llNoNegInf = function(G, theta, X){
  lik = llV(G, theta, data)
  retval = ifelse(lik == -Inf || is.na(lik), big_negative, lik)
  return (retval)
}


llOptimG <- function(G, theta, X){
  return (llNoNegInf(G,theta,X))
}

optimG <- function(g, X, i ){
    
}

# G should be an N x 1 vector
llTheta <- function(G, theta_Lj, theta_Hj, Xj){
  k_j <- Xj + 1
  sum(log(G*theta_Lj[k_j] + (1-G)*theta_Hj[k_j]))
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
optimTheta <- function(theta_L, theta_H, j, G, theta, X, high_flag){
  old_thetas <- theta_L
  
  if(high_flag){
    theta_H <- logistic(theta_H)
  } else {
    theta_L <- logistic(theta_L)
  }

#  stopifnot(!all(theta_H > 0 & theta_L > 0))
  if(!all(theta_H > 0) || !all(theta_L > 0)){
    browser()
    return (big_negative)
  } 
  
  llTheta(G, theta_L, theta_H, X[,j])
}

#2.3
# OPTIMIZE (gomMLE takes like wayyy too long to run)
gomMLE <- function(X, G0, theta0){
  lik <- -Inf
  lik1 <- 0
  
  G <- G0
  theta <- theta0
  
  N <- length(X)
  J <- length(theta)
  
  while(lik != lik1) {
    
    # g_L,n for n = 1,...,N
    loginfo("Running G...")
    res <- optim(par=c(G=G), fn=llOptimG, method="L-BFGS-B", X=X, theta=theta, control=list(fnscale=-1))
    loginfo("Distance: %f", sum((G - res$par)**2))
    G <- res$par    
    loginfo("Result: %f", res$value)
    
    loginfo("Running lows...")
    # theta_l,j for j = 1,...,J
    for(j in 1:J){
      theta_j <- theta[[j]]
      theta_Lj <- theta_j$low
      theta_Hj <- theta_j$high

      res <- optim(par=c(theta_L= rep(.001, length(theta_Lj))), 
                  fn=optimTheta,
                  method="L-BFGS-B",
                  X=X, theta_H=theta_Hj, G=G, j=j, high_flag = FALSE, 
                  lower = rep(-20,length(theta_Lj)), upper = rep(20, length(theta_Lj)),
                  control=list(fnscale=-1))
      dummy <- logistic(res$par)
      theta[[j]]$low <- dummy

      loginfo("On feature: %2d, loglik: %f", j, llV(G, theta, X))
    }

    loginfo("Running highs")
    for(j in 1:J){
      theta_j <- theta[[j]]
      theta_Lj <- theta_j$low
      theta_Hj <- theta_j$high 
      old_theta_Hj <- theta_Hj
      old_val <- llV(G, theta, X)
#      old_val <- llOptimTheta2check(theta_L=theta_Lj,theta_H=theta_Hj,j=j,G=G,X=X)
      old_theta <- unlist(theta)
      res <- optim(par=c(theta_H=rep(.001, length(theta_Hj))), 
                  fn=optimTheta,
                  method="L-BFGS-B",
                  X=X, theta_L=theta_Lj, G=G, j=j, high_flag = TRUE, 
                  lower = rep(-20,length(theta_Hj)), upper = rep(20, length(theta_Hj)),
                  control=list(fnscale=-1))
      dummy <- logistic(res$par)
      
      theta[[j]]$high <- dummy
      if(llV(G, theta, X) < old_val){
#        print(old_theta - unlist(theta))
#        print(theta[[j]]$high)
#        print(old_theta_Hj)
        #print(res$par)
        browser()
      }
                  
      loginfo("On feature: %2d, loglik: %f", j, llV(G, theta, X))
      lik_holder <- llV(G, theta, X)
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
    G0 <- rep(.5,dim(X)[1])
  }
  return (gomMLE(X,G0,theta0))  
}