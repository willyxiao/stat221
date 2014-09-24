rm(list=ls())  # remove current vars in the environment

library(RUnit)
library(logging)
basicConfig()
addHandler(writeToFile, file="gomMLE.log", level="INFO")

sourceMe <- function() {
  source("./keskici_wxiao_ps1_prob2c.R")
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

dataToMatrix = function(data){
  noIds = data[,2:dim(data_area2)[2]]
  t(t(noIds)) # to make it a matrix
}

data = dataToMatrix(data_area2)

get.feature.lengths = function(theta0){
  #For a theta0 list get lengths of each feature
  #theta0 list's categories assumed to be $high and $low
  lengths = c()
  for (i in 1: length(theta0)){
    lengths = c(lengths, length(theta0[[i]]$high))
  }
  return(lengths)
}

feature.len = get.feature.lengths(theta0List) 

num.features = length(feature.len)  # no. of features=49 for 2.c

flatten.thetaList = function(thetaList){
  high = c()
  low = c()
  for (i in 1: length(thetaList)){
    high = c(high, thetaList[[i]]$high)
    low = c(low, thetaList[[i]]$low)
  }
  return(list(high=high, low=low))
}

logistic <- function(v) {
  exp(v) / sum(exp(v))
}

ll <- function(G, theta.low, theta.high, X) {
  # Computes the log-likelihood.
  log.density = c()
  n = length(G)
  if(n != nrow(X)) {
    stop("Not correct data size.")
  }
  for(p in 1:ncol(X)) {
    # for each feature
    # Get the probabilities for the feature.
    theta.Lp = get.aux.theta(theta.low, theta.high, 0, p) # low-risk probs for feature p
    theta.Hp = get.aux.theta(theta.low, theta.high, 1, p) # high-risk probs for feature p
    log.density = c(log.density, log(G * theta.Lp[X[,p]] + (1-G) * theta.Hp[X[,p]]))
  }
  return(sum(log.density))
}

get.feature.thetaIndex <- function(p) {
  # For a specific feature get the vector if indices in the theta vector 
  #   corresponding to that feature.
  #
  # e.g. get.feature.thetaIndex(2) = (4, 5)
  if(! p %in% 1:num.features) {
    stop("Wrong feature id.")
  }
  # Given feature, returns the indices of the feature in a theta vector (either H or L)
  theta.index = c(1 + sum(head(feature.len, p-1)), sum(head(feature.len, p-1)) + feature.len[p])
  theta.p = seq(theta.index[1], theta.index[2])
  return(theta.p)
}

get.aux.theta <- function(theta.L, theta.H, LoHi, p) {
  # Args:
  #   theta.L = vector of Low-risk probs
  #   theta.H = vector of High-risk probs
  #   LoHi = {0 , 1} where 0=Low
  #   p = feature
  # 
  #  Returns the part of the vector of theta that corresponds to that feature.
  #
  theta.p = get.feature.thetaIndex(p)
  if(LoHi == 0) {
    return(theta.L[theta.p])
  } else {
    return(theta.H[theta.p])
  }
}

update.theta <- function(theta, p, theta_p.new) {
  ## Updates a theta vector for a specific feature based on the new 
  # theta vector of that particular feature.
  # Returns a new theta vector.
  # 
  # e.g. update.theta(c(1, 2, 3, 4, 5), p=1, c(0, 0, 0)) = (0, 0, 0, 0.05, 0.95)
  #
  if(length(theta) != sum(feature.len))
    stop("Wrong theta vector.")
  theta.index = get.feature.thetaIndex(p)
  theta_p.new <- theta_p.new / sum(theta_p.new) # make sure to normalize.
  theta[theta.index] <- theta_p.new
  return(theta)
}

gomMLE <- function(X, G0, theta0){
  thetas = flatten.thetaList(theta0)
  theta.low = thetas$low
  theta.high = thetas$high
  num.plots = nrow(X)
  print(ll(G0,theta.low, theta.high, X))
  
  
  optim.fn <- function(par, optimize.for, g.old, thetaL.old, thetaH.old) {
    if(optimize.for=="g") {
      ll(par, thetaL.old, thetaH.old, X)
    } else if(optimize.for=="theta.low") {
      par = c(0, par)
      thetaLp.new = logistic(par)
      thetaL.new = update.theta(thetaL.old, p, thetaLp.new)
      ll(g.old, thetaL.new, thetaH.old, X)
    } else if(optimize.for=="theta.high") {
      par = c(0, par)
      thetaHp.new = logistic(par)
      thetaH.new = update.theta(thetaH.old, p, thetaHp.new)
      ll(g.old, thetaL.old, thetaH.new, X)
    } else {
      stop("Not sure what to optimize for ") 
    }
  }
  
  zeta.lb = c(-10)  # lower-bound for the aux variables (for transformation)
  zeta.ub = c(10)  # upper-bound
  
  
  loginfo("Running G...")
  #res <- optim(par=c(G=G), fn=llOptimG, method="L-BFGS-B", X=X, theta=theta, control=list(fnscale=-1))
  
  # 1. Optimize for the g-vector.
  res <- optim(par=G, fn=optim.fn, optimize.for="g", 
             g.old=g, thetaL.old=theta.low, thetaH.old=theta.high, 
             method="L-BFGS-B", control = list(fnscale=-1), 
             lower=rep(0, num.plots), upper=rep(1, num.plots))
  loginfo("Distance: %f", sum((G - res$par)**2))
  G <- res$par    
  loginfo("Result: %f", res$value)
  
  
  # 2. Optimize for low-risks of feature p.
  loginfo("Running lows...")
  for(p in 1:num.features) {
    len_p = feature.len[p]
    theta.old = theta.low
    
    
    # Use #levels-1 params to optimize because the logit transform is not 1-1.
    
    zetaLp.opt = optim(par=runif(len_p-1, min=zeta.lb, max=zeta.ub), 
                       fn=optim.fn, optimize.for="theta.low",
                       g.old=G, thetaL.old=theta.old, thetaH.old=theta.high, 
                       method="L-BFGS-B", control = list(fnscale=-1), 
                       lower=rep(zeta.lb, len_p-1), 
                       upper=rep(zeta.ub, len_p-1))
    # Transform back into the theta-space (simplex).
    thetaLp = logistic(c(0, zetaLp.opt$par))
    # Update the global parameter for low-risk probs.
    theta.low <<- update.theta(theta.old, p, thetaLp)
    loginfo("On feature: %2d, loglik: %f", p, zetaLp.opt$val)
    
  }
  
  
  
  print(ll(G,theta.low, theta.high, X))
  
}

run.gomMLE <- function(X=data,G0=FALSE,theta0=theta0List) {
  if(!G0){
    G0 <- rep(.5,dim(X)[1])
  }
  return (gomMLE(X,G0,theta0))  
}



