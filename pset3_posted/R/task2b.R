source('distro2b.R')
source('functions.R')

DIMS = 100

# need batch ??
batch = function(data, plot=T){  
  n = nrow(data$X)    
  p = ncol(data$X)

  theta.batch = matrix(0, nrow=p, ncol=1)
  for(i in 1:n){
    left = matrix(0, nrow=p, ncol=p)
    right = rep(0, p)
    for(j in 1:i){
      xi = data$X[j,]
      yi = data$Y[j,]
      left = left + (xi)%*%t(xi)
      right = right + xi*yi
    }
    theta.new = solve(left)*right
    theta.batch = cbind(theta.batch, theta.new)    
  }
  
  theta.batch
}

implicit <- function(data) {
  # check.data(data)
  n = nrow(data$X)
  p = ncol(data$X)
  I = diag(p)
  # matrix of estimates of SGD (p x iters)
  theta.implicit = matrix(0, nrow=p, ncol=1)
  # params for the learning rate seq.
  gamma0 = 1 / (sum(seq(0.01, 1, length.out=p)))
  lambda0 = 0.01
  
  for(i in 1:n) {
    xi = data$X[i, ]
    theta.old = theta.implicit[, i]
    ai = gamma0 / (1 + gamma0 * lambda0 * i)
    # make computations easier.
    theta.new = solve(I + ai*xi%*%t(xi))%*%(theta.old + ai*data$Y[i]*xi)
    theta.implicit = cbind(theta.implicit, theta.new)
  }
  
  theta.implicit
}

asgd <- function(data) {
  # check.data(data)
  n = nrow(data$X)
  p = ncol(data$X)
  I = diag(p)
  # matrix of estimates of SGD (p x iters)
  theta.asgd = matrix(0, nrow=p, ncol=1)
  # params for the learning rate seq.
  gamma0 = 1 / (sum(seq(0.01, 1, length.out=p)))
  lambda0 = 0.01
  
  for(i in 1:n) {
    xi = data$X[i, ]
    theta.old = theta.asgd[, i]
    ai = gamma0 / (1 + gamma0 * lambda0 * i)
    # make computations easier.
    lpred = sum(theta.old * xi)
    theta.new = (1 - 1/i)*theta.old + (1/i)*((theta.old - ai * lpred * xi) + ai * data$Y[i] * xi)
    theta.asgd = cbind(theta.asgd, theta.new)
  }
  
  theta.asgd
}

