# Copyright (c) 2014
# Panos Toulis, ptoulis@fas.harvard.edu
# Helper code for pset3 for Stat 221, Fall 2014
# rm(list=ls())
#
# EXAMPLE run:
#   A = generate.A(p=10)
#   d = sample.data(dim.n=1e5, A)
#   sgd(d, alpha=2)  # should give 'biased' results
#   sgd(d, alpha=100) # better results
#
# Frequentist variance (for a quick demo):
#   run.sgd.many(nreps=500, alpha=100, nlist=c(1e3, 5e3, 1e4), p=5, verbose=T)
library(mvtnorm)

random.orthogonal <- function(p) {
  # Get an orthogonal matrix.
  B = matrix(runif(p^2), nrow=p)
  qr.Q(qr(B))
}

generate.A <- function(p) {
  return (diag(seq(0.01, 1, length.out=p)))
  
  # Create A matrix (variance of the covariates xn)
  Q = random.orthogonal(p)
  lambdas = seq(0.01, 1, length.out=p)
  A = Q %*% diag(lambdas) %*% t(Q)
  return(A)
}

sample.data <- function(dim.n, A, 
                        model="gaussian") {
  # Samples the dataset. Returns a list with (Y, X, A ,true theta)
  dim.p = nrow(A)
  # This call will make the appropriate checks on A.
  X = rmvnorm(dim.n, mean=rep(0, dim.p), sigma=A)
  theta = matrix(1, ncol=1, nrow=dim.p)
  epsilon = rnorm(dim.n, mean=0, sd=1)
  # Data generation
  y = X %*% theta  + epsilon
  
  return(list(Y=y, X=X, A=A, theta=theta))
}

check.data <- function(data) {
  # Do this to check the data object.
  # 
  nx = nrow(data$X)
  ny = length(data$Y)
  p = ncol(data$X)
  stopifnot(nx==ny, p==length(data$theta))
  lambdas = eigen(cov(data$X))$values
  print(lambdas)
  print(mean(data$Y))
  print(var(data$Y))
  print(1 + sum(cov(data$X)))
}

plot.risk <- function(data, est) {
  # est = p x niters 
  est.bias = apply(est, 2, function(colum) 
    log(t(colum-data$theta) %*% data$A %*% (colum-data$theta)))
  
#  est.bias = apply(est, 2, function(colum) 
#    t(colum-data$theta) %*% data$A %*% (colum-data$theta))
  print(sprintf("Risk of first estimate = %.3f Risk of last estimate %.3f", 
                head(est.bias, 1), tail(est.bias, 1)))
  plot(est.bias, type="l", lty=3)
}

lr <- function(alpha, n) {
  alpha / (alpha + n)
}

alg.implicit = function(ai, xi, yi, i, theta.old){
  xi.norm = sum(xi^2)
  lpred = sum(theta.old * xi)
  fi = 1 / (1 + ai * sum(xi^2))
  
  (theta.old  - ai * fi * lpred * xi) +  
    (ai * yi * xi - ai^2 * fi * yi * xi.norm * xi)
}

alg.sgd = function(ai, xi, yi, i, theta.old){
  lpred = sum(theta.old * xi)
  (theta.old - ai * lpred * xi) + ai * yi * xi      
}

alg.asgd = function(ai, xi, yi, i, theta.old){
  (1 - 1/i)*theta.old + (1/i)*alg.sgd(ai, xi, yi, i, theta.old)
}

base.method <- function(data, alpha, alg, plot=T) {
  # check.data(data)

  n = nrow(data$X)
  p = ncol(data$X)
  
  trace = sum(diag(data$A))  # NOTE: data snooping.
  I = diag(p)
  
  theta.old = rep(0, p)
  theta.new = NULL
  
  for(i in 1:n) {
    xi = data$X[i, ]
 
    ai = alpha / (alpha * trace + i)
    
    if (identical(alg, implicit)){
      ai = lr(alpha, i)
    }
    
    yi = data$Y[i]
    
    theta.new = alg(ai, xi, yi, i, theta.old)
    theta.old = theta.new
  }

  theta.new
}

sgd = function(data, alpha, plot=T){
  base.method(data, alpha, alg.sgd, plot)
}

asgd = function(data, alpha, plot=T){
  base.method(data, alpha, alg.asgd, plot)  
}

implicit = function(data, alpha, plot=T){
  base.method(data, alpha, alg.implicit, plot)
}

sqrt.norm <- function(X) {
  sqrt(mean(X^2)  )
}

# Replicate this code on Odyssey
# to get the frequentist variance of SGD.
run.alg.many = function(nreps, 
                        alpha, 
                        alg, 
                        nlist=as.integer(seq(100, 1e5, length.out=10)), 
                        p=100, 
                        verbose=F) {
  ## Run SGD (implicit) with multiple n(=SGD iterations)
  #
  # Example:
  #   run.sgd.many(nreps=1000, alpha=2, p=4)
  #
  # Plots || Empirical variance - Theoretical ||
  #
  A = generate.A(p)
  dist.list = c()  # vector of || . || distances
  bias.list = c()
  for(n in nlist) {
    # matrix of the last iterate theta_n
    last.theta = matrix(NA, nrow=0, ncol=p)
    # Compute theoretical variance
    data0 = sample.data(n, A)
    I = diag(p)
    Sigma.theoretical <- alpha * solve(2 * alpha * A - I) %*% A
    stopifnot(all(eigen(Sigma.theoretical)$values > 0))
    
    # Get many replications for each n
    for(j in 1:nreps) {
      data = sample.data(n, A)
      out = alg(data, alpha, plot=F)
      # Every replication stores the last theta_n
      last.theta <- rbind(last.theta, out)
    }

    theta.avg = colMeans(last.theta)
    print(sprintf("n = %d", n))
    print("Avg. estimates for theta*")
    print(theta.avg)
    print("Bias")
    bias = sqrt.norm(theta.avg - data$theta)
    bias.list = c(bias.list, bias)
    print(bias.list)
    # Store the distance.
    empirical.var = (1 / lr(alpha, n)) * cov(last.theta)
    dist.list <- c(dist.list, sqrt.norm(empirical.var - Sigma.theoretical))
    print("Vector of ||empirical var.  - theoretical var.||")
    print(dist.list)
  }
  list(dist=dist.list, bias=bias.list)
}
