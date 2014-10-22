rm(list=ls())
library(mvtnorm)
library(glmnet)
library(lars)
#####################################
#Functions from paper provided by Panos

# genjerry, genx2 are functions taken from the above paper.
# These functions generate the simulation data.
# NOTE: use function sample.data(..) instead.
genjerry = function(x, snr){
  # generate data according to Friedman's setup
  n=nrow(x)
  p=ncol(x)
  b=((-1)^(1:p))*exp(-2*((1:p)-1)/20)
  # b=sample(c(-0.8, -0.45, 0.45, 0.9, 1.4), size=p, replace=T)
  # ((-1)^(1:p))*(1:p)^{-0.65}#exp(-2*((1:p)-1)/20)
  f=x%*%b
  e=rnorm(n)
  k=sqrt(var(f)/(snr*var(e)))
  y=f+k*e
  return(list(y=y, beta=b))
}

genx2 = function(n,p,rho){
  #    generate x's multivariate normal with equal corr rho
  # Xi = b Z + Wi, and Z, Wi are independent normal.
  # Then Var(Xi) = b^2 + 1
  #  Cov(Xi, Xj) = b^2  and so cor(Xi, Xj) = b^2 / (1+b^2) = rho
  z=rnorm(n)
  if(abs(rho)<1){
    beta=sqrt(rho/(1-rho))
    x0=matrix(rnorm(n*p),ncol=p)
    A = matrix(z, nrow=n, ncol=p, byrow=F)
    x= beta * A + x0
  }
  if(abs(rho)==1){ x=matrix(z,nrow=n,ncol=p,byrow=F)}
  
  return(x)
}

sample.data <- function(dim.n, dim.p, rho=0.0, snr=1) {
  # Samples the dataset according to Friedman et. al.
  #
  # 1. sample covariates
  X = genx2(dim.n, dim.p, rho)
  # 2. ground truth params.
  theta = ((-1)^(1:dim.p))*exp(-2*((1:dim.p)-1)/20)
  
  f= X %*% theta
  e = rnorm(dim.n)
  k= sqrt(var(f)/(snr*var(e)))
  y=f+k*e
  return(list(y=y, X=X, theta=theta))
}

dist <- function(x, y) {
  if(length(x) != length(y))
    stop("MSE should compare vectors of same length")
  sqrt(mean((x-y)^2))
}

##############################################
#3.a
#Adaptation of Panos's function for problem

run.glmnet <- function(dim.n, dim.p, method = "naive",
                       rho.values=c(0.0, 0.1, 0.2, 0.5, 0.9, 0.95),
                       nreps=3, 
                       verbose=F) {
  stopifnot(method == "naive" || method == "cov" || method == "lars")
  ## Runs glmnet() for various param values.
  ##
  niters = 0
  cols = c("method","rho", "rep", "time", "mse")
  timings = matrix(nrow=0, ncol=length(cols))
  colnames(timings) <- cols
  rownames(timings) = NULL
  total.iters = nreps * length(rho.values)
  
  pb = txtProgressBar(style=3)
  
  seeds=sample(1:1e9, size=total.iters)
  for(i in 1:nreps) {
    for(rho in rho.values) {
      niters = niters + 1
      set.seed(seeds[niters])
      # 1. (For every repetition) Sample the dataset
      dataset = sample.data(dim.n=dim.n, dim.p=dim.p, rho=rho, snr=3.0)
      true.theta = dataset$theta
      x = dataset$X
      y = dataset$y
      stopifnot(nrow(x) == dim.n, ncol(x) == dim.p)
      # 1b. Define metrics:
      #   dt = time for the method to finish
      #   mse = Distance (e.g. RMSE) of the estimates to the ground truth.
      #         (q1, q2, q3) representing the quartiles (since glmnet returns grid of estimates)
      #         Implicit has (x, x, x) i.e., the same value in all places.
      new.dt = 0
      new.mse = NA
      # 2. Run the method.
      
      if(method == "lars"){
        new.dt = system.time({ fit = lars(x, y, type = "lasso", normalize=FALSE, use.Gram =FALSE)})[1]
      }else{
        new.dt = system.time({ fit = glmnet(x, y, alpha=1, standardize=FALSE, type.gaussian=method)})[1]
        new.mse = median(apply(fit$beta, 2, function(est) dist(est, true.theta)))
      }
      # 3. Tabulate timings
      timings = rbind(timings, c(method, rho, i, 
                                 new.dt, 
                                 new.mse))
      setTxtProgressBar(pb, niters/total.iters)
    }
    
  }
  return(timings)
}
#use to get desired data
means = function(timings, cors){
  means = rep(0, cors)
  for (i in 1:length(timings[,1])){
    bucket = ((i -1) %% cors) + 1
    means[bucket] = means[bucket] + as.numeric(timings[i,][4])
  }
  means/(length(timings[,1])/cors)
}

#Use stargzer to generate latex tables
get_table = function(c1,c2,c3, title){
  c1 = rbind(c1,c2)
  c1 = rbind(c1,c3)
  colnames(c1) <- c("0", "0.1", "0.2", "0.5", "0.9", "0.95")
  rownames(c1) <- c("glmnet (type = \"naive\")", "glmnet (type = \"cov\")", "lars")
  stargazer(c1, title = title)
}
###These output latex tables that are in our writeup
#Commenting out because the later ones take awhile to run
#get_table(means(run.glmnet(1000,100, "naive"), 6), means(run.glmnet(1000,100, "cov"),6), 
#          means(run.glmnet(1000,100, "lars"),6), "N = 1000, p = 100")

#get_table(means(run.glmnet(5000,100, "naive"), 6), means(run.glmnet(5000,100, "cov"),6), 
#          means(run.glmnet(5000,100, "lars"),6), "N = 5000, p = 100")

#get_table(means(run.glmnet(100,1000, "naive"), 6), means(run.glmnet(100,1000, "cov"),6), 
#          means(run.glmnet(100,1000, "lars"),6), "N = 100, p = 1000")

#get_table(means(run.glmnet(100,5000, "naive"), 6), means(run.glmnet(100,5000, "cov"),6), 
#          means(run.glmnet(100,5000, "lars"),6), "N = 100, p = 5000")

#get_table(means(run.glmnet(100,20000, "naive"), 6), means(run.glmnet(100,20000, "cov"),6), 
#          means(run.glmnet(100,20000, "lars"),6), "N = 100, p = 20000")

#get_table(means(run.glmnet(100,50000, "naive"), 6), means(run.glmnet(100,50000, "cov"),6), 
#          means(run.glmnet(100,50000, "lars"),6), "N = 100, p = 50000")

#3.c

lr <- function(alpha, n) {
  ## learning rate
  alpha / (alpha + n)
}


sgd_implicit <- function(data, alpha) {
  # check.data(data)
  # Implements implicit
  n = nrow(data$X)
  p = ncol(data$X)
  A = cov(data$X)
  # matrix of estimates of SGD (p x iters)
  theta.sgd = matrix(0, nrow=p, ncol=1)
  I = diag(p)
  
  for(i in 1:n) {
    xi = data$X[i, ]
    theta.old = theta.sgd[, i]
    ai = lr(alpha, i)
    
    # make computations easier.
    xi.norm = sum(xi^2)
    lpred = sum(theta.old * xi)
    fi = 1 / (1 + ai * sum(xi^2))
    yi = data$y[i]
    # Implicit SGD
    theta.new = (theta.old  - ai * fi * lpred * xi) +  
      (ai * yi * xi - ai^2 * fi * yi * xi.norm * xi)
    theta.sgd = cbind(theta.sgd, theta.new)
  }
  
  return(theta.sgd[,length(theta.sgd[1,])])
}

sgd <- function(data, alpha) {
  # check.data(data)
  # Implements implicit
  n = nrow(data$X)
  p = ncol(data$X)
  A = cov(data$X) #this is an estimate for A
  trace.A = sum(diag(A)) 
  # matrix of estimates of SGD (p x iters)
  theta.sgd = matrix(0, nrow=p, ncol=1)
  # params for the learning rate seq.
  # gamma0 = 1 / (sum(seq(0.01, 1, length.out=p)))
  #trace = sum(diag(data$A))  # NOTE: data snooping.
  I = diag(p)
  
  for(i in 1:n) {
    xi = data$X[i, ]
    theta.old = theta.sgd[, i]
    ai = alpha / (alpha * trace.A + i)
    
    # make computations easier.
    xi.norm = sum(xi^2)
    lpred = sum(theta.old * xi)
    fi = 1 / (1 + ai * sum(xi^2))
    yi = data$y[i]
    # Standard SGD
    theta.new = (theta.old - ai * lpred * xi) + ai * yi * xi
    theta.sgd = cbind(theta.sgd, theta.new)
  }
  
  return(theta.sgd[,length(theta.sgd[1,])])
}

a = sample.data(10,3,c(0.0, 0.1, 0.2, 0.5, 0.9, 0.95), 3)


