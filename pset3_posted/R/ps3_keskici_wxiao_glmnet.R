rm(list=ls())
library(mvtnorm)
library(glmnet)
library(lars)
library(stargazer)
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
means = function(timings, cors=6){
  means = rep(0, cors)
  for (i in 1:length(timings[,1])){
    bucket = ((i -1) %% cors) + 1
    means[bucket] = means[bucket] + as.numeric(timings[i,][4])
  }
  means/(length(timings[,1])/cors)
}

mses = function(timings, cors=6){
  mses = rep(0, cors)
  for (i in 1:length(timings[,1])){
    bucket = ((i -1) %% cors) + 1
    mses[bucket] = mses[bucket] + as.numeric(timings[i,][5])
  }
  mses/(length(timings[,1])/cors)
}

#Use stargzer to generate latex tables for glmnet
get_table = function(c1,c2,c3, title, lars=TRUE){
  c1 = rbind(c1,c2)
  colnames(c1) <- c("0", "0.1", "0.2", "0.5", "0.9", "0.95")
  if(lars){
    c1 = rbind(c1,c3)
    rownames(c1) <- c("glmnet (type = \"naive\")", "glmnet (type = \"cov\")", "lars")
  }else{
    rownames(c1) <- c("glmnet (type = \"naive\")", "glmnet (type = \"cov\")")
  }
  stargazer(c1, title = title)
}
###These output latex tables that are in our writeup. To see output in normal form just run ai, bi, or ci
#Commenting out because the later ones take awhile to run

a1 = run.glmnet(1000,100, "naive")
b1 = run.glmnet(1000,100, "cov")
c1 = run.glmnet(1000,100, "lars")
get_table(means(a1), means(b1), means(c1), "Times: N = 1000, p = 100")
get_table(mses(a1), mses(b1), mses(c1), "MSEs: N = 1000, p = 100", FALSE)


a2 = run.glmnet(5000,100, "naive")
b2 = run.glmnet(5000,100, "cov")
c2 = run.glmnet(5000,100, "lars")
get_table(means(a2), means(b2), means(c2), "Times: N = 5000, p = 100")
get_table(mses(a2), mses(b2), mses(c2), "MSEs: N = 5000, p = 100", FALSE)

a3 = run.glmnet(100,1000, "naive")
b3 = run.glmnet(100,1000, "cov")
c3 = run.glmnet(100,1000, "lars")
get_table(means(a3), means(b3), means(c3), "Times: N = 100, p = 1000")
get_table(mses(a3), mses(b3), mses(c3), "MSEs: N = 100, p = 1000", FALSE)

#Commenting these ones out because they take a bit to run (in case we want to run the script through)

#a4 = run.glmnet(100,5000, "naive")
#b4 = run.glmnet(100,5000, "cov")
#c4 = run.glmnet(100,5000, "lars")
#get_table(means(a4), means(b4), means(c4), "Times: N = 100, p = 5000")
#get_table(mses(a4), mses(b4), mses(c4), "MSEs: N = 100, p = 5000", FALSE)

#a5 = run.glmnet(100,20000, "naive")
#b5 = run.glmnet(100,20000, "cov")
#c5 = run.glmnet(100,20000, "lars")
#get_table(means(a5), means(b5), means(c5), "Times: N = 100, p = 20000")
#get_table(mses(a5), mses(b5), mses(c5), "MSEs: N = 100, p = 20000", FALSE)

#a6 = run.glmnet(100,50000, "naive")
#b6 = run.glmnet(100,50000, "cov")
#c6 = run.glmnet(100,50000, "lars")
#get_table(means(a6), means(b6), means(c6), "Times: N = 100, p = 50000")
#get_table(mses(a6), mses(b6), mses(c6), "MSEs: N = 100, p = 50000", FALSE)

##########################
#3.c

sgd_implicit <- function(data, alpha) {
  # check.data(data)
  # Implements implicit
  n = nrow(data$X)
  p = ncol(data$X)
  A = cov(data$X)
  eigen.values = eigen(A)$values
  alpha = min(eigen.values)
  trace.A = sum(eigen.values)
  
  # matrix of estimates of SGD (p x iters)
  theta.sgd = matrix(0, nrow=p, ncol=1)
  
  for(i in 1:n) {
    xi = data$X[i, ]
    theta.old = theta.sgd[, i]
    #ai = lr(alpha, i)
    ai = alpha / (alpha * trace.A + i)
    
    
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

sgd <- function(data) {
  # check.data(data)
  # Implements implicit
  n = nrow(data$X)
  p = ncol(data$X)
  A = cov(data$X) #this is an estimate for A
  eigen.values = eigen(A)$values
  alpha = min(eigen.values)
  trace.A = sum(eigen.values)
  # matrix of estimates of SGD (p x iters)
  theta.sgd = matrix(0, nrow=p, ncol=1)
  # params for the learning rate seq.
  
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

run.sgd <- function(dim.n, dim.p, method = "implicit",
                       rho.values=c(0.0, 0.1, 0.2, 0.5, 0.9, 0.95),
                       nreps=3, 
                       verbose=F) {
  ## Runs sgd()/sgd_implicit for various param values.
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

      new.dt = 0
      new.mse = NA
      
      if(method == "implicit"){
        new.dt = system.time({ result = sgd_implicit(dataset)})
        new.mse = dist(result, true.theta)
      }else{
        new.dt = system.time({ result = sgd(dataset)})
        new.mse = dist(result, true.theta)
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

####################################
#Produce tables for 3.c
get_table_sgd = function(c1,c2, c3, c4, title){
  c1 = rbind(c1,c2)
  c1 = rbind(c1,c3)
  c1 = rbind(c1, c4)
  colnames(c1) <- c("0", "0.1", "0.2", "0.5", "0.9", "0.95")
  rownames(c1) <- c("SGD Time", "Implicit SGD Time", "SGD MSE", "Implicit SGD MSE")
  stargazer(c1, title = title)
}

sgd1 = run.sgd(1000, 100, "Normal SGD")
sgd_i1 = run.sgd(1000, 100, "implicit")
get_table_sgd(means(sgd1), means(sgd_i1), mses(sgd1), mses(sgd_i1), "N = 1000, p = 100")
#Commenting these out because they take awhile to run
#sgd2 = run.sgd(5000, 100, "Normal SGD")
#sgd_i2 = run.sgd(5000, 100, "implicit")
#get_table_sgd(means(sgd2), means(sgd_i2), mses(sgd2), mses(sgd_i2), "N = 5000, p = 100")

#sgd3 = run.sgd(100, 1000, "Normal SGD")
#sgd_i3 = run.sgd(100, 1000, "implicit")
#get_table_sgd(means(sgd3), means(sgd_i3), mses(sgd3), mses(sgd_i3), "N = 100, p = 1000")

#sgd4 = run.sgd(100, 5000, "Normal SGD")
#sgd_i4 = run.sgd(100, 5000, "implicit")
#get_table_sgd(means(sgd4), means(sgd_i4), mses(sgd4), mses(sgd_i4), "N = 100, p = 5000")

#sgd5 = run.sgd(100, 20000, "Normal SGD")
#sgd_i5 = run.sgd(100, 20000, "implicit")
#get_table_sgd(means(sgd5), means(sgd_i5), mses(sgd5), mses(sgd_i5), "N = 100, p = 20000")

#sgd6 = run.sgd(100, 50000, "Normal SGD")
#sgd_i6 = run.sgd(100, 50000, "implicit")
#get_table_sgd(means(sgd6), means(sgd_i6), mses(sgd6), mses(sgd_i6), "N = 100, p = 50000")

#3.d
#Code for Part D. SGD parts would take too long to run so it's commented out.
# Really we'd want to make this a SLURM job.
a7 = run.glmnet(100000,1000, "naive")
b7 = run.glmnet(100000,1000, "cov")
get_table(means(a7), means(b7), means(c6), "Times: N = 100000, p = 1000", FALSE)
get_table(mses(a7), mses(b7), mses(c6), "MSEs: N = 100000, p = 1000", FALSE)

#sgd7 = run.sgd(100000, 1000, "Normal SGD")
#sgd_i7 = run.sgd(100000, 1000, "implicit")
#get_table_sgd(means(sgd7), means(sgd_i7), mses(sgd7), mses(sgd_i7), "N = 100000, p = 1000")

