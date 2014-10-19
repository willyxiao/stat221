source('functions.R')

DIMS = 100
NSIMS = 1e4

random.orthogonal <- function(p) {
  # Get an orthogonal matrix.
  B = matrix(runif(p^2), nrow=p)
  qr.Q(qr(B))
}

generate.A <- function(p) {
  # Create A matrix (variance of the covariates xn)
  Q = random.orthogonal(p)
  lambdas = seq(0.01, 1, length.out=p)
  A = Q %*% diag(lambdas) %*% t(Q)
  return(A)
}

find.risk = function(theta.t){
  t(theta.t - rep(1, DIMS))%*%A%*%(theta.t - rep(1, DIMS))
}

batch = function(data, plot=T){
  n = nrow(data$X)
  p = ncol(data$X)
  
  theta.batch = matrix(1, nrow=p, ncol=1)
  risk.list = NULL

  i = 1
  while(i <= n){
    print(i)
    if( i < p ){      
      i = i * 10
      next
    }
    x = data$X[1:i,]
    y = data$Y[1:i,]
    theta.new = lm(y ~ x + 0)$coefficients
    risk.list = c(risk.list, find.risk(theta.new))
    theta.batch = cbind(theta.batch, theta.new)
    i = i * 10
  }
  
  theta.batch
  risk.list
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
  
  risk.list = NULL
  
  for(i in 1:n) {
    xi = data$X[i, ]
    theta.old = theta.implicit[, i]
    ai = gamma0 / (1 + gamma0 * lambda0 * i)
    # make computations easier.
    theta.new = solve(I + ai*xi%*%t(xi))%*%(theta.old + ai*data$Y[i]*xi)
    risk.list = c(risk.list, find.risk(theta.new))
    theta.implicit = cbind(theta.implicit, theta.new)
  }
  
  theta.implicit
  risk.list
}

sgd <- function(data) {
  # check.data(data)
  n = nrow(data$X)
  p = ncol(data$X)
  I = diag(p)
  # matrix of estimates of SGD (p x iters)
  theta.sgd = matrix(0, nrow=p, ncol=1)
  # params for the learning rate seq.
  gamma0 = 1 / (sum(seq(0.01, 1, length.out=p)))
  lambda0 = 0.01
  
  risk.list = NULL
  
  for(i in 1:n) {
    xi = data$X[i, ]
    theta.old = theta.sgd[, i]
    ai = gamma0 / (1 + gamma0 * lambda0 * i)
    # make computations easier.
    lpred = sum(theta.old * xi)
    theta.new = (theta.old - ai * lpred * xi) + ai * data$Y[i] * xi
    risk.list = c(risk.list, find.risk(theta.new))
    theta.sgd = cbind(theta.sgd, theta.new)
  }
  
  theta.sgd
  risk.list
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
  
  risk.list = NULL
  for(i in 1:n) {
    xi = data$X[i, ]
    theta.old = theta.asgd[, i]
    ai = (gamma0 / (1 + gamma0 * lambda0 * i))^(2/3)
    # make computations easier.
    lpred = sum(theta.old * xi)
    theta.new = (1 - 1/i)*theta.old + (1/i)*((theta.old - ai * lpred * xi) + ai * data$Y[i] * xi)
    risk.list = c(risk.list, find.risk(theta.new))
    theta.asgd = cbind(theta.asgd, theta.new)
  }
  
  theta.asgd
  risk.list
}

plot.all = function() {
  A = generate.A(DIMS)
  d = sample.data(NSIMS, A)
  
  x = log10(1:NSIMS)
  x.batch = 2:log10(NSIMS)
  
  plot(x, log10(sgd(d)), "l", ylim = c(-2.5, 1), col="blue", xaxt='n', ann=FALSE)
  lines(x, log10(asgd(d)), col="green")
  lines(x.batch, log10(batch(d)), col="yellow", lwd=2)
  lines(x, log10(implicit(d)))
  
  axis(1, at=0:log10(NSIMS), lab=0:log10(NSIMS))
  title(main="question 2(b)", xlab="Data Points (1e^x)", ylab="Excess Risk (1e^y)")
  legend('bottomleft', legend=c("SGD", "ASGD", "Implicit", "Batch"), lty=c(1,1,1,1,1), col=c('blue', 'green', 'black', 'yellow'))
}
