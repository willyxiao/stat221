source('functions.R')

DIMS = 100
NSIMS = 1e5

gamma0 = 1/sum(seq(.01, 1, length.out=100))
lambda0 = .01
theta.0 = rep(0, DIMS)

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

A = generate.A(DIMS)

gen.data = function(A) {
  x = rmvnorm(1, mean=rep(0, DIMS), sigma=A)
  y = rnorm(1, sum(x), 1)
  list(y=y, x=x)
}

find.risk = function(theta.t){
    t(theta.t - rep(1, DIMS))%*%A%*%(theta.t - rep(1, DIMS))
}

sgd.base = function(t, ai, theta.t){
  data = gen.data(A)
#  2*ai*as.numeric(data$y - data$x%*%theta.t)*t(data$x)
  lpred = sum(theta.t * data$x)
  as.vector((theta.t - ai * lpred * data$x) + ai * data$y * data$x)
}

#2*gamma(y-x*theta)*x

sgd = function(t, theta.t){
#  ai = lambda0 / (lambda0 / gamma0 + t)# * (1 + gamma0 * lambda0 * t)^(-1)  
  ai = gamma0 * (1 + gamma0*lambda0*t)^(-1)
  sgd.base(t, ai, theta.t)  
}

asgd = function(t, theta.t){
  ai = gamma0 * (1 + gamma0 * lambda0 * t)^(-1/3)
  (1 - 1/t)*theta.t + (1/t)*sgd.base(t, ai, theta.t)
}

plot.all = function() {
  x = log10(1:NSIMS)
  
  plot(x, (sim.alg(sgd, theta.0, NSIMS)), "l", col="blue", xaxt='n', ann=FALSE)
  lines(x, (sim.alg(asgd, theta.0, NSIMS)), col="green")
#  lines(x, log10(sim.alg(batch, theta.0, NSIMS)), col="yellow", lwd=2)
#  lines(x, log10(sim.alg(asgd.bad, theta.0, NSIMS)), col="red", lwd=2)
#  lines(x, log10(sim.alg(implicit, theta.0, NSIMS)))
  
  axis(1, at=0:log10(NSIMS), lab=0:log10(NSIMS))
  title(main="question 2(a)", xlab="Data Points (1e^x)", ylab="Excess Risk (1e^y)")
  legend('bottomleft', legend=c("SGD", "ASGD", "ASGD_BAD", "Implicit", "Batch"), lty=c(1,1,1,1,1,1), col=c('blue', 'green', 'red', 'black', 'yellow'))
}

# find.risk = function(theta.t, data){
#   stopifnot(length(theta.t) == ncol(data$X))
#   t(theta.t - data$theta) %*% data$A %*% (theta.t - data$theta)
# }
# 
# batch = function(data, plot=T){
#   n = nrow(data$X)
#   p = ncol(data$X)
#   
#   theta.batch = matrix(1, nrow=p, ncol=1)
#   risk.list = NULL
# 
#   i = 1
#   while(i <= n){
#     print(i)
#     if( i < p ){      
#       i = i * 10
#       next
#     }
#     x = data$X[1:i,]
#     y = data$Y[1:i,]
#     theta.new = lm(y ~ x + 0)$coefficients
#     risk.list = c(risk.list, find.risk(theta.new, data))
#     theta.batch = cbind(theta.batch, theta.new)
#     i = i * 10
#   }
#   
#   theta.batch
#   risk.list
# }
# 
# implicit <- function(data) {
#   # check.data(data)
#   n = nrow(data$X)
#   p = ncol(data$X)
#   I = diag(p)
#   # matrix of estimates of SGD (p x iters)
#   theta.implicit = matrix(0, nrow=p, ncol=1)
#   # params for the learning rate seq.
#   gamma0 = 1 / (sum(seq(0.01, 1, length.out=p)))
#   lambda0 = 0.01
# 
#   theta.old = rep(0, p)
#   theta.new = NULL
#   
#   risk.list = NULL
#   
#   for(i in 1:n) {
#     xi = data$X[i, ]
#     ai = gamma0 / (1 + gamma0 * lambda0 * i)
#     # make computations easier.
#     theta.new = solve(I + ai*xi%*%t(xi))%*%(theta.old + ai*data$Y[i]*xi)
#     risk.list = c(risk.list, find.risk(theta.new, data))
#     theta.old = theta.new
#   }
#   
#   risk.list
# }
# 
# sgd.update = function(ai, xi, yi, i, theta.old){
# 
#   lpred = sum(theta.old * xi)
#   (theta.old - ai * lpred * xi) + ai * yi * xi
# }
# 
# # #sgd <- function(data) {
# # #  # check.data(data)
# # #  n = nrow(data$X)
# #   p = ncol(data$X)
# #   I = diag(p)
# #   # matrix of estimates of SGD (p x iters)
# #   theta.sgd = matrix(0, nrow=p, ncol=1)
# #   # params for the learning rate seq.
# #   gamma0 = 1 / (sum(seq(0.01, 1, length.out=p)))
# #   lambda0 = 0.01
# #   
# #   
# #   theta.old = rep(0, p)
# #   theta.new = NULL
# #   
# #   risk.list = NULL
# #   
# #   for(i in 1:n) {
# #     xi = data$X[i, ]
# #     ai = gamma0 / (1 + gamma0 * lambda0 * i)
# #     # make computations easier.
# #     theta.new = sgd.update(ai, xi, data$Y[i], i, theta.old)
# #     risk.list = c(risk.list, find.risk(theta.new, data))
# #     theta.old = theta.new
# #   }
# #   
# #   risk.list
# # }
# 
# asgd <- function(data) {
#   # check.data(data)
#   n = nrow(data$X)
#   p = ncol(data$X)
#   I = diag(p)
#   # matrix of estimates of SGD (p x iters)
#   #theta.asgd = matrix(0, nrow=p, ncol=1)
#   # params for the learning rate seq.
#   gamma0 = 1 / (sum(seq(.01, 1, length.out=p)))
#   lambda0 = .01
#   
#   theta.old = rep(0, p)
#   theta.new = NULL
#   
#   risk.list = NULL
#   
#   start.avg = FALSE
#   
#   for(i in 1:n) {
#     xi = data$X[i, ]
#     ai = gamma0 * (1 + gamma0 * lambda0 * i)^(-2/3)
#     stopifnot(ai > (gamma0 * (1 + gamma0 * lambda0 * i)^(-1)))
#     # make computations easier.
#     theta.sgd.new = sgd.update(ai, xi, data$Y[i], i, theta.old)
#    
# #    if(!start.avg && i > 1000){
# #      theta.new = (.99)*theta.old + (.01)*theta.sgd.new
# #      if(find.risk(theta.new, data) < find.risk(theta.sgd.new, data)){
# #        start.avg = TRUE
# #        print(i)
# #      }
# #    } else {
#       theta.new = (1 - 1/i)*theta.old + (1/i)*theta.sgd.new
# #    }
#   
#     risk.list = c(risk.list, find.risk(theta.new, data))
#     theta.old = theta.new
#   }
#   
#   risk.list
# }
# 
# plot.all = function() {
#   A = generate.A(DIMS)
#   d = sample.data(NSIMS, A)
#   
#   x = log10(1:NSIMS)
#   x.batch = 2:log10(NSIMS)
#   
#   plot(x, log10(sgd(d)), "l", ylim = c(-2.5, 1), col="blue", xaxt='n', ann=FALSE)
#   lines(x, log10(asgd(d)), col="green")
#   lines(x.batch, log10(batch(d)), col="yellow", lwd=2)
# #  lines(x, log10(implicit(d)))
#   
#   axis(1, at=0:log10(NSIMS), lab=0:log10(NSIMS))
#   title(main="question 2(b)", xlab="Data Points (1e^x)", ylab="Excess Risk (1e^y)")
#   legend('bottomleft', legend=c("SGD", "ASGD", "Implicit", "Batch"), lty=c(1,1,1,1,1), col=c('blue', 'green', 'black', 'yellow'))
# }
