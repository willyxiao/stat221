source('functions.R')
library(MASS)

# Questions
# For implicit, what learning rate should we use? same as SGD?
# What should we use as theta.0
# Unit covariance?
# what is \lambda ? 

DIMS = 100
NSIMS = 2e4

eigen.values = c(1,1,1, rep(.02, 97))
A = diag(eigen.values)
theta.0 = rep(1, DIMS)

sigma = diag(rep(1, DIMS))
sigma.inverse = solve(sigma)
mu = rep(0, DIMS)

gen.x = function() {
  mvrnorm(n=1, mu=mu, Sigma=sigma)
}

find.risk = function(theta.t){
  t(theta.t)%*%A%*%theta.t
}

log.lik.gradient = function(x, theta){
  sigma.inverse%*%(x - theta)
}

sgd = function(t, theta.t){
  theta.t + (1 + .02*t)^(-1)*(log.lik.gradient(gen.x(), theta.t))
}

implicit = function(t, theta.t){
  (1 + (1 + .02*t)^(-1))^(-1)*(theta.t + (1 + .02*t)^(-1)*gen.x())
}

