source('functions.R')
library(MASS)

DIMS = 100

A = diag(1:100/100)
mu = rep(0, DIMS)
alpha.0 = 1/sum(diag(A))

gen.x = function(){
  mvrnorm(1, mu, A)
}

std.learn.rate = function(t){
  alpha.0 / t
}

find.risk = function(theta.t){
  # TODO
}

