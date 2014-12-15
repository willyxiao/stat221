source('functions.R')
library(MASS)

DIMS = 100
NSIMS = 1e4

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

sgd.learn.rate = function(t){
  (1 + .02*t)^(-1)
}

sgd.base = function(t, learn.rate, theta.t){
  theta.t + learn.rate(t)*(log.lik.gradient(gen.x(), theta.t))
}

sgd = function(t, theta.t){
  sgd.base(t, sgd.learn.rate, theta.t)
}

asgd = function(t, theta.t){
  (1 - 1/t)*theta.t + (1/t)*sgd.base(t, function(t) (1+.02*t)^(-2/3), theta.t)
}

implicit = function(t, theta.t){
  (1 + sgd.learn.rate(t))^(-1)*(theta.t + sgd.learn.rate(t)*gen.x())
}

batch = function(t, theta.t){
  (1 - 1/t)*theta.t + (1/t)*gen.x()  
}

asgd.bad = function(t, theta.t){
  (1 - 1/t)*theta.t + (1/t)*sgd.base(t, function(t) (1+t)^(-1/2), theta.t)
}

plot.all = function() {
  x = log10(1:NSIMS)

  plot(x, log10(sim.alg(sgd, theta.0, NSIMS)), "l", ylim = c(-2.5, 1), col="blue", xaxt='n', ann=FALSE)
  lines(x, log10(sim.alg(asgd, theta.0, NSIMS)), col="green")
  lines(x, log10(sim.alg(batch, theta.0, NSIMS)), col="yellow", lwd=2)
  lines(x, log10(sim.alg(asgd.bad, theta.0, NSIMS)), col="red", lwd=2)
  lines(x, log10(sim.alg(implicit, theta.0, NSIMS)))
  
  axis(1, at=0:log10(NSIMS), lab=0:log10(NSIMS))
  title(main="question 2(a)", xlab="Data Points (1e^x)", ylab="Excess Risk (1e^y)")
  legend('bottomleft', legend=c("SGD", "ASGD", "ASGD_BAD", "Implicit", "Batch"), lty=c(1,1,1,1,1,1), col=c('blue', 'green', 'red', 'black', 'yellow'))
}
