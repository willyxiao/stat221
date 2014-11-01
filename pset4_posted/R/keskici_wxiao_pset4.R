kBound = 150
 
plot.chain <- function(mcmc.chain) {
  mcmc.niters = nrow(mcmc.chain)
  burnin = 0.1 * mcmc.niters
  mcmc.chain = mcmc.chain[burnin:mcmc.niters, ]
  f = kde2d(x=mcmc.chain[, 1], y=mcmc.chain[, 2], n=100)
  image(f, xlim=c(0, kBound), ylim=c(0, 1))
}

log.lik <- function(N, theta, Y) {
  sum(dpois(Y, lambda=N*theta, log=T))
}

log.prior <- function(N, theta) {
  1/N
}

log.posterior <- function(N, theta, Y) {
  log.lik(N, theta, Y) + log.prior(N, theta)
}

rgamma.trunc <- function(upper.bound, shape, rate) {
  # Sample from truncated gamma. 
  # TODO: Find a better implementation.
  x <- upper.bound + 10
  while(x > upper.bound) {
    x = rgamma(1, shape=shape, rate=rate)
  }
  x
}

rnorm.trunc = function(mean, sd=1, lower.bound, upper.bound=NA){
  x <- lower.bound - 1
  condition = TRUE
  while(condition){
    x = rnorm(1, mean=mean, sd=sd)
    condition = x < lower.bound
    if(!is.na(upper.bound)){
      condition = condition || upper.bound < x
    }
  }
  x
}

mcmc.mh = function(y, N.start, theta.start, mcmc.niters=1e5){
  S = sum(y)
  n = length(y)
  mcmc.chain <- matrix(0, nrow=mcmc.niters, ncol=2)
  mcmc.chain[1,] = c(N.start, theta.start)
  nacc <- 0
  
  for(i in 2:mcmc.niters) {

    # 1. Current state
    N.old = mcmc.chain[i-1, 1]
    theta.old = mcmc.chain[i-1, 2]
    
    # 2. Propose new state
    N.new = round(rnorm.trunc(mean=N.old, sd=2, lower.bound=max(y)))
    theta.new = rnorm.trunc(mean=theta.old, sd=.05, lower.bound=0, upper.bound=1)
    
    # 3. Ratio
    mh.ratio = min(0, log.posterior(N.new, theta.new, y) - 
                     log.posterior(N.old, theta.old, y))

    if(runif(1) < exp(mh.ratio)) {
      # Accept 
      mcmc.chain[i, ] <- c(N.new, theta.new)
      nacc <- nacc + 1
    } else {
      mcmc.chain[i, ] <- c(N.old, theta.old)
    }
  }
 
  # Cut the burnin period.
  print(sprintf("Acceptance ratio %.2f%%", 100 * nacc / mcmc.niters))
  plot.chain(mcmc.chain)
  mcmc.chain
}