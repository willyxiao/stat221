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

rnorm.trunc = function(lower.bound, mean){
  x <- lower.bound - 1
  while(round(x) < lower.bound){
    x = rnorm(1, mean=mean, sd=3)
  }
  round(x)
}

mcmc.mh = function(y, N.start, theta.start, mcmc.niters=10000){
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
    #   Respect symmetry in (a,b)
    # lambda <- rgamma.trunc(kBound**2, shape=S+1, rate=n)
    
    # N.new = rgeom.trunc(max(y), 1/N)
    N.new = rnorm.trunc(lower.bound=max(y), mean=N.old)
    theta.new = lambda / N.new
    
#    N.new = rgeom.shift(shift=max(y), mean=max(1, N.old - lambda))
#    theta.new = lambda / N.new
    
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