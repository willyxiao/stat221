kBound = 150
impala = c(15, 20, 21, 23, 26)
waterbuck = c(53, 57, 66, 67, 72)

log.lik <- function(N, theta, Y) {
  # Log-likelihood of the data
  sum(dpois(Y, lambda=N*theta, log=T))
}

log.prior <- function(N, theta) {
  return (1/N)
}

log.posterior <- function(N, theta, Y) {
  log.lik(N, theta, Y) + log.prior(N, theta)
}

rgamma.trunc <- function(upper.bound, s, r) {
  # Sample from truncated gamma. 
  # TODO: Find a better implementation.
  x <- upper.bound + 10
  while(x > upper.bound) {
    x = rgamma(1, shape=s, rate=r)
  }
  return(x)
}

mcmc.mh = function(y, mcmc.iters=10000, N.start, theta.start){
  S = sum(y)
  n = length(y)
  mcmc.chain <- matrix(nrow=mcmc.niters, ncol=2)
  mcmc.chain[1,] = c(N.start, theta.start)
  nacc <- 0
  
  for(i in 2:mcmc.niters) {
    lambda <- rgamma.trunc(kBound**2, s=S+1, r=n)
    # 1. Current state
    N.old = mcmc.chain[i-1, 1]
    theta.old = mcmc.chain[i-1, 2]
    # 2. Propose new state
    #   Respect symmetry in (a,b)
    alpha.new = runif(1, min=lambda/kBound, max=kBound)
    beta.new = ldambda / alpha.new
    # 3. Ratio
    mh.ratio = min(0, log.posterior(alpha.new, beta.new, y) - 
                     log.posterior(alpha.old, beta.old, y))
    if(runif(1) < exp(mh.ratio)) {
      # Accept 
      mcmc.chain[i, ] <- c(alpha.new, beta.new)
      nacc <- nacc + 1
    } else {
      mcmc.chain[i, ] <- c(alpha.old, beta.old)
    }
  }
  # Cut the burnin period.
  print(sprintf("Acceptance ratio %.2f%%", 100 * nacc / mcmc.niters))
  plot.chain(mcmc.chain)
  return(mcmc.chain)
  
}

run.impala = function(job.id){
  starting.N = job.id * max(impala)
  mcmc.mh(starting.N, mean(impala)/starting.N, impala)
}

run.waterbuck = function(job.id){
  start.N = max(waterbuck)* (job.id - 10)
  mcmc.mh(start.N, mean(waterbuck)/start.N, waterbuck)
}

run.job = function(job.id){
  if (job.id <=10){
    run.impala(job.id)
  }
  else{
    run.waterbuck(job.id)
  }
}
