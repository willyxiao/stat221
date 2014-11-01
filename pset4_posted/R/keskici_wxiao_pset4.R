kBound = 150
impala = c(15, 20, 21, 23, 26)
waterbuck = c(53, 57, 66, 67, 72)
 
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
  mcmc.chain
}

post.proc = function(job.id, chain){
  if (job.id <=10){
    output.format = "keskici_wxiao_ps4_task_impala_plot%d.png"
    name = sprintf(output.format, job.id)
  }
  else{
    output.format = "keskici_wxiao_ps4_task_waterbuck_plot%d.png"
    name = sprintf(output.format, job.id - 10)
  }
  png(name)
  plot.chain(chain)
  dev.off()
}

niters = 1e4

run.impala = function(job.id){
  starting.N = job.id * max(impala)
  chain = mcmc.mh(impala, starting.N, mean(impala)/starting.N, niters)
  chain = chain[(0.1*niters):niters,]
  post.proc(chain)
}

run.waterbuck = function(job.id){
  start.N = max(waterbuck)* (job.id - 10)
  chain = mcmc.mh(waterbuck, start.N, mean(waterbuck)/start.N, niters)
  chain = chain[0.1*niters:niters,] #burnin period
  post.proc(chain)  
}

run.job = function(job.id){
  if (job.id <= 10){
    run.impala(job.id)
  } else{
    run.waterbuck(job.id)
  }
}

run.test = function(){
  for(i in 1:20){
    print(i)
    run.job(i)
  }
}