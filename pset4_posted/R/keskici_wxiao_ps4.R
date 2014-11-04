library(MASS)
library(mvtnorm)

kBound = 300
impala_bound = 300
waterbuck_bound = 500
impala = c(15, 20, 21, 23, 26)
waterbuck = c(53, 57, 66, 67, 72)
BURNIN = 0.5
RUN_NUMBER = 5
NUM_JOBS = 20
 
plot.chain <- function(mcmc.chain, is.impala=T) {
  mcmc.niters = nrow(mcmc.chain)
  if (is.impala){
    mcmc.chain = mcmc.chain[mcmc.chain[, 1] < impala_bound,]
  }else{
    mcmc.chain = mcmc.chain[mcmc.chain[, 1] < waterbuck_bound,]
  }
  f = kde2d(x=mcmc.chain[, 1], y=mcmc.chain[, 2], n=200)
  image(f, xlim=c(0, kBound), ylim=c(0, 1))
}

log.lik <- function(N, theta, Y) {
  sum(dbinom(Y, size=N, prob=theta, log=T))
}

log.prior <- function(N, theta) {
  1/N
}

log.posterior <- function(N, theta, Y) {
  log.lik(N, theta, Y) + log.prior(N, theta)
}

rgamma.trunc <- function(upper.bound, shape, rate) {
  # Sample from truncated gamma. 
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

rmvnorm.trunc = function(N, theta, lower.bound, sigma){
  res = NULL
  condition = TRUE
  while(condition){
    res = rmvnorm(1, mean=c(N, theta), sigma)
    condition = res[1] < lower.bound
    condition = condition || (res[2] < 0 || 1 < res[2])    
  }
  res
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
    lambda = rnorm(1, mean=theta.old*N.old, sd=2)
    theta.tmp = runif(1)
    N.new = min(max(max(y), round(lambda / theta.tmp)), kBound**2)
    theta.new = min(1, max(0, lambda / N.new))
    
    stopifnot(0 <= theta.new && theta.new <= 1)
    stopifnot(N.new >= max(y))
    
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
  png(name)
  
  if (job.id <=10){
    output.format = "keskici_wxiao_ps4_task_impala_run%d_plot%d.png"
    name = sprintf(output.format, RUN_NUMBER, job.id)
    plot.chain(chain, T)
    
  }
  else{
    output.format = "keskici_wxiao_ps4_task_waterbuck_run%d_plot%d.png"
    name = sprintf(output.format, RUN_NUMBER, job.id - 10)
    plot.chain(chain, F)
    
  }
  dev.off()
  output.format = "keskici_wxiao_ps4_run%d_job%d.RData"
  name = sprintf(output.format, RUN_NUMBER, job.id)
  
  save(chain, file=name)
}

run.impala = function(job.id, niters=1e5){
  starting.N = job.id * max(impala)
  chain = mcmc.mh(impala, starting.N, mean(impala)/starting.N, niters)
  chain = chain[(BURNIN*niters):niters,]
  gc()
  post.proc(job.id, chain)
  chain
}

run.waterbuck = function(job.id, niters=1e5){
  start.N = ceiling(max(waterbuck)* (job.id - 9) /2)
  chain = mcmc.mh(waterbuck, start.N, mean(waterbuck)/start.N, niters)
  chain = chain[(BURNIN*niters):niters,] #burnin period
  gc()
  post.proc(job.id, chain)
  chain
}

run.job = function(job.id, niters=1e6){
  if (job.id <= 10){
    run.impala(job.id, niters)
  } else{
    run.waterbuck(job.id, niters)
  }
}

run.test = function(){
  for(i in 1:20){
    print(i)
    run.job(i, 1e6)
  }
}
