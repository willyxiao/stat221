source("keskici_wxiao_ps2_functions.R")
source("poissonLogN_MCMC.R")

# roughly 15 seconds per simulations (ndraws=1000) ==> 2880 total simulations on 12 nodes for 1 hour

stopifnot(length(mu) == length(sigma))

runSimulation <- function(job.id, num.jobs){
  pair.num = ceiling(job.id / (num.jobs / length(mu)))
  
  num.groups = (num.jobs / length(mu))
  size.groups = (theta.NSIMS / num.groups)
  
  theta.group.num = ((job.id - 1) %% num.groups) + 1
  
  for(theta.group.offset in 1:size.groups){    
    log.theta = simThetagivenMuSigma(mu[pair.num], sigma[pair.num], J)
    Y = simYgivenTheta(exp(log.theta), w, N)
    res = poisson.logn.mcmc(Y, w, mu0=mu[pair.num], sigmasq0=sigma[pair.num]**2)

    theta.num = num.groups * (theta.group.num - 1) + theta.group.offset
    assign(getOutObjectName(pair.num, theta.num), res)
    save(list=c(getOutObjectName(pair.num, theta.num)), file=getOutFileName(pair.num,theta.num))
  }
}
