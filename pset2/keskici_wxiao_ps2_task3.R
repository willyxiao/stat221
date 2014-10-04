source("keskici_wxiao_ps2_functions.R")
source("poissonLogN_MCMC.R")

N           <- 2    # number of draws
J           <- 1000 # length of theta and w vector
#theta.draws <- 30    # theta.nsims * N is the # of total simulations we'll run
#Y.draws     <- 12

theta.draws <- 3
Y.draws <- 1
  
w           <- rep(1, J) # weights are all set to 1
mu          <- c(1.6, 2.5, 5.2, 4.9)
sigma       <- c(0.7, 1.3, 1.3, 1.6)

# roughly 15 seconds per simulations (ndraws=1000) ==> 2880 total simulations on 12 nodes for 1 hour

stopifnot(length(mu) == length(sigma))

runSimulation <- function(job.id, num.jobs){
  pair.num = ceiling(job.id / (num.jobs / length(mu)))
  
  num.groups = (num.jobs / length(mu))
  size.groups = (theta.draws / num.groups)
  
  theta.group.num = ((job.id - 1) %% num.groups) + 1
  
  for(theta.group.offset in 1:size.groups){
    log.theta = simThetagivenMuSigma(mu[pair.num], sigma[pair.num], J)

    for(Y.offset in 1:Y.draws){
      Y = simYgivenTheta(exp(log.theta), w, N)
      
      res = poisson.logn.mcmc(Y, w, mu0=mu[pair.num], sigmasq0=sigma[pair.num]**2)
      res$real.log.theta = log.theta

      theta.num = size.groups * (theta.group.num - 1) + theta.group.offset
      save(res, file=getOutFileName(pair.num,theta.num,Y.offset))
      gc()      
    }
  }
}
