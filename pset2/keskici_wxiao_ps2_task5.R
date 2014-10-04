source("keskici_wxiao_ps2_functions.R")
source("poissonLogN_MCMC.R")
source('rASL.R')

TASK.NUM = 5
WEIGHTS.FILE = 'weights.txt'

N           <- 2    # number of draws
J           <- 1000 # length of theta and w vector
theta.draws <- 3    # theta.nsims * N is the # of total simulations we'll run
Y.draws     <- 2

w           <- read.table(WEIGHTS.FILE)[,1]
x0          <- c(1.6,1.6,1.6,1.6)
m           <- c(0,-0.7,0.7,0)
b           <- c(1.3,1.3,1.3,2.6)


# roughly 15 seconds per simulations (ndraws=1000) ==> 2880 total simulations on 12 nodes for 1 hour

stopifnot(length(x0) == length(m) && length(x0) == length(b))

runSimulation <- function(job.id, num.jobs){
  pair.num = ceiling(job.id / (num.jobs / length(x0)))
  
  num.groups = (num.jobs / length(x0))
  size.groups = (theta.draws / num.groups)
  
  theta.group.num = ((job.id - 1) %% num.groups) + 1
  
  for(theta.group.offset in 1:size.groups){
    log.theta = rASL(J, x0[pair.num], m[pair.num], b[pair.num])

    for(Y.offset in 1:Y.draws){
      theta.num = size.groups * (theta.group.num - 1) + theta.group.offset
      print(paste("Writing", getOutFileName(pair.num, theta.num, Y.offset, TASK.NUM)))

      Y = simYgivenTheta(exp(log.theta), w, N)
      res = poisson.logn.mcmc(Y, w)
      res$real.log.theta = log.theta

      save(res, file=getOutFileName(pair.num,theta.num,Y.offset, TASK.NUM))
      gc()      
    }
  }
}
