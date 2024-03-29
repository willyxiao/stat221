source("keskici_wxiao_ps2_functions.R")
source("poissonLogN_MCMC.R")

TASK.NUM    <- 3
N           <- 2    # number of draws
J           <- 1000 # length of theta and w vector
theta.draws <- 9    # theta.nsims * N is the # of total simulations we'll run
Y.draws     <- 15
  
w           <- rep(1, J) # weights are all set to 1
mu          <- c(1.6, 2.5, 5.2, 4.9)
sigma       <- c(0.7, 1.3, 1.3, 1.6)

stopifnot(length(mu) == length(sigma))

simTheta <- function(pair.num, J){
  simThetagivenMuSigma(J, mu[pair.num], sigma[pair.num])
}

get.mu <- function(pair.num){
  mu[pair.num]
}

get.sigma <- function(pair.num){
  sigma[pair.num]**2
}

run.task3 <- function(job.id, num.jobs) {
  runSimulation(TASK.NUM, 
                job.id, 
                num.jobs, 
                length(mu), 
                theta.draws,
                Y.draws,
                get.mu,
                get.sigma,
                simTheta, 
                w)
}

run.test3 <- function(){
  for(i in 1:12){
    run.task3(i, 12)
  }
}
