source("keskici_wxiao_ps2_functions.R")
source("poissonLogN_MCMC.R")
source('rASL.R')

TASK.NUM = 5
WEIGHTS.FILE = 'weights.txt'

N           <- 2    # number of draws
J           <- 1000 # length of theta and w vector
#theta.draws <- 15    # theta.nsims * N is the # of total simulations we'll run
#Y.draws     <- 24

theta.draws = 9
Y.draws = 12

w           <- read.table(WEIGHTS.FILE)[,1]
x0          <- c(1.6,1.6,1.6,1.6)
m           <- c(0,-0.7,0.7,0)
b           <- c(1.3,1.3,1.3,2.6)

stopifnot(length(x0) == length(m) && length(x0) == length(b))

simTheta <- function(pair.num, J) {
  rASL(J, x0[pair.num], m[pair.num], b[pair.num])
}

get.mu <- function(pair.num){
  0
}

get.sigma <- function(pair.num){
  1
}

run.task5 <- function(job.id, num.jobs) {
  runSimulation(TASK.NUM, 
                job.id, 
                num.jobs, 
                length(x0),
                theta.draws,
                Y.draws,
                get.mu,
                get.sigma,
                simTheta, 
                w, 
                TRUE)
}

run.test5 <- function(){
  for(i in 1:12){
    run.task5(i, 12)
  }
}

