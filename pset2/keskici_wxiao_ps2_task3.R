source("keskici_wxiao_ps2_functions.R")

N           <- 2    # number of draws
J           <- 1000 # length of theta and w vector
theta.NSIMS <- 2    # theta.nsims * N is the # of total simulations we'll run

w           <- rep(1, J) # weights are all set to 1
mu          <- c(1.6, 2.5, 5.2, 4.9)
sigma       <- c(0.7, 1.3, 1.3, 1.6)

object.name.str <- "Y.%d.%d"
file.name.str <- "Y/Y_%d_%d"

for(i in 1:length(mu)){
  for(j in 1:theta.NSIMS) {
    log.theta = simThetagivenMuSigma(mu[i], sigma[i], J)
    Y = simYgivenTheta(exp(log.theta), w, N)
    assign(sprintf(object.name.str, i, j), Y)
    save(list=c(sprintf(object.name.str, i, j)), file=sprintf(file.name.str, i, j))
  }
}
