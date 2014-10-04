N           <- 2    # number of draws
J           <- 1000 # length of theta and w vector
theta.draws <- 30    # theta.nsims * N is the # of total simulations we'll run
Y.draws     <- 12

w           <- rep(1, J) # weights are all set to 1
mu          <- c(1.6, 2.5, 5.2, 4.9)
sigma       <- c(0.7, 1.3, 1.3, 1.6)

OUT_FILE_LOCATION = ""

# theta is length J
# w is length J
# N is scaler
simYgivenTheta <- function(theta, w, N) {
  stopifnot(length(theta) == length(w))  
  J = length(theta)
  matrix(rpois(N*J, w*theta), nrow=J,ncol=N)
}

# all scalars
simThetagivenMuSigma <- function(mu, sigma, J) {
  rnorm(J, mu, sigma)
}

# same as getFileName but in the out file and called "out"
getOutFileName <- function(i,j, k){
  file.name.str <- "out/out_%d_%d_%d.RData"
  sprintf(file.name.str, i, j, k)
}

