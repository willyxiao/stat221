N           <- 2    # number of draws
J           <- 1000 # length of theta and w vector
theta.NSIMS <- 2    # theta.nsims * N is the # of total simulations we'll run

w           <- rep(1, J) # weights are all set to 1
mu          <- c(1.6, 2.5, 5.2, 4.9)
sigma       <- c(0.7, 1.3, 1.3, 1.6)

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

# returns the object to be saved to *.RData given the i'th mu/sigma pair 
# and the j'th simulated theta
getObjectName <- function(i, j){
  object.name.str <- "Y.%d.%d"
  sprintf(object.name.str, i, j)
}

# same as getObjectName, but gets the file name
getFileName <- function(i,j){
  file.name.str <- "Y/Y_%d_%d.RData"
  sprintf(file.name.str, i, j)
}

