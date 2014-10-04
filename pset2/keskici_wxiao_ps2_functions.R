OUT_FILE_LOCATION = ""

simYgivenTheta <- function(theta, w, N) {
  stopifnot(length(theta) == length(w))  
  J = length(theta)
  matrix(rpois(N*J, w*theta), nrow=J,ncol=N)
}

simThetagivenMuSigma <- function(mu, sigma, J) {
  rnorm(J, mu, sigma)
}

# same as getFileName but in the out file and called "out"
getOutFileName <- function(i,j, k){
  file.name.str <- "out/out_%d_%d_%d.RData"
  sprintf(file.name.str, i, j, k)
}

