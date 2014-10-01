
# theta is length J
# w is length J
# N is scaler
simYgivenTheta <- function(theta, w, N) {
  stopifnot(length(theta) == length(w))  
  J = length(theta)
  matrix(rpois(N*J, w*theta), nrow=J,ncol=N)
}

# all scalars
simThetaGivenMuSigma <- function(mu, sigma, J) {
  rnorm(J, mu, sigma)
}