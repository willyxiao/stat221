#simulate data
simulate.data = function(lambda, w, n){
  m = length(lambda)
  #helper function to determine which lambda
  #to use
  data = c()
  for(i in 1:n){
    draw = rpois(1, sum(rmultinom(1, m, w) * lambda))
    data = c(data, draw)
  }
  data
}