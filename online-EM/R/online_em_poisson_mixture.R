#simulate data
simulate.data = function(lambda, w, n){
  m = length(lambda)
  #helper function to determine which lambda
  #to use
  data = c()
  for(i in 1:n){
    draw = rpois(1, sum(rmultinom(1, 1, w) * lambda))
    data = c(data, draw)
  }
  data
}

w.bar = function(y, w, lambda, j){
  num = w[j] * lambda[j]**y * exp(-lambda[j])
  den = 0
  for(i in 1:length(lambda)){
      den = den + w[i] *(lambda[i]**y)*exp(-lambda[i])
  }
  #den = sum(w*(lambda**y)*exp(-lambda))
  num/ den
}

poisson.em = function(data, m){
  #init values
  w = rep(1/m, m)
  lambda = rep(mean(data), m)
  s.hat = cbind(w, lambda)
  
  for(i in 1: length(data)){
    for(j in 1:m){
      w.bar.iter = w.bar(data[i], w, lambda, j)
      if(i == 1){
      print(w.bar.iter)
      }
      s.hat[j,] = s.hat[j,] + (1/i) *(c(w.bar.iter, w.bar.iter*data[i]) - s.hat[j,])
      #print(s.hat[j,])
      w[j] = s.hat[j,1]
      lambda[j] = s.hat[j, 2] / s.hat[j, 1]
    }
  }
  c(lambda, w)
}

#testing time
data = simulate.data(c(1,5,10,100), c(.3, .1, .35, .25), 10000)