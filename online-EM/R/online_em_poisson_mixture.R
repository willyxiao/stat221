#simulate data
simulate.data = function(lambda, w, n){
  m = length(lambda)
  #helper function to determine which lambda
  #to use
  data = rep(NA, n)
  for(i in 1:n){
    data[i] = rpois(1, sum(rmultinom(1, 1, w) * lambda))
  }
  data
}

set.init.lambdas = function(data, m){
  quantiles = seq(1/m, 1, 1/m)
  quantiles = quantiles - 1/(2*m)
  quantile(data, quantiles)
}

w.bar = function(y, w, lambda){
  num.mixtures = length(lambda)
  nums = rep(NaN, num.mixtures)
  for(i in 1:num.mixtures){
    nums[i] = w[i] * dpois(y, lambda[i], log = FALSE)
    if(is.nan(nums[i])){
      browser()
    }
    if(nums[i] == 0){
      nums[i] = .0001
    }
  }
  
  den = sum(nums)
  nums/den
  
}

poisson.em = function(data, m){
  #init values
  w = rep(1/m, m)
  lambda = set.init.lambdas(data, m)
  s.hat = cbind(w, lambda)
  for(i in 1: length(data)){
    w.bar.iter = w.bar(data[i], w, lambda)
    for(j in 1:m){
      s.hat[j,] = s.hat[j,] + (1/(i + 1)) *(c(w.bar.iter[j], w.bar.iter[j]*data[i]) - s.hat[j,])
    }
    w = s.hat[,1]
    lambda = s.hat[,2] / s.hat[,1]
  }
  c(lambda, w)
}

#testing time
#trivial example as sanity check
x = c(1)
data = simulate.data(c(100), x, 10000)
test1 = poisson.em(data, 1)

x = c(.1, .25, .5)
data = simulate.data(c(25,40,10,100), c(x, 1-sum(x)), 100000)
test2 = poisson.em(data, 4)

x = c(.5, .5)
data = simulate.data(c(10, 2), x, 10000)
test3 = poisson.em(data, 2)



