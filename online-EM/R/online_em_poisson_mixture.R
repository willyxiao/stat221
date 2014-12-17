#simulate data
simulate.data = function(lambda, w, n){
  m = length(lambda)
  data = rep(NA, n)
  for(i in 1:n){
    data[i] = rpois(1, sum(rmultinom(1, 1, w) * lambda))
  }
  data
}

#plucks quantiles off data to use
#as initial lambdas
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
test1.probs = c(1)
test1.data = simulate.data(c(100), test1.probs, 100000)
test1 = poisson.em(test1.data, 1)

test2.probs = c(.5, .3, .2)
test2.data = simulate.data(c(2, 15, 30), test2.probs, 100000)
test2 = poisson.em(test2.data, 3)

test3.probs = c(.75, .25)
test3.data = simulate.data(c(10, 2), test3.probs, 100000)
test3 = poisson.em(test3.data, 2)

test.set.init.lambdas = function(){
  stopifnot(set.init.lambdas(seq(0,100),4) == c(12.5, 37.5, 62.5, 87.5))
  stopifnot(set.init.lambdas(seq(0,100),1) == 50)
  stopifnot(set.init.lambdas(seq(0,100),2) == c(25, 75))  
}
test.set.init.lambdas()