rm(list = ls())
source('data_generator.R')

run.test = function(nsamples){
  online.em.each(generate.data(nsamples, 1), generate.A())
}

generate.A = function(){
  #generate A matrix
  A = matrix(0, ncol=16, nrow=7)
  #populate A matrix
  for(i in 1:7){
    if(i <= 4){
      start = 4*(i - 1)
      for(j in 1:4){
        A[i, start + j] = 1
      }
    }
    else{
      start = i - 4
      for(j in 1:4){
        A[i, start + (j-1)*4] = 1
      }
      
    } 
  }
  A
}

online.em.each = function(data, A, verbose=F){
  x.length = 16
  
  v = data[,1]
  x = apply(data, 2, sum)
  
  lambda0 = rep(mean(x/length(data)), x.length)
  phi0 = var(v) / mean(v)
  
  theta.k = c(lambda0, phi0)
  
  s.m = s.m.fun(theta.k, A, data[1,])
  s.R = s.R.fun(theta.k, A)

  for(i in 2:nrow(data)){
    theta.k = exp(optim(log(theta.k), 
                        Q,
                        data=data[i,], A=A, S=list(m=s.m, R=s.R),
                        control=c(fnscale=-1),
                        method='L-BFGS-B',
                        lower=rep(-20, length(theta.k)),
                        upper=rep(20, length(theta.k)))
                   $par)

    lr.rate = lr.fun(i)
    s.m = s.m + lr.rate*(s.m.fun(theta.k, A, data[i,]) - s.m)
    s.R = s.R + lr.rate*(s.R.fun(theta.k, A) - s.R)
  }
  
  list(m=s.m, R=s.R)
}

s.m.fun = function(theta, A, y) {
  s.m.fun.last = function(lambda, sigma, big.multiple, A, y){
    lambda + big.multiple%*%(y - A%*%lambda)
  }
  calculate.stat(s.m.fun.last, theta, A, y)
}

s.R.fun = function(theta, A){
  s.R.fun.last = function(lambda, sigma, big.multiple, A, y){
    sigma - big.multiple%*%A%*%sigma
  }
  calculate.stat(s.R.fun.last, theta, A)
}

calculate.stat = function(stat.last.specific, theta, A, y=NULL){
  x.length = length(theta) - 1
  
  lambda = theta[1:x.length]
  phi = theta[x.length + 1]
  sigma = phi*diag(lambda)
  
  big.multiple = sigma%*%t(A)%*%qr.solve(A%*%sigma%*%t(A))
  
  stat.last.specific(lambda, sigma, big.multiple, A, y)
}

lr.fun = function(i){
  1/i
}

Q = function(theta.log, data, c, A, S){  
  theta = exp(theta.log)
  x.length = length(theta) - 1
  
  lambda = theta[1:x.length]
  phi = theta[x.length + 1]
  
  sigma = phi*diag(lambda)
  sigma.inverse = qr.solve(sigma)
  
  -(1/2)*(log(det(sigma)) + tr(sigma.inverse%*%S$R))
  - (1/2)*t(S$m - lambda)%*%sigma.inverse%*%(S$m - lambda)
}

tr = function(M){
  sum(diag(M))
}
