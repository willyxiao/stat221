

online.em.each = function(data, c, A, verbose=F){
  x.length = ((ncol(data) + 1) / 2)^2
  
  v = data[,1]
  x = apply(data, 2, sum)
  
  s.m = rep(0, x.length)
  s.R = matrix(0, nrow=x.length, ncol=x.length)
  
  lambda0 = rep(mean(x/length(data)), x.length)
  phi0 = var(v) / mean(v)
  
  theta.k = c(lambda0, phi0)

  for(i in 1:ncol(data)){
    theta.k = exp(optim(log(theta.k), 
                         Q,
                         data=data[,i], c=c, A=A, S=list(m=s.m, R=s.R),
                         control=c(fnscale=-1))
                   $par)
    
    lr.rate = lr.fun(i)
    
    s.m = s.m + lr.rate*(s.m.fun(theta.k, A, data[,i]) - s.m)
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
  sigma = phi*diag(lambda^c)
  
  big.multiple = sigma%*%t(A)%*%qr.solve(A%*%sigma%*%t(A))
  
  stat.last.specific(lambda, sigma, big.multiple, A, y)
}

lr.fun = function(i){
  1/(i+1)
}

Q = function(theta.log, data, c, A, S){  
  theta = exp(theta.log)
  x.length = length(theta) - 1
  
  lambda = theta[1:x.length]
  phi = theta[x.length + 1]
  
  sigma = phi*diag(lambda^c)
  sigma.inverse = solve(sigma)
  
  -(1/2)*(log(det(sigma)) + tr(sigma.inverse%*%S$R)) 
  - (1/2)*t(S$m - lambda)%*%sigma.inverse%*%(S$m - lambda)
}