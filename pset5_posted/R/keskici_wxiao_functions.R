rm(list=ls())

# order is f,s,l,c. 
# y = src f, src s, ... dst f, dst s, ... 
# x = f->f, f->s, ... s->f, ... c->l, c->c

locally_iid_EM <- function(data, c, A, w=11){
  # do some shit
}

x.length = 16
tr = function(M){
  sum(diag(M))
}

locally_iid_EM.each = function(data, c, A){
  v = data[,1]
  x = apply(data, 2, sum)
  
  lambda0 = rep(mean(x/length(data)), x.length)
  phi0 = var(v) / mean(v)
  
  theta.k = NULL
  theta.k1 = c(lambda0, phi0)
  
  log.lik.k = NULL
  log.lik.k1 = 0
  
  while(is.null(theta.k) || log.lik.k1 - log.lik > 5){
    theta.k = theta.k1
    theta.k1 = optim(theta.k, Q, c=c, A=A, theta.k=theta.k, control=c(fnscale=-1), method="L-BFGS-B", lower=rep(1e-6, length(theta.k)))$par #, )$par
    
    log.lik.k = log.lik(data, theta.k, c, A)
    log.lik.k1 = log.lik(data, theta.k1, c, A)    
    print(log.lik.k1 - log.kik.k)
  }
  
  theta.k
}

log.lik = function(data, theta, c, A){
  T = length(data)
  
  lambda = theta[1:x.length]
  phi = theta[x.length + 1]
  
  sigma = phi*diag(lambda^c)
  
  applied.sum = sum(apply(data, 1, function(y){
    t(y - A%*%lambda)%*%qr.solve(A%*%sigma%*%t(A))%*%(y - A%*%lambda)
  }))
  
  -(T/2)*log(det(A%*%sigma%*%t(A))) - (1/2)*applied.sum
}

Q = function(theta, theta.k, c, A){  
  lambda = theta[1:x.length]
  phi = theta[x.length + 1]
  lambda.k = theta.k[1:x.length]
  phi.k = theta.k[x.length + 1]
  
  sigma = phi*diag(lambda^c)
  sigma.k = phi.k*diag(lambda.k^c)
  sigma.inverse = qr.solve(sigma)
  big.multiple = sigma.k%*%t(A)%*%qr.solve(A%*%sigma.k%*%t(A))
  r.k = sigma.k - big.multiple%*%A%*%sigma.k
  applied.sum = sum(apply(data, 1, function(row){
    m = lambda.k + big.multiple%*%(row - A%*%lambda.k)
    t(m - lambda)%*%sigma.inverse%*%(m - lambda)
  }))
  res = -(length(data)/2)*(log(det(sigma)) + tr(sigma.inverse%*%r.k)) - (1/2)*applied.sum
   if(res == Inf){
     res = 1e300
   } else if (res == -Inf){
     res = -1e300
   } else{
     res
   }
  res
}

smoothed_EM <- function(data, c, A){
  
  
}