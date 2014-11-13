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
  while(is.null(theta.k) || any(theta.k != theta.k1)){
    theta.k = theta.k1
    theta.k1 = optim(theta.k, Q, c=c, A=A, theta.k=theta.k)$par #, method="L-BFGS-B", lower=rep(1e-6, length(theta.k)))$par
    print(sum((theta.k1 - theta.k)^2))
  }
  
  theta.k
}

Q = function(theta, theta.k, c, A){
  sigma = theta[x.length + 1]*diag(theta[1:x.length]^c)
  sigma.k = theta.k[x.length + 1]*diag(theta.k[1:x.length]^c)
  sigma.inverse = qr.solve(sigma)
  big.multiple = t(A)%*%qr.solve(A%*%sigma.k%*%t(A))
  r.k = sigma.k - (sigma.k%*%big.multiple%*%A)%*%sigma.k
  applied.sum = sum(apply(data, 1, function(row){
    m = theta.k[1:x.length] + sigma.k%*%big.multiple%*%(row - A%*%theta.k[1:x.length])
    t(m - theta[1:x.length])%*%sigma.inverse%*%(m - theta[1:x.length])
  }))
  res = -(length(data)/2)*(log(det(sigma)) + tr(sigma.inverse%*%r.k)) - (1/2)*applied.sum
#   if(res == Inf){
#     1e30
#   } else if (res == -Inf){
#     -1e30
#   } else{
#     res
#   }
  res
}

smoothed_EM <- function(data, c, A){
  
  
}