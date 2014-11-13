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
  
  theta.k = c(lambda0, phi0)
  theta.k1 = NULL
  while(is.null(theta.k1) || any(theta.k != theta.k1)){
    Q = function(theta, theta.k){
      sigma = theta[x.length + 1]*diag(theta[1:x.length]^c)
      sigma.k = theta.k[x.length + 1]*diag(theta.k[1:x.length]^c)
      big.multiple = t(A)%*%solve(A%*%sigma.k%*%t(A))
      r.k = sigma.k - (sigma.k%*%big.multiple%*%A)%*%sigma.k
      applied.sum = sum(apply(data, 1, function(row){
        m = theta.k[1:x.length] + sigma.k%*%big.multiple%*%(row - A%*%theta.k[1:x.length])
        t(m - theta[1:x.length])%*%sigma%*%(m - theta[1:x.length])
      }))
      -(length(data)/2)*(log(det(sigma)) + tr(qr.solve(sigma, r.k))) - (1/2)*applied.sum
    }

    theta.old = theta.k
    if(!any(is.null(theta.k1))){
      theta.k = theta.k1
    }

    theta.k1 = optim(theta.k, Q, theta.k=theta.k)$par
    print(sum((theta.k1 - theta.k)^2))
  }
  
  theta.k
}

smoothed_EM <- function(data, c, A){
  
  
}