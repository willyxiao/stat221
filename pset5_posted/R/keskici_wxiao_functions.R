rm(list=ls())

# order is f,s,l,c. 
# y = src f, src s, ... dst f, dst s, ... 
# x = f->f, f->s, ... s->f, ... c->l, c->c

locally_iid_EM <- function(data, c, A, w=11){
  # do some shit
}

locally_iid_EM.each = function(data, c, A){
  v = data[,1]
  x = apply(data, 2, sum)
  
  lambda0 = mean(x/length(data))
  phi0 = var(v) / mean(v)
  
  theta.k = c(lambda0, phi0)
  theta.k1 = NULL
  
  while(any(theta.k != theta.k1)){
    Q = function(theta){
      sigma = theta[2]*diag(theta[1]^c)
      sigma.k = theta.k[2]*diag(theta.k[1]^c)
      A.inverse = solve(A)
      big.multiple = A.inverse%*%solve(A%*%sigma.k%*%A.inverse)
      r.k = sigma.k - sigma.k%*%big.multiple%*%A%*%sigma.k
      s = sum(apply(data, 1, function(row){
        m = theta.k[1] + sigma.k%*%big.multiple%*%(row - A%*%theta.k[1])
        solve(m - theta[1])%*%sigma%*%(m - theta[1])
      }))
      -(length(data)/2)*(log(sigma) + tr(solve(sigma)%*%r.k)) - (1/2)*s
    }
    theta.old = theta.k
    theta.k = theta.k1
    theta.k1 = optim(theta.old, Q)$par
  }
  
  theta.k
}

smoothed_EM <- function(data, c, A){
  
  
}