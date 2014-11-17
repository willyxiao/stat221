rm(list=ls())

library(mvtnorm)

# order is f,s,l,c. 
# y = src f, src s, ... dst f, dst s, ... 
# x = f->f, f->s, ... s->f, ... c->l, c->c
tr = function(M){
  sum(diag(M))
}

locally_iid_EM = function(data, c, A, w=11){
  x.length = ((ncol(data) + 1) / 2)^2
  
  h = floor(w/2)
  estimates = c()
  for(i in (h+1):(nrow(data))){
    subset = data[max(1, i-h):min(nrow(data), i+h),]
    results = locally_iid_EM.each(subset, c, A)
    print(c(i, mean(results[1:x.length])))
    estimates = rbind(estimates, results)
  }
  estimates
}

smoothed_EM <- function(data, c, A, w=11){
  h = floor(w/2)

  V = 5*diag(lambda0)
  
  eta.t.m1 = log(c(lambda0, phi0))
  sigma.t.m1 = phi0*diag(lambda0^2)
  
  for(t in 1:nrow(data)){
      subset = data[max(1, t-h):min(nrow(data), i+h)]
      sigma.t = sigma.t.m1 + V

#       g = function(eta.t){
#         log(dmvnorm(eta.t, eta.k.m1, sigma.t)) + log.lik(subset, exp(eta.t), c, A)
#       }
# 
#       eta.t = optim(eta.t.m1, g, control=c(fnscale=-1))$par
      eta.t = smoothed_EM.each(subset, c, A)
      eta.t.m1 = eta.t
      
      lambda.k = exp(eta.t)[1:(length(eta.t)-1)]
      phi = eta.t[length(eta.t)]
      
      m.k = 1/(nrow(data))*apply(subset, 1, function(y){
        lambda.k + sigma.t%*%t(A)%*%qr.solve(A%*%sigma.t%*%t(A))%*%(y - A%*%lambda.k)
      })
      b.k = apply(m.k, 2, sum)

      top.left = diag(phi*c^2*lambda.k^(c-1) + 2*(2-c)*lambda.k - 2*(1-c)*b.k)
      left = rbind(top.left, c*lambda.k^c)
      right.col = (2-c)*lambda.k^(1-c)-(1-c)*lambda.k^(-c)*b.k
      second.der = cbind(left, right.col)
      
      sigma.t.m1 = qr.solve(sigma.t) + second.dir

  }
  
  eta.t
}

smoothed_EM.each = function(sigma.t, data, c, A, verbose=F){
  x.length = ((ncol(data) + 1) / 2)^2
  
  v = data[,1]
  x = apply(data, 2, sum)
  
  lambda0 = rep(mean(x/length(data)), x.length)
  phi0 = var(v) / mean(v)
  
  eta.k = NULL
  eta.k1 = c(lambda0, phi0)
  
  post.k = NULL
  post.k.k1 = 0
  
  while(is.null(eta.k) || (post.k.k1 - post.k) > 2){
    
    eta.k = eta.k1
    eta.k1 = optim(eta.k, g, sigma.t=sigma.t, data=data, c=c, A=A, eta.k.m1=eta.k, control=c(fnscale=-1))$par
    
    post.k = post.prob.k.k1
    post.k.k1 = post.prob(eta.k1, eta.k, sigma.t1, data, c, A)
  }

  eta.k
}

g = function(eta, eta.k.m1, sigma.t, data, c, A){
  log(dmvnorm(eta, eta.k.m1, sigma.t)) + Q(eta, exp(eta.k.m1), data, c, A)
}

post.prob = function(eta.t, eta.t.m1, sigma.t, data, c, A){
  log(dmvnorm(eta.t, eta.t.m1, sigma.t)) + log.lik(data, exp(eta.t), c, A)
}

locally_iid_EM.each = function(data, c, A, verbose=F){
  x.length = ((ncol(data) + 1) / 2)^2
  
  v = data[,1]
  x = apply(data, 2, sum)
  
  lambda0 = rep(mean(x/length(data)), x.length)
  phi0 = var(v) / mean(v)
  
  theta.k = NULL
  theta.k1 = c(lambda0, phi0)
  
  log.lik.k = NULL
  log.lik.k1 = 0
  
  while(is.null(theta.k) || (log.lik.k1 - log.lik.k) > 2){
    theta.k = theta.k1
#     theta.k1 = optim(theta.k, Q, data=data, c=c, A=A, theta.k=theta.k, control=c(fnscale=-1), method="L-BFGS-B", lower=rep(1e-6, length(theta.k)))$par #, )$par
    theta.k1 = exp(optim(log(theta.k), Q, data=data, c=c, A=A, theta.k=theta.k, control=c(fnscale=-1))$par) #, )$par
    
    log.lik.k = log.lik(data, theta.k, c, A)
    log.lik.k1 = log.lik(data, theta.k1, c, A)    
    
    if(verbose){
      print(c(log.lik.k, log.lik.k1))
    }
  }
  
  theta.k
}

log.lik = function(data, theta, c, A){
  x.length = length(theta) - 1
  
  T = length(data)
  
  lambda = theta[1:x.length]
  phi = theta[x.length + 1]
  
  sigma = phi*diag(lambda^c)
  
  applied.sum = sum(apply(data, 1, function(y){
    t(y - A%*%lambda)%*%qr.solve(A%*%sigma%*%t(A))%*%(y - A%*%lambda)
  }))
  
  -(T/2)*log(det(A%*%sigma%*%t(A))) - (1/2)*applied.sum
}

Q = function(theta, theta.k, data, c, A){  
  theta = exp(theta)
  x.length = length(theta) - 1
  
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
#    if(res == Inf){
#      res = 1e300
#    } else if (res == -Inf){
#      res = -1e300
#    } else{
#      res
#    }
  res
}

plot.fig5 = function(res, dim, names, indices, filename, ymax){
  pdf(filename)
  par(mfrow=c(dim,dim))
  for(i in 1:length(indices)){
    plot(res[,indices[i]], type='l', ylim = c(0, ymax), main = names[i], 
         xlab = "hour of day", ylab ="bytes/sec")
  }
  dev.off()
}

