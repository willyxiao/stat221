library(matrixcalc)

simulate.data = function(nsamples){
  data = rep(0, nsamples)
  class = rep(0, nsamples)
  Us = rep(0, nsamples)
  
  for(i in 1:nsamples){
    U = runif(1, 0, 10)
    Us[i] = U
    V = rnorm(1, 0, 9)
    
    if(runif(1, 0, 1) > .5){
      class[i] = 1
      data[i] = 5*U + V
    } else {
      class[i] = 2
      data[i] = 15 + 10*U - U**2 + V
    }
  }
  list(data=data, class=class, U=Us)
}

symmetric.matrix = function(dim){
  tmp.matrix = matrix(runif(dim^2, 1, 50), nrow=dim)
  bottom.left = lower.tri(tmp.matrix, diag=FALSE)

  tmp.matrix[bottom.left] = 0
  tmp.matrix + t(tmp.matrix)
}

EPSILON = .01

online.EM = function(lr.fun, data, start.avg=50){ # learning rate...    
#   tmp = runif(1,0,1)
  tmp = .5
  
  s.1 = list(tmp,1 - tmp)
  s.2 = list(0, 0)
  s.3 = list(0, 0)
  s.4 = list(0, 0)

  A = list(0, 0)

  A[[1]] = matrix(runif(16, -10, 10), nrow=4)
  A[[2]] = A[[1]] + matrix(rnorm(16, 0, .5), nrow=4)

  for(j in 1:2){
    pos.def = t(A[[j]])%*%A[[j]]
    s.2[[j]] = pos.def[1:3, 4]
    s.3[[j]] = pos.def[1:3, 1:3]
    s.4[[j]] = pos.def[4,4]
  }
  
  w = rep(0,2)
  beta = matrix(0,nrow=2,ncol=3)
  sigma.sq = rep(0,2)
  
  for(j in 1:2){
    w[j] = s.1[[j]]
    beta[j,] = solve(s.3[[j]])%*%s.2[[j]]
    sigma.sq[j] = (s.4[[j]] - t(beta[j,])%*%s.2[[j]])/s.1[[j]]    
  }

  for(i in 1:length(data$data)){    
    r = data$data[i]
    z = c(1, data$U[i], (data$U[i]^2)/10)
    
    tmp = rep(0,2)
    for(j in 1:2){
      tmp[j] = (((w[j])/sqrt(sigma.sq[j]))
                *exp(-(1/2)
                     *((r - t(beta[j,])%*%z)^2)
                     /sigma.sq[j]))
    }

    w[1] = tmp[1] / sum(tmp)
    w[1] = min(1 - EPSILON, max(EPSILON, w[1]))
    w[2] = 1 - w[1]
    
    lr.rate = lr.fun(i)
    
    for(j in 1:2) {
      s.1[[j]] = s.1[[j]] + lr.rate*(w[j] - s.1[[j]])
      s.2[[j]] = s.2[[j]] + lr.rate*(w[j]*r*z - s.2[[j]])
      s.3[[j]] = s.3[[j]] + lr.rate*(w[j]*r*z%*%t(z) - s.3[[j]])
      s.4[[j]] = s.4[[j]] + lr.rate*(w[j]*r^2 - s.4[[j]])

      beta[j,] = solve(s.3[[j]])%*%s.2[[j]]
      sigma.sq[j] = (s.4[[j]] - t(beta[j,])%*%s.2[[j]])/s.1[[j]]
    }
  }
      
  list(w=w, beta=beta, sigma.sq=sigma.sq)
}

batch.EM = function(){
  # hannah...
}