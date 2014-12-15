rm(list=ls())

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

INHIBITION = 20

online.EM = function(lr.fun, data, start.avg=50){
  tmp = .5
  
  s.1 = list(tmp,1 - tmp)
  s.2 = list(0, 0)
  s.3 = list(0, 0)
  s.4 = list(0, 0)

  A = list(0, 0)

  A[[1]] = matrix(runif(16, 0, 3), nrow=4)
  A[[2]] = matrix(runif(16, 0, 3), nrow=4)

  for(j in 1:2){
    pos.def = t(A[[j]])%*%A[[j]]
    s.2[[j]] = pos.def[1:3, 4]
    s.3[[j]] = pos.def[1:3, 1:3]
    s.4[[j]] = pos.def[4,4]
  }

  w = rep(.5,2)
  beta = rbind(c(0, 4, 0), c(10, 10, 10))
  sigma.sq = rep(200,2)
  
#   for(j in 1:2){
#     w[j] = s.1[[j]]
#     beta[j,] = solve(s.3[[j]])%*%s.2[[j]]
#     sigma.sq[j] = (s.4[[j]] - t(beta[j,])%*%s.2[[j]])/s.1[[j]]
#   }

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
    w[2] = 1 - w[1]
    
    lr.rate = lr.fun(i)
    
    for(j in 1:2) {
      s.1[[j]] = (1 - lr.rate)*s.1[[j]] + lr.rate*w[j]
      s.2[[j]] = (1 - lr.rate)*s.2[[j]] + lr.rate*(w[j]*r*z)
      s.3[[j]] = (1 - lr.rate)*s.3[[j]] + lr.rate*(w[j]*z%*%t(z))
      s.4[[j]] = (1 - lr.rate)*s.4[[j]] + lr.rate*(w[j]*r^2)

      if(i >= INHIBITION) {
        w[j] = s.1[[j]]
        beta[j,] = solve(s.3[[j]])%*%s.2[[j]]
        sigma.sq[j] = (s.4[[j]] - t(beta[j,])%*%s.2[[j]])/s.1[[j]]
      }
    }
  }

  list(w=w, beta=beta, sigma.sq=sigma.sq)
}

plot.figure.1 = function(){
  data = simulate.data(500)
  plot(data$U[data$class == 1], data$data[data$class ==1], xlab="U", ylab="Z")
  points(data$U[data$class == 2], data$data[data$class ==2], pch=4, col='red')
}
