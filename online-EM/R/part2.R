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

online.EM = function(lr.fun, data, start.avg=50){ # learning rate...
  w = c(.5,.5)
  w.new = w
  
  r = data$data[1]
  z = c(1,data$U[1], (data$U[1]^2)/10)
  
  s.1 = list(.2,.8)
  s.2 = list(.5*r*z, .5*r*z)
  s.3 = list(matrix(runif(9,1,10), nrow=3), matrix(runif(9,1,10), nrow=3))
  s.4 = list(.5*r^2, .5*r^2)
  
  for(i in 1:20){
    r = data$data[i]
    z = c(1, data$U[i], (data$U[i]^2)/10)
    
    for(j in 1:2){
      lr.rate = lr.fun(i)
      s.1[[j]] = s.1[[j]] + lr.rate*(w[j] - s.1[[j]])
      s.2[[j]] = s.2[[j]] + lr.rate*(w[j]*r*z - s.2[[j]])
      s.3[[j]] = s.3[[j]] + lr.rate*(w[j]*r*z%*%t(z) - s.3[[j]])
      s.4[[j]] = s.4[[j]] + lr.rate*(w[j]*r^2 - s.4[[j]])
    }
  }

  beta = matrix(0,nrow=2,ncol=3)
  sigma.sq = rep(0,2)
  
  for(j in 1:2){
    w[j] = s.1[[j]]
    beta[j,] = solve(s.3[[j]])%*%s.2[[j]]
    sigma.sq[j] = (s.4[[j]] - t(beta[j,])%*%s.2[[j]])/w[j]    
  }
  
  beta.new = beta
  sigma.sq.new = sigma.sq
    
  for(i in 20:length(data$data)){
    r = data$data[i]
    z = c(1, data$U[i], (data$U[i]^2)/10)
    
    w.1 = (((w[1])/sigma.sq[1])
            *exp(-(1/2)
                 *(r - t(beta[1,])%*%z)
                 /sigma.sq[1]))
    w.2 = (((w[2])/sigma.sq[2])
           *exp(-(1/2)
                *(r - t(beta[2,])%*%z)
                /sigma.sq[2]))
    
    w.new[1] = w.1 / (w.1 + w.2)
    
    if(w.new[1] == 1){
      w.new[1] = .99
    } else if(w.new[1] == 0){
      w.new[1] = .01
    }
    
    w.new[2] = 1 - w.new[1]
    
    for(j in 1:2) { # only 2 classes
      lr.rate = lr.fun(i)
      s.1[[j]] = s.1[[j]] + lr.rate*(w.new[j] - s.1[[j]])
      s.2[[j]] = s.2[[j]] + lr.rate*(w.new[j]*r*z - s.2[[j]])
      s.3[[j]] = s.3[[j]] + lr.rate*(w.new[j]*r*z%*%t(z) - s.3[[j]])
      s.4[[j]] = s.4[[j]] + lr.rate*(w.new[j]*r^2 - s.4[[j]])

      beta.new[j,] = solve(s.3[[j]])%*%s.2[[j]]
      sigma.sq.new[j] = (s.4[[j]] - t(beta.new[j,])%*%s.2[[j]])/w.new[j]
    }
    
    if(i >= start.avg && FALSE){      
#       lr.rate = lr.fun(i)
#       w = w + lr.rate*(w.new - w)
#       beta = beta + lr.rate*(beta.new - beta)
#       sigma.sq = sigma.sq + lr.rate*(sigma.sq.new - sigma.sq)
    } else {
      w = w.new
      beta = beta.new
      sigma.sq = sigma.sq.new
    }
  }  
  
  list(w=w, beta=beta, sigma.sq=sigma.sq)
}

batch.EM = function(){
  # hannah...
}



# plot(tmp$U[temp$class == 1], tmp$data[tmp$class == 1])
# points(tmp$U[tmp$class ==2], tmp$data[tmp$class == 2])