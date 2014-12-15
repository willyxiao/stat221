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

online.EM = function(lr.fun, data, start.avg=NULL){
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
  beta = rbind(c(0, 4, 0), c(30, 0, -2))
  sigma.sq = rep(200,2)

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
    
    tmp[1] = tmp[1] / sum(tmp)
    tmp[2] = 1 - tmp[1]
    
    lr.rate = lr.fun(i)
    
    for(j in 1:2) {
      s.1[[j]] = (1 - lr.rate)*s.1[[j]] + lr.rate*tmp[j]
      s.2[[j]] = (1 - lr.rate)*s.2[[j]] + lr.rate*(tmp[j]*r*z)
      s.3[[j]] = (1 - lr.rate)*s.3[[j]] + lr.rate*(tmp[j]*z%*%t(z))
      s.4[[j]] = (1 - lr.rate)*s.4[[j]] + lr.rate*(tmp[j]*r^2)

      if(i >= INHIBITION) {
        if(is.null(start.avg) || i < start.avg){
          w[j] = s.1[[j]]
          beta[j,] = solve(s.3[[j]])%*%s.2[[j]]
          sigma.sq[j] = (s.4[[j]] - t(beta[j,])%*%s.2[[j]])/s.1[[j]]
        } else {
          avg.rate = 1/(i+1)^.5
          w[j] = (1 - avg.rate)*w[j] + avg.rate*s.1[[j]]
          beta[j,] = (1 - avg.rate)*beta[j,] + avg.rate*(solve(s.3[[j]])%*%s.2[[j]])
          sigma.sq[j] = (1 - avg.rate)*sigma.sq[j] + avg.rate*((s.4[[j]] - t(beta[j,])%*%s.2[[j]])/s.1[[j]])
        }
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

generate.betas = function(nsamples=100){
  nruns = 500
  true.beta = c(15, 10, -10)
  
  beta.1 = matrix(0, nrow=nruns, ncol=3)
  beta.2 = matrix(0, nrow=nruns, ncol=3)
  beta.3 = matrix(0, nrow=nruns, ncol=3)
  
  lr.funs = c(function(i){1/(i+1)}, function(i){1/(i+1)^.6})
  for(f.i in 1:length(lr.funs)){
    for(i in 1:500){
      tmp = matrix(NaN, nrow=1, ncol=1)
      while(is.nan(tmp[1,1])){
        tmp = online.EM(lr.funs[[f.i]], simulate.data(nsamples))$beta
      }
      res = find.best.res(tmp)
      beta.1[i,f.i] = res[1]
      beta.2[i,f.i] = res[2]
      beta.3[i,f.i] = res[3]
      print(sprintf("Completed %d for function %d", i, f.i))
    }
  }
  
  for(i in 1:500){
    tmp = matrix(NaN, nrow=1, ncol=1)
    while(is.nan(tmp[1,1])){
      tmp = online.EM(function(i){1/(i+1)^.6}, simulate.data(nsamples), 50)$beta  
    }    
    res = find.best.res(tmp)
    beta.1[i,3] = res[1]
    beta.2[i,3] = res[2]
    beta.3[i,3] = res[3]
    print(sprintf("Completed %d for OL06a", i))
  }
  
  list(beta.1 = beta.1, beta.2 = beta.2, beta.3 = beta.3)
}

plot.figure.2 = function(nsamples=100){
  x = generate.betas(nsamples)
  produce.plots(x$beta.1, x$beta.2, x$beta.3)
  x
}

produce.plots = function(beta.1, beta.2, beta.3){
  names = c("OL1", "OL06", "OL06a")
  boxplot(beta.1, names=names, ylab="beta_2(1)")
  boxplot(beta.2, names=names, ylab="beta_2(2)")
  boxplot(beta.3, names=names, ylab="beta_2(3)")
}

find.best.res = function(res){
  true.beta = c(15,10,-10)
  if(sum((res[1,] - true.beta)^2) < sum((res[2,] - true.beta)^2)){
    res[1,]
  } else{
    res[2,]
  }
}