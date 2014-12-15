rm(list=ls())
library(matrixcalc)

simulate.data = function(nsamples){
  data = rep(0, nsamples)
  class = rep(0, nsamples)
  Us = rep(0, nsamples)
  
  for(i in 1:nsamples){
    U = runif(1, 0, 10)
    V = rnorm(1, 0, 9)

    Us[i] = U
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

online.EM = function(lr.fun, data, start.avg=NULL, inhibition=20){  
  s = initialize.s()  
  s.1 = s$s.1
  s.2 = s$s.2
  s.3 = s$s.3
  s.4 = s$s.4

  theta = initialize.theta()
  w = theta$w
  beta = theta$beta
  sigma.sq = theta$sigma.sq

  # to store betas later
  betas = list(matrix(0, nrow=length(data$data), ncol=3), 
               matrix(0, nrow=length(data$data), ncol=3))
  
  for(i in 1:length(data$data)){  
    r = data$data[i]
    z = c(1, data$U[i], (data$U[i]^2)/10)
    
    w.hat = rep(0,2)
    for(j in 1:2){
      w.hat[j] = w.hat.fun(w[j], beta[j,], sigma.sq[j], r, z)
    }
    w.hat[1] = w.hat[1] / sum(w.hat)
    w.hat[2] = 1 - w.hat[1]
    
    lr.rate = lr.fun(i)
    for(j in 1:2) {
      s.hat = s.fun(w.hat[j], r, z)
      new.s = stochastic.update(lr.rate, 
                                s.1[[j]], s.hat$s.1,
                                s.2[[j]], s.hat$s.2,
                                s.3[[j]], s.hat$s.3,
                                s.4[[j]], s.hat$s.4)
      s.1[[j]] = new.s[[1]]
      s.2[[j]] = new.s[[2]]
      s.3[[j]] = new.s[[3]]
      s.4[[j]] = new.s[[4]]

      if(i >= inhibition) {
        theta = theta.fun(s.1[[j]], s.2[[j]], s.3[[j]], s.4[[j]])
        avg.rate = ifelse(is.null(start.avg) || i < start.avg, 1, 1/(i+1)^.5)
        new.theta = stochastic.update(avg.rate,
                                      w[j], theta$w,
                                      beta[j,], theta$beta,
                                      sigma.sq[j], theta$sigma.sq)
        
        w[j] = new.theta[[1]]
        beta[j,] = new.theta[[2]]
        sigma.sq[j] = new.theta[[3]]
      }
      
      betas[[j]][i,] = beta[j,]
    }
  }

  list(w=w, beta=beta, sigma.sq=sigma.sq, betas=betas)
}

# online-EM helpers
stochastic.update = function(lr.rate, ...){
  elements = list(...)
  results = list()
  stopifnot(length(elements) %% 2 == 0)
  odd.indices = seq(1, length(elements), 2)  
  for(i in odd.indices){
    results[[floor(i/2)+1]] = (1 - lr.rate)*elements[[i]] + lr.rate*elements[[i+1]]
  }
  results
}

initialize.theta = function(){
  list(w=rep(.5,2), beta=rbind(c(0, 4, 0), c(30, 0, -2)), sigma.sq=rep(200,2))
}

initialize.s = function(){
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
  
  list(s.1=s.1,s.2=s.2,s.3=s.3,s.4=s.4)
}

w.hat.fun = function(w,beta,sigma.sq,r,z){
  (w/sqrt(sigma.sq)
   *exp(-(1/2)
        *(r - t(beta)%*%z)^2
        /sigma.sq))
}

s.fun = function(w.hat, r, z){
  list(s.1=w.hat, s.2=w.hat*r*z, s.3=w.hat*z%*%t(z), s.4=w.hat*r^2)  
}

theta.fun = function(s.1, s.2, s.3, s.4){
  w = s.1
  beta = solve(s.3)%*%s.2
  sigma.sq = (s.4 - t(beta)%*%s.2)/s.1
  
  list(w=w, beta=beta, sigma.sq=sigma.sq)
}

# PLOTTING

plot.figure.1 = function(){
  data = simulate.data(500)
  plot(data$U[data$class == 1], data$data[data$class ==1], xlab="U", ylab="Z")
  points(data$U[data$class == 2], data$data[data$class ==2], pch=4, col='red')
}

plot.figure.2 = function(nsamples=100){
  x = generate.betas(nsamples)
  produce.plots.fig.2(x$beta.1, x$beta.2, x$beta.3)
  x
}

plot.figure.4 = function(nsamples=5000){
  beta.1 = matrix(0, nrow=nsamples, ncol=3)
  beta.2 = matrix(0, nrow=nsamples, ncol=3)
  beta.3 = matrix(0, nrow=nsamples, ncol=3)
  
  lr.funs = c(function(i){1/(i+1)}, function(i){1/(i+1)^.6})
  
  for(f.i in 1:length(lr.funs)){
    res = online.EM(lr.funs[[f.i]], simulate.data(nsamples))$betas
    index = find.best.res.i(res[[1]][nsamples,], res[[2]][nsamples,])    
    beta.1[,f.i] = res[[index]][,1]
    beta.2[,f.i] = res[[index]][,2]
    beta.3[,f.i] = res[[index]][,3]
  }
  
  res = online.EM(function(i){1/(i+1)^.6}, simulate.data(nsamples), 1000)$betas
  index = find.best.res.i(res[[1]][nsamples,], res[[2]][nsamples,])
  beta.1[,3] = res[[index]][,1]
  beta.2[,3] = res[[index]][,2]
  beta.3[,3] = res[[index]][,3]
  
  betas = list(beta.1 = beta.1, beta.2 = beta.2, beta.3 = beta.3)
  
  plot(1:nsamples, beta.1[,1], type="l", col="red", lty=2, ylim=c(0,30))
  lines(1:nsamples, beta.1[,2], lty=3)
  lines(1:nsamples, beta.1[,3], col="blue")
  
  plot(1:nsamples, beta.2[,1], type="l", col="red", lty=2, ylim=c(0,20))
  lines(1:nsamples, beta.2[,2], lty=3)
  lines(1:nsamples, beta.2[,3], col="blue")
  
  plot(1:nsamples, beta.3[,1], type="l", col="red", lty=2, ylim=c(-20,0))
  lines(1:nsamples, beta.3[,2], lty=3)
  lines(1:nsamples, beta.3[,3], col="blue")
  
  betas
}

# PLOTTING HELPERS

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

produce.plots.fig.2 = function(beta.1, beta.2, beta.3){
  names = c("OL1", "OL06", "OL06a")
  boxplot(beta.1, names=names, ylab="beta_2(1)")
  boxplot(beta.2, names=names, ylab="beta_2(2)")
  boxplot(beta.3, names=names, ylab="beta_2(3)")
}

find.best.res.i = function(a, b){
  true.beta = c(15,10,-10)
  ifelse(sum((a - true.beta)^2) < sum((b - true.beta)^2), 1, 2)
}

find.best.res = function(res){
  res[find.best.res.i(res[1,], res[2,]),]
}