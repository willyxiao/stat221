source('online_em_gaussian_mixture.R')

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
