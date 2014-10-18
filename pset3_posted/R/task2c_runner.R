source('task2b.R')

methods = list(sgd, implicit, asgd)

a.length = 10
amin = .05
amax = 10
a.tests = seq(amin, amax, length.out=a.length)

m = 120

run.task = function(job.id, num.ids){
  if(job.id <= 10){
    # run.sgd
    a.id = job.id %% 10
    
    run.sgd(m, a)

  } else if(job.id <= 20){
    # run.asgd
    a.id = job.id %% 10
    
  } else {
    # run.implicit
    
  }
}

run.method = function(method, alpha){
  # check.data(data)
  n = nrow(data$X)
  p = ncol(data$X)

  theta.list = matrix(0, nrow=p, ncol=1)
  
  for(i in 1:n) {
    xi = data$X[i, ]
    ai = alpha / (alpha + i)
    yi = data$Y[i]
    theta.old = theta.list[, i]

    theta.new = method(ai, xi, yi, i, theta.old)

    theta.list = cbind(theta.list, theta.new)
  }
  
  theta.list
}

sgd = function(ai, xi, yi, i, theta.old){
  lpred = sum(theta.old * xi)
  (theta.old - ai * lpred * xi) + (ai * yi * xi)
}

asgd = function(ai, xi, yi, i, theta.old){
  (1 - 1/i) * theta.old + (1/i) * sgd(ai, xi, yi, i, theta.old)
}

implicit = function(ai, xi, yi, i, theta.old){
  solve(diag(length(xi)) + ai*xi%*%t(xi))%*%(theta.old + ai*yi*xi)
}