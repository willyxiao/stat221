source('distro2b.R')

NSIMS = 1e4
DIMS = 100

a.length = 10
amin = .05
amax = 10
a.tests = seq(amin, amax, length.out=a.length)

jobs.each = 10

theta.run = 2

m = 80

run.test = function(){
  for(i in 1:15){
    run.task(i, 15)
  }
}

run.task = function(job.id, num.ids){
  if(job.id <= jobs.each){
    # run.sgd
    a.id = job.id %% jobs.each + 1
    
    run.job(m, a.id, sgd, 'sgd')
  } else if(job.id <= 2*jobs.each){
    # run.asgd
    a.id = job.id %% jobs.each + 1
    run.job(m, a.id, asgd, 'asgd')
  
  } else {
    # run.implicit
    a.id = job.id %% jobs.each + 1
    
    if(job.id <= 3*jobs.each){
      m.group = 1
    } else{
      m.group = 2
    }

    run.job(m/2, a.id, implicit, 'implicit', m.group)
  }
}

run.job = function(m, a.id, alg, alg.name, m.group = 1){
  theta.list = as.list(rep(NA, m))
  
  for(i in 1:m){
    d = sample.data(NSIMS, DIMS)
    theta = run.method(alg, a.tests[a.id], d)
    theta.list[[i]] = theta
  }
  
  save(theta.list, file=file.name(alg.name, a.id, m.group))
}

file.name = function(alg.name, a.id, m.group){
  sprintf("out/theta%d_%s_%d_%d.RData", theta.run, alg.name, a.id, m.group)  
}

run.method = function(method, alpha, data){
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