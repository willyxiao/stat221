source('distro2b.R')

a.length = 10
amin = 51
amax = 100
a.tests = seq(amin, amax, length.out=a.length)

theta.run = 4

nlist = c(1e3, 1e4, 1e5, 1e6)

nreps = 10

run.test = function(){
  for(i in 1:30){
    run.task(i, 30)
  }
}

run.task = function(job.id, num.ids){
  jobs.each = num.ids / 3
  
  if(job.id <= jobs.each){
    # run.sgd
    a.id = (job.id - 1) %% jobs.each + 1
    run.job(a.id, sgd, 'sgd')
  } else if(job.id <= 2*jobs.each){
    # run.asgd
    a.id = (job.id - 1) %% jobs.each + 1
    run.job(a.id, asgd, 'asgd')
  } else {
    # run.implicit
    a.id = (job.id - 1) %% jobs.each + 1
    run.job(a.id, implicit, 'implicit')
  }
}

run.job = function(a.id, alg, alg.name){
  res = run.alg.many(nreps, a.tests[a.id], alg, nlist=nlis)#, 1e5, 1e6, 5e6))  
  save(res, file=file.name(alg.name, a.id))
}

file.name = function(alg.name, a.id){
  sprintf("out/theta%d_%s_%d.RData", theta.run, alg.name, a.id)
}
