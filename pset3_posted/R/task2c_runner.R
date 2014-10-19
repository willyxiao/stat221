source('distro2b.R')

a.length = 10
amin = 50
amax = 200
a.tests = seq(amin, amax, length.out=a.length)

theta.run = 2

nreps = 2

run.test = function(){
  for(i in 1:30){
    run.task(i, 30)
  }
}

run.task = function(job.id, num.ids){
  jobs.each = num.ids / 3
  
  if(job.id <= jobs.each){
    # run.sgd
    a.id = job.id %% jobs.each + 1
    run.job(a.id, sgd, 'sgd')
  } else if(job.id <= 2*jobs.each){
    # run.asgd
    a.id = job.id %% jobs.each + 1
    run.job(a.id, asgd, 'asgd')
  } else {
    # run.implicit
    a.id = job.id %% jobs.each + 1
    run.job(a.id, implicit, 'implicit')
  }
}

run.job = function(a.id, alg, alg.name){
  dist.list = run.alg.many(nreps, a.tests[a.id], alg, nlist=c(1e2, 1e3, 1e4, 1e5))  
  save(dist.list, file=file.name(alg.name, a.id))
}

file.name = function(alg.name, a.id){
  sprintf("out/theta%d_%s_%d.RData", theta.run, alg.name, a.id)
}
