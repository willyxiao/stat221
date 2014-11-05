source('keskici_wxiao_ps4.R')
CUTOFF = 300

log.posterior.analytic = function(N, Y){
  S = sum(Y)
  n = length(Y)
  
  numerator = sum(log(seq(1,S))) + sum(log(seq(1, n*N-S)))
  denominator = sum(log(seq(1,n*N+1))) + log(N)
  
  (numerator - denominator) + sum(log(vapply(Y, choose, 1, n=N)))
}

process.data = function(num.jobs, run.n, is.impala=TRUE){
    tt = NA
    
    for(i in 1:10){
      if(!is.impala){
        index = i + 10  
      } else{
        index = i
      }
      
      load(sprintf("keskici_wxiao_ps4_run%d_job%d.RData", run.n, index))
      r = range(chain[,1])
      tmp = table(chain[,1])
      if(length(tt) == 1 && is.na(tt)){
        tt = tmp
      } else{
        tt = c(tt, tmp)
        tt = tapply(tt, names(tt), sum)
      }
      tmp = NA
      r = NA
      gc()
    }
    
    tt
}

optim.fun = function(log.norm.const, log.analytic, log.real){
  stopifnot(length(log.analytic) == length(log.real))
  stopifnot(length(log.norm.const) == 1)
  sum(abs(log.analytic + log.norm.const - log.real))
}

run.post = function(is.impala=TRUE){
  res = process.data(NUM_JOBS, RUN_NUMBER, is.impala)
  find.probs(res)
}

find.probs = function(res){
  prop.table(res[vapply(sort(as.numeric(names(res))), toString, "hi")])
}

run.normalizing = function(res, y, cutoff=CUTOFF){
  log.analytic = vapply(as.numeric(names(res[1:cutoff])), log.posterior.analytic, .5, Y=y)
  optim(20, optim.fun, log.analytic=log.analytic, log.real=log(res[1:cutoff]), method="Brent", lower=0, upper=30)
}

analytical.cdf = function(N, Y, norm.const){
  p = 0
  for(i in max(Y):N){
    p = p + exp(log.posterior.analytic(i, Y))
  }
  p * norm.const
}

mcmc.cdf = function(N, probs){
  sum(probs[(as.numeric(names(probs)) <= N)])
}

find.mcmc.cdfs = function(run.n, is.impala=TRUE){
  cdfs = c()
  
  for(i in 1:10){
    if(!is.impala){
      index = i + 10  
    } else{
      index = i
    }
    
    load(sprintf("keskici_wxiao_ps4_run%d_job%d.RData", run.n, index))
    cdfs = c(cdfs, mcmc.cdf(100, find.probs(table(chain[,1]))))
    gc()
  }
  
  cdfs
}

compute.N.100 = function(run.n){
  impala.analytic = 1 - analytical.cdf(100, impala, exp(run.normalizing(run.post(), impala)$par))
  impala.mcmc = 1 - find.mcmc.cdfs(run.n)
  
  waterbuck.analytic = 1 - analytical.cdf(100, waterbuck, exp(run.normalizing(run.post(is.impala=FALSE), waterbuck)$par))
  waterbuck.mcmc = 1 - find.mcmc.cdfs(run.n, is.impala=FALSE)
  
  list(impala=list(analytic=impala.analytic, mcmc=impala.mcmc), waterbuck=list(analytic=waterbuck.analytic, mcmc=waterbuck.mcmc))
}

compute.N.100(RUN_NUMBER)