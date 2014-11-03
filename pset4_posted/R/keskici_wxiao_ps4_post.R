source('keskici_wxiao_ps4.R')

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

CUTOFF = 700

optim.fun = function(log.norm.const, log.analytic, log.real){
  stopifnot(length(log.analytic) == length(log.real))
  stopifnot(length(log.norm.const) == 1)
  sum((log.analytic + log.norm.const - log.real)^2)
}

run.post = function(is.impala=TRUE){
  res = process.data(NUM_JOBS, RUN_NUMBER, is.impala)
  res = prop.table(res[vapply(sort(as.numeric(names(res))), toString, "hi")])
  barplot(res)
  res
}

run.normalizing = function(res, y, cutoff=CUTOFF){
  log.analytic = vapply(as.numeric(names(res[1:cutoff])), log.posterior.analytic, .5, Y=y)
  optim(20, optim.fun, log.analytic=log.analytic, log.real=log(res[1:cutoff]), method="Brent", lower=0, upper=30)  
}


