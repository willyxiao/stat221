source('keskici_wxiao_ps4.R')

log.posterior.analytic = function(N, Y){
  S = sum(Y)
  n = length(Y)
  
  numerator = sum(log(seq(1,S))) + sum(log(seq(1, n*N-S)))
  denominator = sum(log(seq(1,n*N+1))) + N
  
  (numerator - denominator) + sum(log(vapply(Y, choose, 1, n=N)))
}

optim.fun = function(norm.const, log.analytic, log.real){
  sum((log.analytic + log(norm.const) - log.real)^2)
}

process.data = function(num.jobs, run.n){
    tt = NA
    
    for(i in 2:10){
      load(sprintf("keskici_wxiao_ps4_run%d_job%d.RData", run.n, i))
      r = range(chain[,1])
      tmp = table(chain[,1])
      if(is.na(tt)){
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

res = process.data(NUM_JOBS, RUN_NUMBER)
res = prop.table(res[vapply(sort(as.numeric(names(res))), toString, "hi")])
barplot(res)
tmp = vapply(as.numeric(names(res)), log.posterior.analytic, .5, Y=impala)
optim(10, optim.fun, log.analytic=tmp, log.real=log(res), method="Brent", lower=0, upper=1e5)
