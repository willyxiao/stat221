source('keskici_wxiao_ps4.R')

log.posterior.analytic = function(N, Y){
  S = sum(Y)
  n = length(Y)
  
  numerator = sum(log(seq(1,S))) + sum(log(seq(1, n*N-S)))
  denominator = sum(log(seq(1,n*N+1))) + log(N)
  
  (numerator - denominator) + sum(log(vapply(Y, choose, 1, n=N)))
}

optim.fun = function(log.norm.const, log.analytic, log.real){
  stopifnot(length(log.analytic) == length(log.real))
  stopifnot(length(log.norm.const) == 1)
  sum((log.analytic + log.norm.const - log.real)^2)
}

process.data = function(num.jobs, run.n){
    tt = NA
    
    for(i in 11:20){
      load(sprintf("keskici_wxiao_ps4_run%d_job%d.RData", run.n, i))
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

CUTOFF = 400

run.post = function(){
  res = process.data(NUM_JOBS, RUN_NUMBER)
  res = prop.table(res[vapply(sort(as.numeric(names(res))), toString, "hi")])
  barplot(res)
  res
}

run.normalizing = function(res, y){
  tmp = vapply(as.numeric(names(res[1:CUTOFF])), log.posterior.analytic, .5, Y=y)
  barplot(tmp)
  optim(20, optim.fun, log.analytic=tmp, log.real=log(res[1:CUTOFF]), method="Brent", lower=0, upper=20)  
}


