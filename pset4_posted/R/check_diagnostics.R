load('keskici_wxiao_ps4_run3_job12.Rdata')


test.chain <- function(mh.draws, show.diagnostics=F) {  
  ##  Do the diagnostics here.
  if (show.diagnostics) {
    library(coda)
    mcmc.chain <- mcmc(mh.draws)
    fns <- c("summary", "plot", "autocorr.plot", "rejectionRate")
    for (fn in fns) {
      readline(sprintf("Press [ENTER] for %s", fn))
      do.call(fn, args=list(mcmc.chain))
    }
  }
  print("Geweke diagnostic")
  print(geweke.diag(mh.draws))
  
  print("Raftery diagnostic")
  print(raftery.diag(mh.draws, r=0.005))
  
  print("Heidelberg diagnostic")
  print(heidel.diag(mh.draws))
}
test.chain(chain, T)