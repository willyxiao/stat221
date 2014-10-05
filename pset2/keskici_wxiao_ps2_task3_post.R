#! /usr/bin/env Rscript

source('keskici_wxiao_ps2_functions.R')

TASK.NUM = 3

J           <- 1000 # length of theta and w vector
theta.draws <- 30    # theta.nsims * N is the # of total simulations we'll run
Y.draws     <- 12

w           <- rep(1, J) # weights are all set to 1
mu          <- c(1.6, 2.5, 5.2, 4.9)
sigma       <- c(0.7, 1.3, 1.3, 1.6)

stopifnot(length(mu) == length(sigma))

is.covered <- function(lower, higher, x){
  as.numeric(lower <= x && x <= higher)
}

for(pair in 1:length(mu)){
  for(theta.draw in 1:theta.draws){
    print(sprintf("Writing par%d, theta%d", pair, theta.draw))
   
    is.covered.95 = NULL
    is.covered.68 = NULL
    
    real.log.theta = NULL
    
    for(Y.draw in 1:Y.draws){
      load(getOutFileName(pair,theta.draw,Y.draw))
      
      CI.95.lower = apply(res$logTheta, 1, quantile, probs=.025, names=FALSE)
      CI.68.lower = apply(res$logTheta, 1, quantile, probs=.16, names=FALSE)
      
      CI.95.higher = apply(res$logTheta, 1, quantile, probs=.975, names=FALSE)
      CI.68.higher = apply(res$logTheta, 1, quantile, probs=.84, names=FALSE)
      
      is.covered.95.i = mapply(is.covered, CI.95.lower, CI.95.higher, res$real.log.theta)
      is.covered.68.i = mapply(is.covered, CI.68.lower, CI.68.higher, res$real.log.theta)
            
      if (is.null(is.covered.95)) {
        stopifnot(is.null(is.covered.68) && is.null(real.log.theta))
        
        real.log.theta = res$real.log.theta
        is.covered.95 = is.covered.95.i
        is.covered.68 = is.covered.68.i
      } else {
        is.covered.95 = cbind(is.covered.95, is.covered.95.i)
        is.covered.68 = cbind(is.covered.68, is.covered.68.i)
      }
    }
    
    cover.95 = apply(is.covered.95, 1, mean)
    cover.68 = apply(is.covered.68, 1, mean)
    
    to.write = cbind(real.log.theta, cover.95, cover.68)
    write.table(to.write, 
                file=sprintf("keskici_wxiao_ps2_task%d_par%d_theta.dat", TASK.NUM, pair), 
                append=TRUE,
                row.names=FALSE,
                col.names=FALSE)
  }
  plot.coverage(sprintf("keskici_wxiao_ps2_task%d_par%d_theta.dat", TASK.NUM, pair), TASK.NUM, pair)
  
}
