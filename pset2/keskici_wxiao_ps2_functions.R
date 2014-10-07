simYgivenTheta <- function(theta, w, N) {
  stopifnot(length(theta) == length(w))  
  J = length(theta)
  matrix(rpois(N*J, w*theta), nrow=J,ncol=N)
}

simThetagivenMuSigma <- function(J, mu, sigma) {
  rnorm(J, mu, sigma)
}

# same as getFileName but in the out file and called "out"
getOutFileName <- function(theta.num, pair.num, task.num, w = FALSE){
    file.name.str = ""
    if(w){
      file.name.str <- "out/wout%d_%d_%d.RData"
    } else{
      file.name.str <- "out/out%d_%d_%d.RData"      
    }
    sprintf(file.name.str, task.num, theta.num, pair.num)
}


is.covered <- function(lower, higher, x){
  as.numeric(lower <= x && x <= higher)
}

print.mean.sd <- function(log.thetas){
  print(mean(log.thetas), sd(log.thetas))
}

runSimulation <- function(task.num, job.id, num.jobs, pair.nums, theta.draws, Y.draws, get.mu, get.sigma, simTheta, w, write.w = FALSE){
  pair.num = ceiling(job.id / (num.jobs / pair.nums))
  
  num.groups = (num.jobs / pair.nums)
  size.groups = (theta.draws / num.groups)
  
  theta.group.num = ((job.id - 1) %% num.groups) + 1
  
  for(theta.group.offset in 1:size.groups){

    log.theta = simTheta(pair.num, J)
    theta.num = size.groups * (theta.group.num - 1) + theta.group.offset    

    print(sprintf("Running pair%d theta%d", pair.num, theta.num))
    
    is.covered.95 = NULL
    is.covered.68 = NULL
    
    for(Y.offset in 1:Y.draws){
      Y = simYgivenTheta(exp(log.theta), w, N)
      res = poisson.logn.mcmc(Y, w, mu0=get.mu(pair.num), sigmasq0=get.sigma(pair.num))

      apply(res$logTheta, print.mean.sd)
      
      CI.95.lower = apply(res$logTheta, 1, quantile, probs=.025, names=FALSE)
      CI.68.lower = apply(res$logTheta, 1, quantile, probs=.16, names=FALSE)

      CI.95.higher = apply(res$logTheta, 1, quantile, probs=.975, names=FALSE)
      CI.68.higher = apply(res$logTheta, 1, quantile, probs=.84, names=FALSE)

      is.covered.95.i = mapply(is.covered, CI.95.lower, CI.95.higher, log.theta)
      is.covered.68.i = mapply(is.covered, CI.68.lower, CI.68.higher, log.theta)

      if (is.null(is.covered.95)) {
        stopifnot(is.null(is.covered.68))
        is.covered.95 = is.covered.95.i
        is.covered.68 = is.covered.68.i
      } else {
        is.covered.95 = cbind(is.covered.95, is.covered.95.i)
        is.covered.68 = cbind(is.covered.68, is.covered.68.i)
      }
    }
    
    cover.95 = apply(is.covered.95, 1, mean)
    cover.68 = apply(is.covered.68, 1, mean)
    
    theta.cover = cbind(log.theta, cover.95, cover.68)
    save(theta.cover, file=getOutFileName(pair.num, theta.num, task.num))
    
    if(write.w){
      w.cover = cbind(w, cover.95, cover.68)
      save(w.cover, file=getOutFileName(pair.num, theta.num, task.num, TRUE))
    }

    gc()      
  }
}

aggregate.cover <- function(pair.nums, theta.draws, task.num, write.w = FALSE) {
  for(pair.num in 1:pair.nums){
    for(theta.draw in 1:theta.draws){
      load(getOutFileName(pair,theta.draw,task.num))
      write.table(theta.cover, 
                  file=sprintf("keskici_wxiao_ps2_task%d_par%d_theta.dat", task.num, pair.num), 
                  append=TRUE,
                  row.names=FALSE,
                  col.names=FALSE)
      
      if(write.w){
        load(getOutFileName(pair,theta.draw,task.num,TRUE))
        write.table(w.cover,
                    file=sprintf("keskici_wxiao_ps2_task%d_par%d_w.dat", task.num, pair.num), 
                    append=TRUE,
                    row.names=FALSE,
                    col.names=FALSE)
      }
    }
    plot.coverage(sprintf("keskici_wxiao_ps2_task%d_par%d_theta.dat", task.num, pair.num), task.num, pair.num)
    if(write.w){
      plot.coverage(sprintf("keskici_wxiao_ps2_task%d_par%d_theta.dat", task.num, pair.num), task.num, 2*(pair.num - 1) + 2, 
                    type = 'log w_j')      
    }
  }
}

# Plot coverage graphs
plot.coverage = function(file, tasknum, plotnum, type= "log theta_j"){
  a = as.matrix(read.table(file))
  if (type == 'log w_j'){
    a[,1] = log(a[,1])
  }
  a.68 = cbind(a[,1], a[,3])
  a.95 = cbind(a[,1],a[,2])
  smoothingSpline = smooth.spline(a[,1], a[,3],spar = 1)
  smoothingSpline95 = smooth.spline(a[,1], a[,2], spar = 1)
  output.format = "keskici_wxiao_ps2_task%d_plot%d.png"
  name = sprintf(output.format, tasknum, plotnum)
  png(name)
  
  par(mfcol=c(2,2))
  #Plot 68% CI Coverage
  plot(a.68, xlab = type, ylab = "Coverage", main = paste("68% CI Coverage for", type))
  lines(smoothingSpline, col='red', lwd=2)
  #Plot 95% CI Coverage
  plot(a.95, xlab = type, ylab = "Coverage", main = paste("95% CI Coverage for ", type))
  lines(smoothingSpline95, col='red', lwd=2)
  dev.off()
}
