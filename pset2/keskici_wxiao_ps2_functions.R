simYgivenTheta <- function(theta, w, N) {
  stopifnot(length(theta) == length(w))  
  J = length(theta)
  matrix(rpois(N*J, w*theta), nrow=J,ncol=N)
}

simThetagivenMuSigma <- function(mu, sigma, J) {
  rnorm(J, mu, sigma)
}

# same as getFileName but in the out file and called "out"
getOutFileName <- function(i,j, k, task = NULL){
  if(is.null(task)){
    file.name.str <- "out/out_%d_%d_%d.RData"    
    sprintf(file.name.str, i, j, k)
  } else {
    file.name.str <- "out/out%d_%d_%d_%d.RData"
    sprintf(file.name.str, task, i, j, k)
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
