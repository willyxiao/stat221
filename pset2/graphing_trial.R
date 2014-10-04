library(tools)

plot.coverage = function(file, type= "log theta_j"){
  a = as.matrix(read.table(file))
  a.68 = cbind(a[,1], a[,3])
  a.95 = cbind(a[,1],a[,2])
  
  smoothingSpline = smooth.spline(a[,1], a[,3],spar = 1)
  smoothingSpline95 = smooth.spline(a[,1], a[,2], spar = 1)
  name = paste(file_path_sans_ext(file),".png", sep = "")
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

plot.coverage('keskici_wxiao_ps2_task3_par2_theta.dat')