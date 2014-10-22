source('task2c_runner.R')

for(alg.name in c('sgd', 'asgd', 'implicit')){
  plot.name = paste(alg.name, toString(theta.run), "plot", sep="_")
  print(plot.name)
  pdf(paste(plot.name, "pdf", sep="."))
  load(file.name(alg.name, 1))
  plot(nlist, log10(res$bias), "l", lwd=2, col=colors()[25], ann=FALSE)
  
  for(a.id in 2:length(a.tests)){
    load(file.name(alg.name, a.id))
    lines(nlist, log10(res$bias), lwd=2, col=colors()[(a.id-1)*25])
  }

  title(main=plot.name, 
        xlab="n", 
        ylab="bias")
  
  legend('topright', legend=a.tests, 
         lty=c(1,1,1,1,1,1), 
         col=colors()[(1:length(a.tests)*25)])
  
  dev.off()
}
