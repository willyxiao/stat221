source("keskici_wxiao_ps2_functions.R")

for(i in 1:length(mu)){
  for(j in 1:theta.NSIMS) {
    log.theta = simThetagivenMuSigma(mu[i], sigma[i], J)
    Y = simYgivenTheta(exp(log.theta), w, N)
    assign(getObjectName(i,j), Y)
    save(list=c(getObjectName(i,j)), file=getFileName(i,j))
  }
}
