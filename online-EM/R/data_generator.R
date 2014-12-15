library(mvtnorm)
source("ps5_functions.R")
dat = read.csv("1router_allcount.dat")
#get levels
src.fddi = dat[,3][17]
src.switch = dat[,3][18]
src.local = dat[,3][19]
src.corp = dat[,3][20]

dst.fddi = dat[,3][21]
dst.switch = dat[,3][22]
dst.local = dat[,3][23]
dst.corp = dat[,3][24]

fddi.fddi = dat[,3][1]
fddi.switch = dat[,3][2]
fddi.local = dat[,3][3]
fddi.corp = dat[,3][4]
switch.fddi = dat[,3][5]
switch.switch = dat[,3][6]
switch.local = dat[,3][7]
switch.corp = dat[,3][8]
local.fddi = dat[,3][9]
local.switch = dat[,3][10]
local.local = dat[,3][11]
local.corp = dat[,3][12]
corp.fddi = dat[,3][13]
corp.switch = dat[,3][14]
corp.local = dat[,3][15]
corp.corp = dat[,3][16]


x1 = (dat[dat[,3] == src.fddi,])[,'value']
x2 = (dat[dat[,3] == src.switch,])[,'value']
x3 = (dat[dat[,3] == src.local,])[,'value']
x4 = (dat[dat[,3] == src.corp,])[,'value']
x5 = (dat[dat[,3] == dst.fddi,])[,'value']
x6 = (dat[dat[,3] == dst.switch,])[,'value']
x7 = (dat[dat[,3] == dst.local,])[,'value']
x8 = (dat[dat[,3] == dst.corp,])[,'value']

y1 = (dat[dat[,3] == fddi.fddi,])[,'value']
y2 = (dat[dat[,3] == fddi.switch,])[,'value']
y3 = (dat[dat[,3] == fddi.local,])[,'value']
y4 = (dat[dat[,3] == fddi.corp,])[,'value']
y5 = (dat[dat[,3] == switch.fddi,])[,'value']
y6 = (dat[dat[,3] == switch.switch,])[,'value']
y7 = (dat[dat[,3] == switch.local,])[,'value']
y8 = (dat[dat[,3] == switch.corp,])[,'value']
y9 = (dat[dat[,3] == local.fddi,])[,'value']
y10 = (dat[dat[,3] == local.switch,])[,'value']
y11 = (dat[dat[,3] == local.local,])[,'value']
y12 = (dat[dat[,3] == local.corp,])[,'value']
y13 = (dat[dat[,3] == corp.fddi,])[,'value']
y14 = (dat[dat[,3] == corp.switch,])[,'value']
y15 = (dat[dat[,3] == corp.local,])[,'value']
y16 = (dat[dat[,3] == corp.corp,])[,'value']

#splits data vector into 26 chunks
split.dat = function(router.dat){
  router.dat = router.dat[1:286]
  dat.matrix = matrix(router.dat, nrow = 11, ncol = 26)
  means = colMeans(dat.matrix)
  vars = apply(dat.matrix, 2, var)
  return(list(means = means, vars = vars))
}

actual.window.means = function(iter=0){
  means = split.dat(y1)$mean
  for(i in 2:16){
    y = get(sprintf("y%d", i))
    res = split.dat(y)$means
    means = cbind(means, res)
  }
  if (iter == 0){
    means
  }
  else{
    means[iter,]
  }
}

#takes in generated router data (i.e. fddi->fddi)
#and returns total (i.e. src fddi, dest corp)
get.totals = function(router.counts){
  res = router.counts[,1] + router.counts[,2] + 
    router.counts[,3] + router.counts[,4]
  for(i in 1:3){
    res = cbind(res, router.counts[,(i*4) + 1] + router.counts[,(i*4) + 2] + 
                  router.counts[,(i*4) + 3] + router.counts[,(i*4) + 4])
  }
  for(i in 1:3){
    res = cbind(res, router.counts[,i] + router.counts[,i + 4] + 
                  router.counts[,i + 8] + router.counts[,i + 12])
  }
  res
}

#if we wanted to generate latent router and sum
#for totals
generate.data = function(n, iter){
  mu = c()
  sigma = diag(NaN, 16)
  for(i in 1:16){
    y = get(sprintf("y%d", i))
    res = split.dat(y)
    mu = c(mu, res$means[iter])
    sigma[i,i] = res$vars[iter]
  }  
  res = rmvnorm(n, mu, sigma)
  res[res<0] = 0
  get.totals(res)
}

#generate A matrix
A = matrix(0, ncol=16, nrow=7)
#populate A matrix
for(i in 1:7){
  if(i <= 4){
    start = 4*(i - 1)
    for(j in 1:4){
      A[i, start + j] = 1
    }
  }
  else{
    start = i - 4
    for(j in 1:4){
      A[i, start + (j-1)*4] = 1
    }
    
  }
  
}

#run normal EM on window "iteration"
#SLURM job calls this function
run.normal.em = function(iteration = 1){
  data = generate.data(n=500, iter = iteration)
  em.fig5.dat = locally_iid_EM(data, 2, A)
  filename = sprintf("normal_em_window%d.RData", iteration)
  save(em.fig5.dat, file=filename)
}

#generates data for iteration iter (1-26)
#dst.crp is not sampled b/c it's determined by
#others
# generate.data = function(n, iter){
#   mu = c()
#   sigma = diag(NaN, 7)
#   
#   for(i in 1:7){
#     x = get(sprintf("x%d", i))
#     res = split.dat(x)
#     mu = c(mu, res$means[iter])
#     sigma[i,i] = res$vars[iter]
#   }  
#   res = rmvnorm(n, mu, sigma)
#   res[res<0] = 0
#   res
# }