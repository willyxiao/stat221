library(mvtnorm)
dat = read.csv("1router_allcount.dat")

src.fddi = dat[,3][17]
src.switch = dat[,3][18]
src.local = dat[,3][19]
src.copr = dat[,3][20]

dst.fddi = dat[,3][21]
dst.switch = dat[,3][22]
dst.local = dat[,3][23]
dst.copr = dat[,3][24]

x1 = (dat[dat[,3] == src.fddi,])[,'value']
x2 = (dat[dat[,3] == src.switch,])[,'value']
x3 = (dat[dat[,3] == src.local,])[,'value']
x4 = (dat[dat[,3] == src.copr,])[,'value']
x5 = (dat[dat[,3] == dst.fddi,])[,'value']
x6 = (dat[dat[,3] == dst.switch,])[,'value']
x7 = (dat[dat[,3] == dst.local,])[,'value']
x8 = (dat[dat[,3] == dst.copr,])[,'value']

#splits data vector into 26 chunks
split.dat = function(router.dat){
  router.dat = router.dat[1:286]
  dat.matrix = matrix(router.dat, nrow = 11, ncol = 26)
  means = colMeans(dat.matrix)
  vars = apply(dat.matrix, 2, var)
  return(list(means = means, vars = vars))
}

#generates data for iteration iter (1-26)
generate.data = function(n, iter){
  mu = c()
  sigma = diag(NaN, 8)
  
  for(i in 1:8){
    x = get(sprintf("x%d", i))
    res = split.dat(x)
    mu = c(mu, res$means[iter])
    sigma[i,i] = res$vars[iter]
  }  
  rmvnorm(n, mu, sigma)
}