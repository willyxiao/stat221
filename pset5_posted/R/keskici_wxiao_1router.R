source('keskici_wxiao_functions.R')
dat = read.csv("1router_allcount.dat")

#1.1
src.fddi = dat[,3][17]
src.switch = dat[,3][18]
src.local = dat[,3][19]
src.copr = dat[,3][20]

dst.fddi = dat[,3][21]
dst.switch = dat[,3][22]
dst.local = dat[,3][23]
dst.copr = dat[,3][24]

x1 = (dat[dat[,3] == src.fddi,])[,c('time', 'value')]
x2 = (dat[dat[,3] == src.switch,])[,c('time', 'value')]
x3 = (dat[dat[,3] == src.local,])[,c('time', 'value')]
x4 = (dat[dat[,3] == src.copr,])[,c('time', 'value')]
x5 = (dat[dat[,3] == dst.fddi,])[,c('time', 'value')]
x6 = (dat[dat[,3] == dst.switch,])[,c('time', 'value')]
x7 = (dat[dat[,3] == dst.local,])[,c('time', 'value')]
x8 = (dat[dat[,3] == dst.copr,])[,c('time', 'value')]

#pdf("keskici_wxiao_fig2.pdf")
#par(mfrow=c(4,1))
#plot(x4, ylim=c(0,1e6), main="corp", ylab="bytes/sec")
#lines(x4[,1], x4[,2])
#lines(x8[,1], x8[,2])
#plot(x3, ylim=c(0,1e6), main="local", ylab="bytes/sec")
#lines(x3[,1], x3[,2])
#lines(x7[,1], x7[,2])
#plot(x2, ylim=c(0,1e6), main="switch", ylab="bytes/sec")
#lines(x2[,1], x2[,2])
#lines(x6[,1], x6[,2])
#plot(x1, ylim=c(0,1e6), main="fddi", ylab="bytes/sec")
#lines(x1[,1], x1[,2])
#lines(x5[,1], x5[,2])
#dev.off()

#1.2
#via inspection we know relevent indices are
#134 135 136 137 138 139 140 141 142 143 144
#for the window around 11:32
#and   182 183 184 185 186 187 188 189 190 191 192
#for the window around 15:32
means1130 = c()
vars1130 = c()
means1530 = c()
vars1530 = c()
for(i in 1:8){
  item = get(sprintf("x%d", i))
  means1130 = c(means1130, log(mean(item[,2][134:144]), base=10))
  vars1130 = c(vars1130, log(var(item[,2][134:144]), base = 10))
  
  means1530 = c(means1530, log(mean(item[,2][182:192]), base=10))
  vars1530 = c(vars1530, log(var(item[,2][182:192]), base = 10))
  
}
# pdf("keskici_wxiao_fig4_1router.pdf")
# par(mfrow=c(1,2))
# plot(means1130, vars1130, xlab="log10(mean)", ylab="log10(var)", main="Time 11:30")
# abline(lm(vars1130~means1130))
# 
# plot(means1530, vars1530, xlab="log10(mean)", ylab="log10(var)", main="Time 15:30")
# abline(lm(vars1530~means1530))
# dev.off()

#1.4
data = x1[,2]
for(i in 2:7){
  item = get(sprintf("x%d", i))
  data = cbind(data, item[,2])
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

locally_iid_EM.each(data, 2, A)
  
