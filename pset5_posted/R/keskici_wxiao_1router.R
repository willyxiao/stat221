source('keskici_wxaio_functions.R')
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

pdf("keskici_wxiao_fig2.pdf")
par(mfrow=c(4,1))
plot(x4, ylim=c(0,1e6), main="corp", ylab="bytes/sec")
lines(x4[,1], x4[,2])
lines(x8[,1], x8[,2])
plot(x3, ylim=c(0,1e6), main="local", ylab="bytes/sec")
lines(x3[,1], x3[,2])
lines(x7[,1], x7[,2])
plot(x2, ylim=c(0,1e6), main="switch", ylab="bytes/sec")
lines(x2[,1], x2[,2])
lines(x6[,1], x6[,2])
plot(x1, ylim=c(0,1e6), main="fddi", ylab="bytes/sec")
lines(x1[,1], x1[,2])
lines(x5[,1], x5[,2])
dev.off()

#1.2