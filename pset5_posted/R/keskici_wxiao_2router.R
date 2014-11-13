source('keskici_wxiao_functions.R')
dat = read.csv("2router_linkcount.dat")

ori.router5 = dat[,3][1]
ori.local = dat[,3][2]
ori.switch = dat[,3][3]
ori.r4_others = dat[,3][4]
ori.gw1 = dat[,3][5]
ori.gw2 = dat[,3][6]
ori.gw3 = dat[,3][7]
ori.gw_others = dat[,3][8]
dest.router5 = dat[,3][9]
dest.local = dat[,3][10]
dest.switch = dat[,3][11]
dest.r4_others = dat[,3][12]
dest.gw1 = dat[,3][13]
dest.gw2 = dat[,3][14]
dest.gw3 = dat[,3][15]
dest.gw_others = dat[,3][16]

x1 = (dat[dat[,3] == ori.router5,])[,c('time', 'value')]
x2 = (dat[dat[,3] == ori.local,])[,c('time', 'value')]
x3 = (dat[dat[,3] == ori.switch,])[,c('time', 'value')]
x4 = (dat[dat[,3] == ori.r4_others,])[,c('time', 'value')]
x5 = (dat[dat[,3] == ori.gw1,])[,c('time', 'value')]
x6 = (dat[dat[,3] == ori.gw2,])[,c('time', 'value')]
x7 = (dat[dat[,3] == ori.gw3,])[,c('time', 'value')]
x8 = (dat[dat[,3] == ori.gw_others,])[,c('time', 'value')]
x9 = (dat[dat[,3] == dest.router5,])[,c('time', 'value')]
x10 = (dat[dat[,3] == dest.local,])[,c('time', 'value')]
x11 = (dat[dat[,3] == dest.switch,])[,c('time', 'value')]
x12 = (dat[dat[,3] == dest.r4_others,])[,c('time', 'value')]
x13 = (dat[dat[,3] == dest.gw1,])[,c('time', 'value')]
x14 = (dat[dat[,3] == dest.gw2,])[,c('time', 'value')]
x15 = (dat[dat[,3] == dest.gw3,])[,c('time', 'value')]
x16 = (dat[dat[,3] == dest.gw_others,])[,c('time', 'value')]

#1.2
#via inspection we know relevent indices are
#134 135 136 137 138 139 140 141 142 143 144
#for the window around 11:30
#and   182 183 184 185 186 187 188 189 190 191 192
#for the window around 15:30
means1130 = c()
vars1130 = c()
means1530 = c()
vars1530 = c()
for(i in 1:16){
  item = get(sprintf("x%d", i))
  means1130 = c(means1130, log(mean(item[,2][134:144]), base=10))
  vars1130 = c(vars1130, log(var(item[,2][134:144]), base = 10))
  
  means1530 = c(means1530, log(mean(item[,2][182:192]), base=10))
  vars1530 = c(vars1530, log(var(item[,2][182:192]), base = 10))
  
}
pdf("keskici_wxiao_fig4_2router.pdf")
par(mfrow=c(1,2))
plot(means1130, vars1130, xlab="log10(mean)", ylab="log10(var)", main="Time 11:30")
abline(lm(vars1130~means1130))

plot(means1530, vars1530, xlab="log10(mean)", ylab="log10(var)", main="Time 15:30")
abline(lm(vars1530~means1530))
dev.off()





