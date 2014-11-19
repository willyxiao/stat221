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
# means1130 = c()
# vars1130 = c()
# means1530 = c()
# vars1530 = c()
# for(i in 1:16){
#   item = get(sprintf("x%d", i))
#   means1130 = c(means1130, log(mean(item[,2][134:144]), base=10))
#   vars1130 = c(vars1130, log(var(item[,2][134:144]), base = 10))
#   
#   means1530 = c(means1530, log(mean(item[,2][182:192]), base=10))
#   vars1530 = c(vars1530, log(var(item[,2][182:192]), base = 10))
#   
# }
# pdf("keskici_wxiao_fig4_2router.pdf")
# par(mfrow=c(1,2))
# plot(means1130, vars1130, xlab="log10(mean)", ylab="log10(var)", main="Time 11:30")
# abline(lm(vars1130~means1130))
# 
# plot(means1530, vars1530, xlab="log10(mean)", ylab="log10(var)", main="Time 15:30")
# abline(lm(vars1530~means1530))
# dev.off()


names2 = c("dst router5", "dst r4-local", "dst switch", "dst r4-others",
          "dst gw1", "dst gw2", "dst gw3", "dst gw-others", "total",
          "gw-others->router5", "gw-others->r4-local", "gw-others->switch",
          "gw-others->r4-others", "gw-others->gw1", "gw-others->gw2",
          "gw-others->gw3", "gw-others->gw-others", "origin gw-others",
          "gw3->router5", "gw3->r4-local", "gw3->switch", "gw3->r4-others",
          "gw3->gw1", "gw3->gw2", "gw3->gw3", "gw3->gw-others", "origin gw-others",
          "gw2->router5", "gw2->r4-local", "gw2->switch", "gw2->r4-others",
          "gw2->gw1", "gw2->gw2", "gw2->gw3", "gw2->gw-others", "origin gw2",
          "gw1->router5", "gw1->r4-local", "gw1->switch", "gw1->r4-others",
          "gw1->gw1", "gw1->gw2", "gw1->gw3", "gw1->gw-others", "origin gw1",
          "r4-others->router5", "r4-others->r4-local", "r4-others->switch", "r4-others->r4-others",
          "r4-others->gw1", "r4-others->gw2", "r4-others->gw3", "r4-others->gw-others", 
          "origin r4-others", "switch->router5", "switch->r4-local", "switch->switch", "switch->r4-others",
          "switch->gw1", "switch->gw2", "switch->gw3", "switch->gw-others", "origin switch",
          "r4-local->router5", "r4-local->r4-local", "r4-local->switch", "r4-local->r4-others",
          "r4-local->gw1", "r4-local->gw2", "r4-local->gw3", "r4-local->gw-others", "origin r4-local",
          "router5->router5", "router5->r4-local", "router5->switch", "router5->r4-others", 
          "router5->gw1", "router5->gw2", "router5->gw3", "router5->gw-others", "origin router5")

#generate A matrix
A = matrix(0, ncol=64, nrow=15)
#populate A matrix
for(i in 1:15){
  if(i <= 8){
    start = 8*(i - 1)
    for(j in 1:8){
      A[i, start + j] = 1
    }
  }
  else{
    start = i - 8
    for(j in 1:8){
      A[i, start + (j-1)*8] = 1
    }
    
  } 
}

#1.9
data = x1[,2]
for(i in 2:15){
  item = get(sprintf("x%d", i))
  data = cbind(data, item[,2])
}


#router2.fig5.dat = locally_iid_EM(data, 2, A)
#save(router2.fig5.dat, file="r2fig5.RData")

router2.fig6.dat = smoothed_EM(data, 2, A)
save(router2.fig6.dat, file="r2fig6.RData")

