source('keskici_wxiao_functions.R')
dat = read.csv("2router_linkcount.dat")

ori.router5 = dat[,3][1]
ori.local = dat[,3][2]
ori.switch = dat[,3][3]
ori.r4-others = dat[,3][4]
ori.gw1 = dat[,3][5]
ori.gw2 = dat[,3][6]
ori.gw3 = dat[,3][7]
ori.gw-others = dat[,3][8]

dest.router5 = dat[,3][9]
dest.local = dat[,3][10]
dest.switch = dat[,3][11]
dest.r4-others = dat[,3][12]
dest.gw1 = dat[,3][13]
dest.gw2 = dat[,3][14]
dest.gw3 = dat[,3][15]
dest.gw-others = dat[,3][16]






