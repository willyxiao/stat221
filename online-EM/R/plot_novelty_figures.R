source("ps5_functions.R")

res = c()
for(i in 1:26){
  filename = sprintf("normal_em_window%d.RData", i)
  load(filename) #loads variable em.fig5.dat into memory
  res = rbind(res, colMeans(em.fig5.dat)) #add relevant values
}

#add in dest totals
for(i in 1:4){
  res = cbind(res, res[,i] + res[,i + 4] + res[,i + 8] + res[,i + 12])
}

#add in total
total = res[,1]
for(i in 2:16){
  total = total + res[, i]
}
res = cbind(res, total)

#add in origin totals
for(i in c(13, 9, 5, 1)){
  res = cbind(res, res[,i] + res[,i + 1] + res[,i + 2] + res[,i + 3])
}

# names = c("destination fddi", "destination switch", "destination local",
#           "destination corp", "total", "corp->fddi", "corp->switch", "corp->local",
#           "corp->corp", "origin corp", "local->fddi", "local->switch", "local->local", "local->corp",
#           "origin local", "switch->fddi", "switch->switch", "switch->local", "switch->corp",
#           "origin switch", "fddi->fddi", "fddi->switch", "fddi->local", "fddi->corp", "origin fddi")
# indices = c(18, 19, 20, 21, 22, 13, 14, 15, 16, 23, 9, 10, 11, 12, 24, 5, 6, 7, 8, 25, 1, 2, 3, 4, 26)

#plot normal em
plot.fig(res, 5, names, indices, "novelty_normal_em.pdf", 1e6)

