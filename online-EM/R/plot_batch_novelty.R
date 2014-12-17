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

#plot normal em
plot.fig(res, 5, names, indices, "novelty_normal_em.pdf", 1e6)