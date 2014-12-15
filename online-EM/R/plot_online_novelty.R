source('part3.R')

res = c()
for(i in 1:26){
  estimates = online.em.each(generate.data(500, i), A)$m
  res = rbind(res, c(estimates, NaN)) #add relevant values
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

#plot online em results
plot.fig(res, 5, names, indices, "novelty_online_em_500iters.pdf", 1e6)
