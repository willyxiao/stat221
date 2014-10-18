sim.alg = function(alg, theta.0, nsims=2e3){
  risk.list = c()
  theta.t = theta.0
  
  for(t in 1:nsims){
    theta.t.1 = alg(t, theta.t)
    theta.t = theta.t.1
    risk = find.risk(theta.t.1)
    risk.list = c(risk.list, risk)    
  }

  risk.list
}
