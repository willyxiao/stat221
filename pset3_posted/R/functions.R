sim.alg = function(alg){
  risk.list = c()
  theta.t = theta.0
  
  for(t in 1:NSIMS){
    theta.t.1 = alg(t, theta.t)
    theta.t = theta.t.1
    risk = find.risk(theta.t.1)
    print(risk)
    risk.list = c(risk.list, risk)    
  }

  risk.list
}
