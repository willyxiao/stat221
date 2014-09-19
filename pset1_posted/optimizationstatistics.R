## Optimization examples from R's optim() documentation.
# Panos Toulis ptoulis@fas.harvard.edu
rm(list=ls())
require(graphics)

# Runs some optimization methods one by one and prints the iterates path.
# We use Rosenbrock's "banana" function. This is an adaptation of the 
# the optim() documentation.
#
# f(x) = M * (x2 - x1^2)^2 + (1-x1)^2
#
# M = some scale parameter. 
# init.point.scale = higher values push the starting point further away.
banana <- function(M=100, init.point.scale=1) {
  
  # XYpath = nx2 matrix where n=#iterates, row i has the (x1, x2)_i
  # i.e., the iterate (x1, x2) at i
  # Note: Using XYpath <<- ...   makes it a global variable
  XYpath <<- 0
  
  # Clears up the path matrix.
  clear.path <- function() {
    XYpath <<- matrix(0, nrow=0, ncol=2)
  }
  # Plots the path
  plot.path <- function(str) {
    m = nrow(XYpath)
    x0 = XYpath[,1]
    y0 = XYpath[,2]
    
    x1 = c(tail(x0, m-1), x0[m]+1e-5)
    y1 = c(tail(y0, m-1), y0[m] + 1e-5)
    
    arrows(x0=x0, y0=y0, 
           x1=x1, y1=y1, 
           lty=2, col="magenta", lwd=0.5, length=0.1)
  }
  
  # Optimization function (naming taken from optim doc)
  fr <- function(x) {   ## Rosenbrock Banana function
    x1 <- x[1]
    x2 <- x[2]
    M * (x2 - x1 * x1)^2 + (1 - x1)^2
  }
  # Same function but also adds point in the path.
  fr.addPath <- function(x) {   ## Rosenbrock Banana function
    XYpath <<- rbind(XYpath, x)
    return(fr(x))
  }
  
  # Calculate the contours.
  x = seq(-2, 2, length.out=100) # x1
  y = x  # actually x2
  # expand.grid gets all combinations of (x, y) i.e., (x1, x2)
  grid = expand.grid(x, y)
  df = as.data.frame(grid)
  # z = has values f(x1, x2)
  z = apply(df, 1, fr)
  z = matrix(z, nrow=100, byrow=F)
  
  # Create 2x3 plot area.
  par(mfrow=c(2,3))
  # Gradient of the banana function..
  grr <- function(x) { ## Gradient of 'fr'
    x1 <- x[1]
    x2 <- x[2]
    c(-4 * M * x1 * (x2 - x1 * x1) - 2 * (1 - x1),
      2 * M *      (x2 - x1 * x1))
  }
  
  # 1. Plot the contours.
  print(sprintf("Contours of banana function. min=%.2f", fr(c(1, 1))))
  contour(x, y, z, levels=c(1, 10, 20, 30, 40, 50, 100, 150, 200, 400))
  points(c(1), c(1), pch=3, col="magenta", cex=2, lwd=2)
  readline("Press KEY to run optim().")
  
  # 2. Check and plot the contours.
  p = c(-0.1, -1)
  grad.p = -0.8 * grr(p) / sqrt(sum(grr(p)^2))
  new.vec = c(p[1] + grad.p[1], p[2] + grad.p[2])
  #   abline(v=p[1])
  #   abline(h=p[2])
  arrows(p[1], p[2], new.vec[1], new.vec[2], length=0.15,col="red", lwd=2)
  # arrows(0, 0, grad.p[1], grad.p[2], length=0.15,col="red", lwd=2)
  
  point0 = init.point.scale * c(-2, -2)
  # 3. Run the optimization method, info about iterations --plot path.
  run.method <- function(str) {
    print(sprintf("Running %s", str))
    clear.path()
    t = NA
    if(str=="BFGS+gradient") {
      t = system.time({out = optim(point0, 
                                   fr.addPath, gr=grr, method = "BFGS")})
      
    } else {
    t = system.time({out = optim(point0, 
                                 fr.addPath, method = str)})
    }

    print(sprintf("Iterations = %d Convergence=%s", 
                  out$counts["function"], out$convergence==0))
    print(sprintf("Estimates = %s", paste(round(out$par, 2), collapse=", ")))
    contour(x, y, z, main=sprintf("Method=%s", str),
            levels=c(1, 10, 20, 30, 40, 50, 100, 150, 200, 400))
    plot.path()
    readline("Press KEY for next method.")
  }
  
  run.method("Nelder-Mead")
  run.method("BFGS")
#   run.method("BFGS+gradient")
  run.method("L-BFGS-B")
  run.method("CG")
  run.method("SANN")
}

# Wild function.
run.wild <- function() {
  ## "wild" function , global minimum at about -15.81515
  fw <- function (x)
    10*sin(0.3*x)*sin(1.3*x^2) + 0.00001*x^4 + 0.2*x+80
  plot(fw, -50, 50, n=1000, main = "optim() minimising 'wild function'")
  
  res.annealing <- optim(50, fw, method="SANN",
                         control=list(maxit=20000, temp=20, parscale=20))
  res.bfgs <- optim(50, fw, method="BFGS",
                         control=list(maxit=20000, temp=20, parscale=20))
  res.annealing
  res.bfgs
}

# Travelling-salesman problem
run.TS <- function() {
  library(stats) # normally loaded
  
  eurodistmat <- as.matrix(eurodist)
  
  distance <- function(sq) {  # Target function
    sq2 <- embed(sq, 2)
    sum(eurodistmat[cbind(sq2[,2],sq2[,1])])
  }
  
  genseq <- function(sq) {  # Generate new candidate sequence
    idx <- seq(2, NROW(eurodistmat)-1)
    changepoints <- sample(idx, size=2, replace=FALSE)
    tmp <- sq[changepoints[1]]
    sq[changepoints[1]] <- sq[changepoints[2]]
    sq[changepoints[2]] <- tmp
    sq
  }
  
  sq <- c(1:nrow(eurodistmat), 1)  # Initial sequence: alphabetic
  distance(sq)
  # rotate for conventional orientation
  loc <- -cmdscale(eurodist, add=TRUE)$points
  x <- loc[,1]; y <- loc[,2]
  s <- seq_len(nrow(eurodistmat))
  tspinit <- loc[sq,]
  
  plot(x, y, type="n", asp=1, xlab="", ylab="",
       main="initial solution of traveling salesman problem", axes = FALSE)
  arrows(tspinit[s,1], tspinit[s,2], tspinit[s+1,1], tspinit[s+1,2],
         angle=10, col="green")
  text(x, y, labels(eurodist), cex=0.8)
  
  set.seed(123) # chosen to get a good soln relatively quickly
  res <- optim(sq, distance, genseq, method = "SANN",
               control = list(maxit = 30000, temp = 2000, trace = TRUE,
                              REPORT = 500))
  res  # Near optimum distance around 12842
  
  tspres <- loc[res$par,]
  plot(x, y, type="n", asp=1, xlab="", ylab="",
       main="optim() 'solving' traveling salesman problem", axes = FALSE)
  arrows(tspres[s,1], tspres[s,2], tspres[s+1,1], tspres[s+1,2],
         angle=10, col="red")
  text(x, y, labels(eurodist), cex=0.8)
}
