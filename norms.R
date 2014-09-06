library(RUnit)
library(logging)
basicConfig()
addHandler(writeToFile, file="sample.log", level="WARNING")

# norms.R
# Functions to compute the L-2 norm of a column vector
# ||x|| = sqrt(x_1^2 + x_2^2 + ...)
# First implementation: Iterates over all elements
normOne <- function(x) {
  n = length(x)
  if(n==0) return(0)
  sum = 0
  if(length(unique(x)) == 1)
    logwarn("List of length %d has same elements", length(x))
  for(i in 1:n)
    sum = sum + x[i]^2
  return(sqrt(sum))
}

# Second implementation: recursive definition
normTwo <- function(x) {
  n = length(x)
  if(n==0) return(0)
  if(n==1)
    return(abs(x[1]))
  sqrt(normTwo(head(x, n-1))^2 + x[n]^2)
}
# Third implementation: Uses vector multiplication x’x
normThree <- function(x) {
  n = length(x)
  if(n==0) return(0)
  sqrt(t(x) %*% x)[1][1]
}

# Function that tests a "norm" function
norm.checks <- function(norm) {
  print("Starting tests..")
  checkTrue(is_near(norm(c()), 0))
  checkTrue(is_near(norm(c(3,4)), 5))
  x = rnorm(20)
  checkTrue(is_near(norm(2 * x), 2 * norm(x), tol=1e-6))
}
is_near <- function(x, x0, tol=1e-9) {
  return(abs(x-x0) < tol)
}

test.norm.one <- function() {
  norm.checks(normOne)
}
test.norm.two <- function() {
  norm.checks(normTwo)
}
test.norm.three <- function() {
  norm.checks(normThree)
}

run.tests <- function() {
  test.suite <- defineTestSuite(name="norm", dirs=".", testFileRegexp="norms.R")
  test.result <- runTestSuite(test.suite)
  printTextProtocol(test.result)
}

profile.norms <- function() {
  Rprof(filename="profile.txt",memory.profiling=T)
  for (i in 1:1000) {
    x <- rnorm(100)
    normOne(x)
    normTwo(x)
    normThree(x)
  }
  Rprof(NULL)
}

