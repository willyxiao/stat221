#! /usr/bin/env Rscript
#SBATCH -n 1                            # (Max) number of tasks per job, for R usually 1
#SBATCH -o out.txt                      # File for the standard output
#SBATCH -e err.txt                      # File for the standard error
#SBATCH --open-mode=append              # Append to standard output and error files
#SBATCH -p serial_requeue               # Partition to use
#SBATCH --mem-per-cpu=4096              # Memory required per CPU, in MegaBytes
#SBATCH --mail-user=willy@chenxiao.us   # Where to send mail
#SBATCH --mail-type=ALL                 # When to send mail

noisy_cor <- function(n, true_cor, noise_level, to_file = FALSE) {
  
  ## random seed
  seed <- fracSec()
  set.seed(seed)
  
  sigma <- matrix(c(1,true_cor,true_cor,1), nrow=2)
  D <- rmvnorm(n,sigma=sigma)
  D[,1] <- D[,1] + rnorm(n, sd=sqrt(noise_level))
  D[,2] <- D[,2] + rnorm(n, sd=sqrt(noise_level))
  
  res = cor(D[,1], D[,2])
  if(to_file) {
    outputfile <- "test_output.txt"
    save(res, seed, file=outputfile)
  }
}
