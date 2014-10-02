#! /usr/bin/env Rscript

#FIXME How do we request 12 machines? Is it based on number of tasks per job?

#SBATCH -n 1                            # (Max) number of tasks per job, for R usually 1
#SBATCH -t 3
#SBATCH -o %A_%a.out                     # File for the standard output
#SBATCH -e %A_%a.err                    # File for the standard error
#SBATCH --open-mode=append              # Append to standard output and error files
#SBATCH -p stats               # Partition to use
#SBATCH --mem-per-cpu=4096              # Memory required per CPU, in MegaBytes
#SBATCH --mail-user=willy@chenxiao.us   # Where to send mail
#SBATCH --mail-type=ALL                 # When to send mail
#SBATCH -a 1-4                          # Test

source("keskici_wxiao_ps2_functions.R")

mcmc.NSIMS = 2

stopif(length(mu) != length(sigma))

if(Sys.getenv("SLURM_JOB_ID") != ""){

} else { # local testing
  for(i in 1:length(mu)){
    for(j in 1:theta.NSIMS) {
      load(getFileName(i,j))
      res = poisson.logn.mcmc(getObjectName(i,j), w, mu0=mu[i], sigmasq0=sigma[i]**2)
      browser()  
  }
  }  
}

