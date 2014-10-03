#! /usr/bin/env Rscript

#SBATCH -n 1                            # (Max) number of tasks per job, for R usually 1
#SBATCH -t 3
#SBATCH -o %A_%a.out                     # File for the standard output
#SBATCH -e %A_%a.err                    # File for the standard error
#SBATCH --open-mode=append              # Append to standard output and error files
#SBATCH -p stats               # Partition to use
#SBATCH --mem-per-cpu=4096              # Memory required per CPU, in MegaBytes
#SBATCH --mail-user=willy@chenxiao.us   # Where to send mail
#SBATCH --mail-type=ALL                 # When to send mail
#SBATCH -a 1-12                          # Test

source("keskici_wxiao_ps2_functions.R")
source("poissonLogN_MCMC.R")

# roughly 15 seconds per simulations (ndraws=1000) ==> 2880 total simulations on 12 nodes for 1 hour

stopifnot(length(mu) == length(sigma))

JOB_ARRAYS = 12

if(Sys.getenv("SLURM_JOB_ID") != ""){
  id = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))


  pair.num = floor((id - 1) / length(mu)) + 1

  num.groups = (JOB_ARRAYS / length(mu))
  size.groups = (theta.NSIMS / num.groups)

  theta.group.num = ((id - 1) %% num.groups) + 1

  for(theta.group.offset in 1:size.groups){

    log.theta = simThetagivenMuSigma(mu[pair.num], sigma[pair.num], J)
    Y = simYgivenTheta(exp(log.theta), w, N)

    theta.num = (theta.NSIMS / size.groups) * theta.group.num + theta.group.offset
#    load(getFileName(pair.num, theta.num))
#    Y = get(getObjectName(pair.num, theta.num))
    res = poisson.logn.mcmc(Y, w, mu0=mu[pair.num], sigmasq0=sigma[pair.num]**2)
    assign(getOutObjectName(pair.num, theta.num), res)
    save(list=c(getOutObjectName(pair.num, theta.num)), file=getFileName(pair.num,theta.num))
  }

} else { # local testing
  for(i in 1:length(mu)){
    for(j in 1:theta.NSIMS) {
      load(getFileName(i,j))
      Y = get(getObjectName(i,j))
      res = poisson.logn.mcmc(Y, w, mu0=mu[i], sigmasq0=sigma[i]**2)
    }
  }  
}

