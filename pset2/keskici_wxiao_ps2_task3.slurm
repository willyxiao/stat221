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

if(Sys.getenv("SLURM_JOB_ID") != ""){
  library(pset2)

  mu <-     c(1,6, 2.5, 5.2, 4.9)
  sigma <-  c(0.7, 1.3, 1.3, 1.6)

  id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
  
  hello.world(mu[id], sigma[id])
}

