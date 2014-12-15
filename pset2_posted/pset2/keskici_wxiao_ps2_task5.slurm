#! /usr/bin/env Rscript

#SBATCH -n 1                            # (Max) number of tasks per job, for R usually 1
#SBATCH -t 60
#SBATCH -o dump/%A_%a.out                     # File for the standard output
#SBATCH -e dump/%A_%a.err                    # File for the standard error
#SBATCH --open-mode=append              # Append to standard output and error files
#SBATCH -p stats               # Partition to use
#SBATCH --mem-per-cpu=4096              # Memory required per CPU, in MegaBytes
#SBATCH --mail-user=willy@chenxiao.us   # Where to send mail
#SBATCH --mail-type=ALL                 # When to send mail
#SBATCH -a 1-12                          # Test

source("keskici_wxiao_ps2_task5.R")

JOB_ARRAYS = 12

if(Sys.getenv("SLURM_JOB_ID") != ""){
  job.id = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
  run.task5(job.id, JOB_ARRAYS)
}
