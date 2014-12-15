#! /usr/bin/env Rscript

#SBATCH -n 1
#SBATCH -t 60
#SBATCH -o dump/%A_%a.out
#SBATCH -e dump/%A_%a.err
#SBATCH -p stats
#SBATCH --mem-per-cpu=1000
#SBATCH --mail-user=willy@chenxiao.us
#SBATCH --mail-type=ALL
#SBATCH -a 1-26

source("data_generator.R")

job.id = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
run.normal.em(job.id)

