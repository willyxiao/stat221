#!/bin/bash
#SBATCH -J parallel-testing # name for job array
#SBATCH -o dump  #Standard output
#SBATCH -e all.err #Standard error
#SBATCH -p stats  #Partition
#SBATCH -t 24:0:00 #Running time of one hour
#SBATCH --mem-per-cpu 5000 #Memory request
#SBATCH -n 1 #Number of cores
#SBATCH -N 1 #All cores on one machine
#SBATCH -a 1-50 #array of 50 jobs with ids 1-50

Rscript projectBatchEM.R $SLURM_ARRAY_TASK_ID
