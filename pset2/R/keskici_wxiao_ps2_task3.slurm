#! /usr/bin/env Rscript
#SBATCH -n 1                            # (Max) number of tasks per job, for R usually 1
#SBATCH -o out.txt                      # File for the standard output
#SBATCH -e err.txt                      # File for the standard error
#SBATCH --open-mode=append              # Append to standard output and error files
#SBATCH -p serial_requeue               # Partition to use
#SBATCH --mem-per-cpu=4096              # Memory required per CPU, in MegaBytes
#SBATCH --mail-user=willy@chenxiao.us   # Where to send mail
#SBATCH --mail-type=ALL                 # When to send mail


