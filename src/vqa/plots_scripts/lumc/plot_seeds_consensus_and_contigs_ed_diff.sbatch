#!/bin/sh
#SBATCH --partition=general,insy # Request partition. Default is 'general' 
#SBATCH --qos=long         # Request Quality of Service. Default is 'short' (maximum run time: 4 hours)
#SBATCH --time=168:00:00      # Request run time (wall-clock). Default is 1 minute
#SBATCH --ntasks=1          # Request number of parallel tasks per job. Default is 1
#SBATCH --cpus-per-task=1   # Request number of CPUs (threads) per task. Default is 1 (note: CPUs are always allocated to jobs per 2).
#SBATCH --mem=5120          # Request memory (MB) per node. Default is 1024MB (1GB). For multiple tasks, specify --mem-per-cpu instead
#SBATCH --output=slurm_%j.out # Set name of output log. %j is the Slurm jobId
#SBATCH --error=slurm_%j.err # Set name of error log. %j is the Slurm jobId
#SBATCH --gres=gpu:0 # Request 1 GPU

srun python plots_scripts/lumc/plot_seeds_consensus_and_contigs_ed_diff.py
