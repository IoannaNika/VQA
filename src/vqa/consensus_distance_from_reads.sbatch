#!/bin/sh

#SBATCH --partition=general,insy
#SBATCH --qos=long
#SBATCH --time=168:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=50GB
#SBATCH --output=slurm_lumc_analysis_%j.out # Set name of output log. %j is the Slurm jobId
#SBATCH --error=slurm_lumc_analysis_%j.err # Set name of error log. %j is the Slurm jobId
#SBATCH --gres=gpu:0 

srun python consensus_distance_from_reads.py

