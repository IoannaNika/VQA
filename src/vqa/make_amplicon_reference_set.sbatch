#!/bin/sh

#SBATCH --partition=general
#SBATCH --qos=short
#SBATCH --time=3:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=5120
#SBATCH --mail-type=END
#SBATCH --gres=gpu:0 

python make_amplicon_reference_set.py
