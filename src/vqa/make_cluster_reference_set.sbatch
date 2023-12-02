#!/bin/sh

#SBATCH --partition=general
#SBATCH --qos=medium
#SBATCH --time=3:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=40960
#SBATCH --mail-type=END

python make_cluster_reference_set.py
