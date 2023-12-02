#!/bin/sh

#SBATCH --partition=general
#SBATCH --qos=medium
#SBATCH --time=7:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=4096
#SBATCH --mail-type=END

python make_siamese_reference_set.py
