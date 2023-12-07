#!/bin/sh

#SBATCH --partition=general
#SBATCH --qos=medium
#SBATCH --time=36:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=15360
#SBATCH --mail-type=END
#SBATCH --gres=gpu:0 

python make_siamese_reference_set.py
