#!/bin/sh
#SBATCH --partition=general
#SBATCH --qos=short
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1024
#SBATCH --mail-type=END
#SBATCH --gres=gpu:0 


python plot_ed_distribution.py --ref_set /tudelft.net/staff-umbrella/ViralQuasispecies/inika/Read_simulators/data/triplet_dataset_test_100000/triplets.tsv --directory /tudelft.net/staff-umbrella/ViralQuasispecies/inika/Read_simulators/data/triplet_dataset_test_100000/  --triplet True
