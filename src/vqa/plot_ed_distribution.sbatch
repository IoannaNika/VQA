#!/bin/sh
#SBATCH --partition=general
#SBATCH --qos=medium
#SBATCH --time=05:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1024
#SBATCH --mail-type=END
#SBATCH --gres=gpu:0 

# echo "srun python plot_ed_distribution.py --ref_set /tudelft.net/staff-umbrella/ViralQuasispecies/inika/Read_simulators/data/triplet_dataset_primers_template_excl_neg_identicals_327896/samples.tsv --directory /tudelft.net/staff-umbrella/ViralQuasispecies/inika/Read_simulators/data/triplet_dataset_primers_template_excl_neg_identicals_327896/"
# srun python plot_ed_distribution.py --ref_set /tudelft.net/staff-umbrella/ViralQuasispecies/inika/Read_simulators/data/triplet_dataset_primers_template_excl_neg_identicals_327896/samples.tsv --directory /tudelft.net/staff-umbrella/ViralQuasispecies/inika/Read_simulators/data/triplet_dataset_primers_template_excl_neg_identicals_327896/
# echo "srun python models/NaiveBayesClassifier.py"
# srun python models/NaiveBayesClassifier.py



# echo "srun python run_naive_model_predictions.py --virus_name SARS-CoV-2 --dataset /tudelft.net/staff-umbrella/ViralQuasispecies/inika/Read_simulators/data/lumc_data/lumc_dataset.tsv"
# srun python run_naive_model_predictions.py --virus_name SARS-CoV-2 --dataset /tudelft.net/staff-umbrella/ViralQuasispecies/inika/Read_simulators/data/lumc_data/lumc_dataset.tsv

########## REV COMPLEMENT CORRECTED ##############

# echo "srun python plot_ed_distribution.py --ref_set /tudelft.net/staff-umbrella/ViralQuasispecies/inika/Read_simulators/data/tuples_pacbio_sars_cov_2_rev_compl_more/dataset/samples.tsv --directory /tudelft.net/staff-umbrella/ViralQuasispecies/inika/Read_simulators/data/tuples_pacbio_sars_cov_2_rev_compl_more/dataset"
# srun python plot_ed_distribution.py --ref_set /tudelft.net/staff-umbrella/ViralQuasispecies/inika/Read_simulators/data/tuples_pacbio_sars_cov_2_rev_compl_more/dataset/samples.tsv --directory /tudelft.net/staff-umbrella/ViralQuasispecies/inika/Read_simulators/data/tuples_pacbio_sars_cov_2_rev_compl_more/dataset
# echo "srun python models/NaiveBayesClassifier.py"
# srun python models/NaiveBayesClassifier.py

echo "srun python run_naive_model_predictions.py --virus_name SARS-CoV-2 --dataset /tudelft.net/staff-umbrella/ViralQuasispecies/inika/Read_simulators/data/lumc_rev_compl/all_test_pairs.tsv"
srun python run_naive_model_predictions.py --virus_name SARS-CoV-2 --dataset /tudelft.net/staff-umbrella/ViralQuasispecies/inika/Read_simulators/data/lumc_rev_compl/all_test_pairs.tsv
