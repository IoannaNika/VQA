python plots_scripts/HCV-1b/100/PB/90/plot_edit_distance_to_consensus.py --input_dir Experiments/HCV-1b/PB/100/90  --output_dir  Experiments/HCV-1b/PB/100/90

python plots_scripts/HCV-1b/100/PB/90/plot_relative_abundance_error.py --input_dir Experiments/HCV-1b/PB/100/90 --gt_dir /tudelft.net/staff-umbrella/ViralQuasispecies/inika/Read_simulators/Experiments_data_planned/HCV-1b/PacBio-Hifi/100x/90  --output_dir  Experiments/HCV-1b/PB/100/90

python plots_scripts/HCV-1b/100/PB/90/plot_community_evaluation_results.py --input_dir Experiments/HCV-1b/PB/100/90 --output_dir  Experiments/HCV-1b/PB/100/90

python plots_scripts/HCV-1b/100/PB/90/make_summary_statistics.py --input_dir Experiments/HCV-1b/PB/100/90 --gt_dir /tudelft.net/staff-umbrella/ViralQuasispecies/inika/Read_simulators/Experiments_data_planned/HCV-1b/PacBio-Hifi/100x/90  --output_dir  Experiments/HCV-1b/PB/100/90

python plots_scripts/make_summary_statistics.py --sequence_ids HCV-1b_90 --input_dir Experiments/HCV-1b/PB/100/90 --gt_input_dir /tudelft.net/staff-umbrella/ViralQuasispecies/inika/Read_simulators/Experiments_data_planned/HCV-1b/PacBio-Hifi/100x/90  --output_dir  Experiments/HCV-1b/PB/100/90