

python plots_scripts/HIV-1/PB/1000/plot_abundance_evaluations.py --input_dir Experiments/HIV-1/PacBio-Hifi/1000 --gt_dir /tudelft.net/staff-umbrella/ViralQuasispecies/inika/Read_simulators/Experiments_data_planned/HIV-1/PacBio-Hifi/99/1000 --output_dir  Experiments/HIV-1/PacBio-Hifi/1000

python plots_scripts/HIV-1/PB/1000/plot_community_evaluation_results.py --input_dir Experiments/HIV-1/PacBio-Hifi/1000 --output_dir  Experiments/HIV-1/PacBio-Hifi/1000

python plots_scripts/HIV-1/PB/1000/plot_edit_distance_to_consensus.py --input_dir Experiments/HIV-1/PacBio-Hifi/1000  --output_dir  Experiments/HIV-1/PacBio-Hifi/1000

python plots_scripts/HIV-1/PB/1000/plot_relative_abundance_error.py --input_dir Experiments/HIV-1/PacBio-Hifi/1000 --gt_dir /tudelft.net/staff-umbrella/ViralQuasispecies/inika/Read_simulators/Experiments_data_planned/HIV-1/PacBio-Hifi/99/1000 --output_dir  Experiments/HIV-1/PacBio-Hifi/1000

python plots_scripts/HIV-1/PB/1000/make_summary_statistics.py --input_dir Experiments/HIV-1/PacBio-Hifi/1000 --gt_input_dir /tudelft.net/staff-umbrella/ViralQuasispecies/inika/Read_simulators/Experiments_data_planned/HIV-1/PacBio-Hifi/99/1000 --output_dir  Experiments/HIV-1/PacBio-Hifi/1000

python plots_scripts/make_summary_statistics.py --sequence_ids HIV-1 --input_dir Experiments/HIV-1/PacBio-Hifi/1000 --gt_input_dir /tudelft.net/staff-umbrella/ViralQuasispecies/inika/Read_simulators/Experiments_data_planned/HIV-1/PacBio-Hifi/99/1000 --output_dir  Experiments/HIV-1/PacBio-Hifi/1000