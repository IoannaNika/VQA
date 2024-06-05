# python write_consensus_to_fasta.py --dir Experiments/sars_cov_2/PacBio-Hifi/100/20
# python plots_scripts/sars-cov-2/quast.py --references_file /tudelft.net/staff-umbrella/ViralQuasispecies/inika/Read_simulators/Experiments_data_planned/sars_cov_2/data/20/sequences.fasta --data_dir Experiments/sars_cov_2/PacBio-Hifi/100/20 


for haplotype in 20 15 10 5 1
do
    python write_consensus_to_fasta.py --dir Experiments/sars_cov_2/PacBio-Hifi/100/$haplotype
    python plots_scripts/sars-cov-2/evaluation_v2.py --template_dir /tudelft.net/staff-umbrella/ViralQuasispecies/inika/Read_simulators/Experiments_data_planned/sars_cov_2/data/$haplotype/NCBI_processed --data_dir /tudelft.net/staff-umbrella/ViralQuasispecies/inika/VQA/src/vqa/Experiments/sars_cov_2/PacBio-Hifi/100/$haplotype/ --n_true_haplotypes $haplotype
done

python plots_scripts/sars-cov-2/plot_evaluation_v2.py --directory Experiments/sars_cov_2/PacBio-Hifi/100


# python plots_scripts/sars-cov-2/plot_edit_distances.py --data_dir /tudelft.net/staff-umbrella/ViralQuasispecies/inika/Read_simulators/Experiments_data_planned/sars_cov_2/data/  --outdir Experiments/sars_cov_2/PacBio-Hifi/100

python plots_scripts/sars-cov-2/make_summary_statistics.py --data_dir /tudelft.net/staff-umbrella/ViralQuasispecies/inika/Read_simulators/Experiments_data_planned/sars_cov_2/data/  --outdir Experiments/sars_cov_2/PacBio-Hifi/100 --directory Experiments/sars_cov_2/PacBio-Hifi/100

python plots_scripts/make_summary_statistics.py --sequence_ids sars-cov-2 --input_dir Experiments/sars_cov_2/ONT --gt_input_dir /tudelft.net/staff-umbrella/ViralQuasispecies/inika/Read_simulators/Experiments_data_planned/sars_cov_2/ONT --output_dir Experiments/sars_cov_2/ONT
