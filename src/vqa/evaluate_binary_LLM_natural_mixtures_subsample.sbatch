#!/bin/sh
#SBATCH --partition=general,insy # Request partition. Default is 'general' 
#SBATCH --qos=medium         # Request Quality of Service. Default is 'short' (maximum run time: 4 hours)
#SBATCH --time=15:00:00      # Request run time (wall-clock). Default is 1 minute
#SBATCH --ntasks=1         # Request number of parallel tasks per job. Default is 1
#SBATCH --cpus-per-task=1   # Request number of CPUs (threads) per task. Default is 1 (note: CPUs are always allocated to jobs per 2).
#SBATCH --mem-per-cpu=20GB          # Request memory (MB) per node. Default is 1024MB (1GB). For multiple tasks, specify --mem-per-cpu instead
#SBATCH --output=slurm_nm_sub_%j.out # Set name of output log. %j is the Slurm jobId
#SBATCH --error=slurm_nm_sub_%j.err # Set name of error log. %j is the Slurm jobId
#SBATCH --gres=gpu:a40:1 # Request 1 GPU

/usr/bin/nvidia-smi # Check sbatch settings are working (it should show the GPU that you requested)

# Measure GPU usage of your job (initialization)
previous=$(/usr/bin/nvidia-smi --query-accounted-apps='gpu_utilization,mem_utilization,max_memory_usage,time' --format='csv' | /usr/bin/tail -n '+2') 
  
module load cuda/12.2

# for mix in 07_98
# do  
#     echo "srun python evaluate_binary_LLM_command_line_nm.py --devices 1 --virus_name SARS-CoV-2 --skip_grs 21 --path_to_dataset /tudelft.net/staff-umbrella/ViralQuasispecies/inika/Read_simulators/data/lumc_data/natural_mixtures/$mix/subsampled_reads_1000 --outdir Experiments/lumc_subsample/1000/$mix"
#     srun python evaluate_binary_LLM_command_line_nm.py --devices 1 --virus_name SARS-CoV-2 --skip_grs 26 --path_to_dataset /tudelft.net/staff-umbrella/ViralQuasispecies/inika/Read_simulators/data/lumc_data/natural_mixtures/$mix/subsampled_reads_1000 --outdir Experiments/lumc_subsample/1000/$mix
# done

# Record the start time
start_time=$(date +%s)

for mix in 03_50
do  
    for seed in 100
    do
        mkdir Experiments/lumc_subsample/$mix/gr_25712/seed_$seed  
        srun python evaluate_binary_LLM_command_line_nm.py --devices 1 --virus_name SARS-CoV-2 --path_to_dataset /tudelft.net/staff-umbrella/ViralQuasispecies/inika/Read_simulators/data/lumc_data/natural_mixtures/$mix/seed_$seed/subsampled_reads_1000_25712 --outdir Experiments/lumc_subsample/$mix/gr_25712/seed_$seed  
    done
done

# Record the end time
end_time=$(date +%s)

# Calculate the duration
duration=$((end_time - start_time))

# Display the duration
echo "Duration: $duration seconds"

# Measure GPU usage of your job (result)
/usr/bin/nvidia-smi --query-accounted-apps='gpu_utilization,mem_utilization,max_memory_usage,time' --format='csv' | /usr/bin/grep -v -F "$previous"