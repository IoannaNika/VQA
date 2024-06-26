#!/bin/sh
#SBATCH --partition=general,insy # Request partition. Default is 'general' 
#SBATCH --qos=long         # Request Quality of Service. Default is 'short' (maximum run time: 4 hours)
#SBATCH --time=168:00:00      # Request run time (wall-clock). Default is 1 minute
#SBATCH --ntasks=1         # Request number of parallel tasks per job. Default is 1
#SBATCH --cpus-per-task=1   # Request number of CPUs (threads) per task. Default is 1 (note: CPUs are always allocated to jobs per 2).
#SBATCH --mem-per-cpu=5120          # Request memory (MB) per node. Default is 1024MB (1GB). For multiple tasks, specify --mem-per-cpu instead
#SBATCH --output=slurm_deletions_%j.out # Set name of output log. %j is the Slurm jobId
#SBATCH --error=slurm_deletions_%j.err # Set name of error log. %j is the Slurm jobId
#SBATCH --gres=gpu:1# Request 1 GPU

/usr/bin/nvidia-smi # Check sbatch settings are working (it should show the GPU that you requested)

# Measure GPU usage of your job (initialization)
previous=$(/usr/bin/nvidia-smi --query-accounted-apps='gpu_utilization,mem_utilization,max_memory_usage,time' --format='csv' | /usr/bin/tail -n '+2') 

module load cuda/12.2

for n in 60 70 80 90 100
do 
    mkdir Experiments/deletions/$n
    echo "srun python evaluate_binary_LLM_command_line_indelm.py --devices 1 --virus_name SARS-CoV-2 --path_to_dataset /tudelft.net/staff-umbrella/ViralQuasispecies/inika/Read_simulators/Experiments_data_planned/deletion_datasets/$n/dataset --outdir Experiments/deletions/$n"
    srun python evaluate_binary_LLM_command_line_indelm.py --devices 1 --virus_name SARS-CoV-2 --path_to_dataset /tudelft.net/staff-umbrella/ViralQuasispecies/inika/Read_simulators/Experiments_data_planned/deletion_datasets/$n/dataset --outdir Experiments/deletions/$n
done


# Measure GPU usage of your job (result)
/usr/bin/nvidia-smi --query-accounted-apps='gpu_utilization,mem_utilization,max_memory_usage,time' --format='csv' | /usr/bin/grep -v -F "$previous"