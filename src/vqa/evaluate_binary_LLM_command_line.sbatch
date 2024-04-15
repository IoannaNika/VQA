#!/bin/sh
#SBATCH --partition=general,insy # Request partition. Default is 'general' 
#SBATCH --qos=long         # Request Quality of Service. Default is 'short' (maximum run time: 4 hours)
#SBATCH --time=168:00:00      # Request run time (wall-clock). Default is 1 minute
#SBATCH --ntasks=1         # Request number of parallel tasks per job. Default is 1
#SBATCH --cpus-per-task=1   # Request number of CPUs (threads) per task. Default is 1 (note: CPUs are always allocated to jobs per 2).
#SBATCH --mem-per-cpu=20GB          # Request memory (MB) per node. Default is 1024MB (1GB). For multiple tasks, specify --mem-per-cpu instead
#SBATCH --output=slurm_all%j.out # Set name of output log. %j is the Slurm jobId
#SBATCH --error=slurm_all%j.err # Set name of error log. %j is the Slurm jobId
#SBATCH --gres=gpu:3 # Request 1 GPU

/usr/bin/nvidia-smi # Check sbatch settings are working (it should show the GPU that you requested)

# Measure GPU usage of your job (initialization)
previous=$(/usr/bin/nvidia-smi --query-accounted-apps='gpu_utilization,mem_utilization,max_memory_usage,time' --format='csv' | /usr/bin/tail -n '+2') 

module load cuda/12.2

# # HIV
# echo "srun python evaluate_binary_LLM_command_line.py  --virus_name HIV-1 --path_to_dataset /tudelft.net/staff-umbrella/ViralQuasispecies/inika/Read_simulators/Experiments_data/HIV-1/99/cov_100x_mixture_pb_hifi_094003003/dataset --outdir Experiments/HIV-1/99/cov_100x_mixture_pb_hifi_094003003/"
# srun python evaluate_binary_LLM_command_line.py  --virus_name HIV-1 --path_to_dataset /tudelft.net/staff-umbrella/ViralQuasispecies/inika/Read_simulators/Experiments_data/HIV-1/99/cov_100x_mixture_pb_hifi_094003003/dataset --outdir Experiments/HIV-1/99/cov_100x_mixture_pb_hifi_094003003/

# echo "srun python evaluate_binary_LLM_command_line.py  --virus_name HIV-1 --path_to_dataset /tudelft.net/staff-umbrella/ViralQuasispecies/inika/Read_simulators/Experiments_data/HIV-1/99/cov_100x_mixture_pb_hifi_040303_reverse_complement/dataset --outdir Experiments/HIV-1/99/cov_100x_mixture_pb_hifi_040303_reverse_complement/"
# srun python evaluate_binary_LLM_command_line.py  --virus_name HIV-1 --path_to_dataset /tudelft.net/staff-umbrella/ViralQuasispecies/inika/Read_simulators/Experiments_data/HIV-1/99/cov_100x_mixture_pb_hifi_040303_reverse_complement/dataset --outdir Experiments/HIV-1/99/cov_100x_mixture_pb_hifi_040303_reverse_complement/


# cerbaresearch 
# echo "srun python evaluate_binary_LLM_command_line.py --path_to_dataset /tudelft.net/staff-umbrella/ViralQuasispecies/inika/Read_simulators/Experiments_data/cerbaresearch/Pacbio-nsp12-fragment/1005/dataset --outdir /tudelft.net/staff-umbrella/ViralQuasispecies/inika/VQA/src/vqa/Experiments/cerbaresearch/Pacbio-nsp12-fragment/1005/ --virus_name SARS-CoV-2"
# srun python evaluate_binary_LLM_command_line.py --path_to_dataset /tudelft.net/staff-umbrella/ViralQuasispecies/inika/Read_simulators/Experiments_data/cerbaresearch/Pacbio-nsp12-fragment/1005/dataset --outdir /tudelft.net/staff-umbrella/ViralQuasispecies/inika/VQA/src/vqa/Experiments/cerbaresearch/Pacbio-nsp12-fragment/1005/ --virus_name SARS-CoV-2

# echo "srun python evaluate_binary_LLM_command_line.py --path_to_dataset /tudelft.net/staff-umbrella/ViralQuasispecies/inika/Read_simulators/Experiments_data/cerbaresearch/Pacbio-nsp12-fragment/1007/dataset --outdir /tudelft.net/staff-umbrella/ViralQuasispecies/inika/VQA/src/vqa/Experiments/cerbaresearch/Pacbio-nsp12-fragment/1007/ --virus_name SARS-CoV-2"
# srun python evaluate_binary_LLM_command_line.py --path_to_dataset /tudelft.net/staff-umbrella/ViralQuasispecies/inika/Read_simulators/Experiments_data/cerbaresearch/Pacbio-nsp12-fragment/1007/dataset --outdir /tudelft.net/staff-umbrella/ViralQuasispecies/inika/VQA/src/vqa/Experiments/cerbaresearch/Pacbio-nsp12-fragment/1007/ --virus_name SARS-CoV-2

# LUMC 
# echo "srun python evaluate_binary_LLM_command_line.py --path_to_dataset /tudelft.net/staff-umbrella/ViralQuasispecies/inika/Read_simulators/data/lumc_rev_compl --outdir /tudelft.net/staff-umbrella/ViralQuasispecies/inika/VQA/src/vqa/Experiments/lumc/rev_compl --virus_name SARS-CoV-2"
# srun python evaluate_binary_LLM_command_line.py --path_to_dataset /tudelft.net/staff-umbrella/ViralQuasispecies/inika/Read_simulators/data/lumc_rev_compl --outdir /tudelft.net/staff-umbrella/ViralQuasispecies/inika/VQA/src/vqa/Experiments/lumc/test --virus_name SARS-CoV-2

echo "srun python evaluate_binary_LLM_command_line.py --path_to_dataset /tudelft.net/staff-umbrella/ViralQuasispecies/inika/Read_simulators/Experiments_data/lumc/dataset --outdir Experiments/lumc/mixtures/artificial_all --virus_name SARS-CoV-2 --lumc False --devices 3"
srun python evaluate_binary_LLM_command_line.py --path_to_dataset /tudelft.net/staff-umbrella/ViralQuasispecies/inika/Read_simulators/Experiments_data/lumc/dataset --outdir Experiments/lumc/mixtures/artificial_all --virus_name SARS-CoV-2 --lumc False --devices 3


# Measure GPU usage of your job (result)
/usr/bin/nvidia-smi --query-accounted-apps='gpu_utilization,mem_utilization,max_memory_usage,time' --format='csv' | /usr/bin/grep -v -F "$previous"