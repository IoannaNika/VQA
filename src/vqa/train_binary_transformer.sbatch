#!/bin/sh
#SBATCH --partition=general,insy # Request partition. Default is 'general' 
#SBATCH --qos=long         # Request Quality of Service. Default is 'short' (maximum run time: 4 hours)
#SBATCH --time=168:00:00      # Request run time (wall-clock). Default is 1 minute
#SBATCH --ntasks=1          # Request number of parallel tasks per job. Default is 1
#SBATCH --cpus-per-task=4   # Request number of CPUs (threads) per task. Default is 1 (note: CPUs are always allocated to jobs per 2).
#SBATCH --mem=20GB          # Request memory (MB) per node. Default is 1024MB (1GB). For multiple tasks, specify --mem-per-cpu instead
#SBATCH --output=train_slurm_%j.out # Set name of output log. %j is the Slurm jobId
#SBATCH --error=train_slurm_%j.err # Set name of error log. %j is the Slurm jobId
#SBATCH --gres=gpu:a40:2 # Request 1 GPU

/usr/bin/nvidia-smi # Check sbatch settings are working (it should show the GPU that you requested)

# Measure GPU usage of your job (initialization)
previous=$(/usr/bin/nvidia-smi --query-accounted-apps='gpu_utilization,mem_utilization,max_memory_usage,time' --format='csv' | /usr/bin/tail -n '+2') 

module load cuda/12.2

# echo "srun python train_binary_transformer.py --devices 2 --checkpoint_dir final_checkpoints/PB --pb True"
# srun python train_binary_transformer.py --devices 2 --checkpoint_dir final_checkpoints/PB --pb True

# echo "srun python train_binary_transformer.py --devices 2 --checkpoint_dir final_checkpoints/ONT --pb False"
# srun python train_binary_transformer.py --devices 2 --checkpoint_dir final_checkpoints/ONT --pb False

echo "srun python train_binary_transformer.py --devices 2 --checkpoint_dir final_checkpoints/PB_random_lr_2_acc10_r2 --pb True --init_random True"
srun python train_binary_transformer.py --devices 2 --checkpoint_dir final_checkpoints/PB_random_lr_2_acc10_r2 --pb True --init_random True

# Measure GPU usage of your job (result)
/usr/bin/nvidia-smi --query-accounted-apps='gpu_utilization,mem_utilization,max_memory_usage,time' --format='csv' | /usr/bin/grep -v -F "$previous"