#!/bin/bash
#SBATCH --job-name=betaVAE_test
#SBATCH --output=/ddn_exa/campbell/sli/methylcdm-project/logs/betaVAE_test_%j.out
#SBATCH --error=/ddn_exa/campbell/sli/methylcdm-project/logs/betaVAE_test_%j.err
#SBATCH --time=24:30:00
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --nodelist=gpu3

echo "=========================================="
echo "Job ID:       $SLURM_JOB_ID"
echo "Node:         $SLURMD_NODENAME"
echo "GPUs:         $CUDA_VISIBLE_DEVICES"
echo "Start time:   $(date)"
echo "=========================================="

# Activate the project virtual environment
source $HOME/envs/methylcdm/bin/activate

# Navigate to the project directory
cd /ddn_exa/campbell/sli/methylcdm-project

srun python slurm/betaVAE_test.py \
    --config configs/betaVAE.yaml \
    --run_name "test_baseline"
