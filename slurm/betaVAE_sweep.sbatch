#!/bin/bash
#SBATCH --job-name=betaVAE_sweep
#SBATCH --output=/ddn_exa/campbell/sli/methylcdm-project/logs/betaVAE_sweep_%j.out
#SBATCH --error=/ddn_exa/campbell/sli/methylcdm-project/logs/betaVAE_sweep_%j.err
#SBATCH --time=72:00:00
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --nodelist=gpu1

# Activate Conda Environment
source ~/miniforge3/etc/profile.d/conda.sh
conda activate methylcdm-env

# Navigate to the project directory
cd /ddn_exa/campbell/sli/methylcdm-project

# Install the project as a package to resolve imports
pip install -e .

# Run the sweep
python scripts/train_betaVAE.py \
    --config_pipeline pipeline.yaml \
    --config_train betaVAE.yaml \
    --verbose True