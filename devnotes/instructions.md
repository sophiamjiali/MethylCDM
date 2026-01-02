# Instructions

Instructions for an oblivious developer (me)

## 1. Download and process DNA methylation data

```
conda activate methylcdm-env
pip install -e .
python ./scripts/process_methylation.py \
  --project {TCGA-###} \
  --config_pipeline pipeline.yaml \
  --config_preproc methylation_preproc.yaml \
  --verbose True
```

## 2. Prepare DNA methylation data for model training
```
conda activate methylcdm-env
pip install -e .
python ./scripts/prepare_data.py \
  --config_pipeline pipeline.yaml \
  --config_betaVAE betaVAE.yaml \
  --verbose True
```

## 3. Train Beta-VAE by initializing a hyperparameter sweep

If running locally:
```
conda activate methylcdm-env
pip install -e .
python ./scripts/train_betaVAE.py \
  --config_pipeline pipeline.yaml \
  --config_train betaVAE.yaml \
  --verbose True
```

If importing to HPC for the first time:
```
# Import files to HPC
cd /Volumes/FBI_Drive/MethylCDM-project
scp -r . sophiali@dev1gpu.mshri.on.ca:/ddn_exa/campbell/sli/methylcdm-project

# Log into HPC
ssh sophiali@galen.mshri.on.ca
cd /ddn_exa/campbell/sli/methylcdm-project

# Download Miniforge (Conda) if necessary
curl -L -O https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash Miniforge3-Linux-x86_64.sh
source ~/miniforge3/bin/activate              

# Create the virtual environment
conda env create -f environment.yaml
```

To initiate a hyperparameter sweep on the HPC:
```
# Install the package in the environment
ssh sophiali@gpu1.galen.mshri.on.ca
cd /ddn_exa/campbell/sli/methylcdm-project
source ~/miniforge3/bin/activate 
conda env create -f environment.yaml
pip install -e .

# Queue the sweep to the SLURM scheduler
sbatch slurm/betaVAE_sweep.sbatch

# Monitor job progress
squeue -u $USER
```