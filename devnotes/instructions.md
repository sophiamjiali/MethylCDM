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