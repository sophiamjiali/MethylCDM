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