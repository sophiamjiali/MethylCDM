## Data
- data held in /ddn_exa/campbell/datasets/sophia2025
- just download and preprocess locally, move AnnData only


(To run a script)
conda activate methylcdm-env
pip install -e .
python ./scripts/process_methylation.py \
  --project TCGA-CHOL \
  --config_pipeline ./config/pipeline.yaml \
  --config_preproc ./config/methylation_preproc.yaml




(below fetched from RNA-CDM paper)

Project Code | Cancer Type                                            | Number of Samples
-------------------------------------------------------------------------------------------
TCGA-LUAD    | Lung adenocarcinoma                                    | 520
TCGA-KIRP    | Kidney renal papillary cell carcinoma                  | 298
TCGA-COAD    | Colon adenocarcinoma                                   | 289
TCGA-CESC    | Cervical squamous cell carcinoma & endocervical adeno. | 277
TCGA-GBM     | Glioblastoma multiforme                                | 212
TCGA-PAAD    | Pancreatic adenocarcinoma                              | 202
TCGA-ESCA    | Oesophageal carcinoma                                  | 156
TCGA-OV      | Ovarian serous cystadenocarcinoma                      | 83
TCGA-UVM     | Uveal melanoma                                         | 80
TCGA-CHOL    | Cholangiocarcinoma                                     | 36


Downloaded Status
[] TCGA-LUAD
[] TCGA-KIRP
[] TCGA-COAD
[] TCGA-CESC
[] TCGA-GBM
[] TCGA-PAAD
[] TCGA-ESCA
[] TCGA-OV 
[] TCGA-UVM
[] TCGA-CHOL
