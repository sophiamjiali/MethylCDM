## Data
- data held in /ddn_exa/campbell/datasets/sophia2025
- just download and preprocess locally, move AnnData only


(To run a script)
conda activate methylcdm-env
pip install -e .
python ./scripts/process_methylation.py \
  --project TCGA-UVM \
  --config_pipeline pipeline.yaml \
  --config_preproc methylation_preproc.yaml




Final shape: 

(below fetched from RNA-CDM paper)

Project Code | Cancer Type                                            | # of Samples | # Downloaded | Post-Processing
---------------------------------------------------------------------------------------------------------------------
TCGA-LUAD    | Lung adenocarcinoma                                    | 520          |
TCGA-KIRP    | Kidney renal papillary cell carcinoma                  | 298          |
TCGA-COAD    | Colon adenocarcinoma                                   | 289          | 555          | 347
TCGA-CESC    | Cervical squamous cell carcinoma & endocervical adeno. | 277          | 312          | 308
TCGA-GBM     | Glioblastoma multiforme                                | 212          | 450          | 153
TCGA-PAAD    | Pancreatic adenocarcinoma                              | 202          | 195          | 189
TCGA-ESCA    | Oesophageal carcinoma                                  | 156          | 202          | 201
TCGA-OV      | Ovarian serous cystadenocarcinoma                      | 83           | 623          | 10
TCGA-UVM     | Uveal melanoma                                         | 80           | 80           | 78
TCGA-CHOL    | Cholangiocarcinoma                                     | 36           | 45           | 45


Downloaded Status
[] TCGA-LUAD
[] TCGA-KIRP (currently downloading)
[X] TCGA-COAD
[X] TCGA-CESC
[X] TCGA-GBM
[X] TCGA-PAAD
[X] TCGA-ESCA
[X] TCGA-OV 
[X] TCGA-UVM
[X] TCGA-CHOL
