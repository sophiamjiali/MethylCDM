# MethylCDM
Generative Modeling of Tumour Histology from Differential Methylation

# Repository Structure

```
MethylCDM/
├── config/             # Pipeline Configuration Files
├── data/               # Input & Output Data
├── experiments/        # Training Logs, Metrics, & Figures
├── models/             # Model Checkpoints
├── notebooks/          # Exploratory Notebooks for Development
├── scripts/            # Pipeline Entry-points
├── src/                # Core Reusable Modules
├── tools/              # Third-party tools & APIs
└── workflow/           # SnakeMake Pipeline Orchestration
```

## Scripts
The `scripts` folder holds the entry-points to the pipeline, using the logic defined in the `src` folder. Shell wrappers exist to run the pipeline for cluster/HPC.
```
└── script/
    ├── preprocess_methylation.py   # Preprocess methylation data
    ├── preprocess_wsi.py           # Preprocess WSI data
    ├── train_betaVAE.py            # Train β-VAE (start Optuna sweep)
    ├── train_diffusion.py          # Train CDM
    ├── eval_betaVAE.py             # Evaluate β-VAE
    └── eval_diffusion.py           # Evaluate CDM
```

## SRC
The `src` folder holds the core logic and main stateless implementation of the workflow. The functions defined are used in the `scripts` folder to run the full pipeline.
```
└── src/
    ├── data/                           # Downloading & Datamodule Logic
    │   ├── load_methylation.py         # - Load & download methylation data
    │   └── methylation_datamodule.py   # - Methylation DataModule & Dataset classes
    ├── preprocessing/                  # Methylation & WSI Preprocessing Logic
    │   ├── process_methylation.py      # - Methylation quality control & preprocessing
    │   └── reconcile_methylation.py    # - Methylation cohort aggregation
    ├── models/                         # Model Definitions
    │   └── betaVAE.py                  # - β-VAE LightningModule class
    ├── training/                       # Training Logic
    │   └── betaVAE_objective.py        # - β-VAE Optuna objective function
    ├── evaluation/                     # Evaluation Logic
    │   └── ...
    ├── utils/                          # General Utilities
    │   ├── training_utils.py           # - Callback & Logging utilities
    │   └── utils.py                    # - General utilities
    └── constants.py                    # Global constants

```
