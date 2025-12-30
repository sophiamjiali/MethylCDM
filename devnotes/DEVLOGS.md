# Developer Logs

## Current TO-DO
- consolidate anndatas together: filter for the common genes, append metadata
- design beta architecture
- begin training



### RNA-CDM Beta-VAE notes
- considered 12 cancer types for training; five used during training and validation
    - LUAD, GBM, KIRP, CESC, COAD
- NaN values were removed
- selected common genes between all cancer types (17,655 genes)
- gene expression log-transformed and normalized with z-score from training set values

- final architecture: empirically determined
    - two hidden layers of 6k and 4k neurons each for both the encoder and decoder
    - size of 200 for the latent dimension
    - batch norm btween the layers
    - LeakyReLU as the activation function
    - beta value of 0.005 was used in the loss function
    - Adam optimizer for training, lr = 3x10-3
    - warm-up and cosine learning rate scheduler
    - mean square error as loss function
    - trained for 250 epochs
    - early stopping based on validation set loss
    - batch size of 128
- 60-20-20% training, validation, test stratified splits

#### Considerations
- RNA-seq modeled after log/CPM transforms with continuous, heavy-tailed distributions
- beta values are bounded to [0, 1] and usually bimodal; may affect choice of reconstruction loss
    - e.g. Gaussian versus beta/bernoulli-like
    - output activation to sigmoid?



### 12/30/2025 - Log 13:
- split manifest into temporary chunks to avoid API crashes; retrying donwloading

### 12/29/2025 - Log 12:
- successfully ran downloading and preprocessing workflow for all twelve projects
    - deleted raw files to save space, seed is set for workflow so preprocessing is reproducible 
- next, reconcile gene-level matrix into one AnnData object; stratify into training-testing splits
- save stratified mass-object, import to HPC
- design architecture, aligning with Beta-VAE notes (align with RNA-CDM first)

### 12/28/2025 - Log 11:
- changed to only preserve Illumina 450K assay, removed 27K assay
- as first pass, preserve the richest coverage; harmonization would require keeping only the overlapping probes
- finished full preprocesing pipeline; outputted and saved processed AnnData object
- starting with the RNA-CDM beta-VAE architecture to start
- listed data to download in `data_info.md`

### 12/25/2025 - Log 10: My chungus life, working on Christmas Day
- removed non-standard probes during probe quality control preprocessing
- refactored beta value parquet loading to parallelization to fix bottleneck
- moved sample QC to after probe QC to avoid inflating missing values
- need to refactor preprocessing to subset to only the common CpG probes

### 12/23/2025 - Log 9:
- consolidated cross-reactive, SNP, and Illumina annotations in `resources/`
- bugfixed `process_array_methylation()`, need to do a final runthrough of the for-loop in `process_methylation()`
- need to bugfix anndata initializing, then move onto beta-VAE architecture and training logic

### 12/22/2025 - Log 8:
- refactored the preprocessing pipeline to separate by array type
- adjusted preprocessing pipeline to save processed matrix as a gene-level matrix

### 12/16/2025 - Log 7:
- finished bugfixing `download_methylation()`, `clean_methylation_data()`, and `merge_cohort()` for downloading and cleaning logic
- next, move onto `process_methylation()`

### 12/06/2025 - Log 6: 
- created resource folder for package-provided data; downloaded and added DNA methylation annotation manifests to the package
- filtered methylation for common probes across manifests, then used larger manifest (i.e. 450 versus 27) as a proxy to get annotations if missing (i.e. SNP)
- finished implementing preprocessing pipeline, needs debugging
- added cleaning step to parse downloaded data and convert to just parquets

### 11/25/2025 - Log 5:
- downloaded TCGA-BRCA methylation data (~12 GB)
- removed metadata fields from configurations, standardized to `src/MethylCDM/constants.py`
    - too complicated to normalize the nested fields, better for downstream analysis
- user will have to provide SNP, cross-reactive, and sex chromosome probe annotations
- built structure for preprocessing workflow; need to download probe annotation manifests

### 11/22/2025 - Log 4:  gdc-client blew up my laptop
- renamed `preprocess_*` to just `process_*`
- moved `src/*` to `src/MethylCDM` for better align to importing standards 
- subsequently renamed project root to `MethylCDM-project` to prevent namespace shadowingW
- changed all pathing to support both absolute and relative; added `src/MethylCDM/constants.py`
- specified that users should download the GDC Data Transfer Tool and place it in `tools/`
    - need to update README to provide instructions for this

### 11/18/2025 - Log 3: MethylCDM now exists in a USB called 'FBI Drive'
- refactored DNA methylation downloading, loading, preprocessing, and entry-point script workflow to be performed per project for better integration with SnakeMake
- will either need to have an alternate workflow/file that aggregates all data together and have it exist under a 'project name' like 'pan-cancer' so the preprocessing (?) and gene clustering analysis can work on both a single project and the pan-cancer approach when training the beta-VAE
    - propagate this to the later workflow (CDM training) to facilitate pan-cancer approach eventually
- next: complete `load_raw_methylation()` and bugtest downloading logic

### 11/14/2025 - Log 2: Downloading all files and blowing up the H4H cluster
- Began defining GDC API query logic in `src/data/load_methylation.py`
- began populating the `methylation_preproc.yaml` configurations file

### 11/13/2025 - Log 1: Too many config files did not kill the pope
- Initialized repository and structure