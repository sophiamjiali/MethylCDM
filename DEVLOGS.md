# Developer Logs

## DNA Methylation Workflow
- pandas, numpy, parquet, anndata, pycombat (combat batch correction)
1. raw per-file TCGA downloads
2. load all files
3. harmonize sample barcodes
4. concatenate into a CpG x sample matrix, save
5. Perform probe QC / sample QC on the full matrix
6. Perform beta value imputation
7. Normalize beta values (optional)
8. Save final cohort-level matrix


### 11/30/2025 - Log 6:

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