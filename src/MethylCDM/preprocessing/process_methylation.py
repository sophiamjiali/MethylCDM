# ==============================================================================
# Script:           process_methylation.py
# Purpose:          Performs quality control and preprocesses methylation data
# Author:           Sophia Li
# Affiliation:      CCG Lab, Princess Margaret Cancer Center, UHN, UofT
# Date:             11/25/2025
# ==============================================================================

from MethylCDM.utils.utils import load_annotations

# =====| Preprocessing Wrapper |================================================

def process_methylation(cpg_matrix, metadata, config):
    """
    Preprocesses and returns DNA methylation data in the form of a CpG x Sample
    ID matrix. Performs sample-level and probe-level quality control (QC) and
    imputation using k-nearest neighbours (KNN) if toggled in the associated 
    configuration file.

    The CpG matrix provided was already filtered for the common set of probes 
    present across all datapoints to control for if multiple manifests were
    used (i.e. Illumina 450K, 27K, EPIC). The largest of the manifests present
    in the CpG matrix will be used as the canonical annotations for the common
    set of probes (EPIC > 450K > 27K).

    The following preprocessing is applied if toggled:
        - Sample-level QC:  remove samples with > 5% CpGs missing
        - Probe-level QC:   remove CpGs with > 5% values missing, SNP-affected, 
                            cross-reactive, and sex chromosome probes
        - Imputation:       perform probe-wise mean imputation on missing values

    Parameters
    ----------
    cpg_matrix (DataFrame): a CpG x Sample ID matrix of beta values
    metadata (DataFrame): a metadata matrix corresponding to the beta values
    config (dict): a configuration object containing:
        - toggle_sample_filtering (boolean): perform sample filtering
        - toggle_probe_filtering (boolean): perform probe filtering
        - togge_imputation (boolean): perform imputation

    Returns
    ------
    (DataFrame): processed cohort-level CpG x Sample ID matrix of 
                            beta values of a given project

    Raises
    ------
    FileNotFoundError: if the probe annotation files do not exist or is not a 
                       file
    """

    # Fetch the toggles and preprocessing thresholds from the configurations
    preproc_cfg = config.get('preprocess', {})
    max_missing_sample = preproc_cfg.get('max_missing_sample', 0)
    max_missing_probe = preproc_cfg.get('max_missing_probe', 0)
    n_neighbours = preproc_cfg.get('n_neighbours', 0)

    # Fetch the dominant annotation based on manifests present
    manifests_present = set(metadata['platform'])
    annotation = load_annotations(manifests_present)

    # 1. Sample quality control
    if config.get('toggle_sample_filtering', False):
        cpg_matrix = sample_qc(cpg_matrix, metadata, max_missing_sample)
    
    # 2. Probe quality control
    if config.get('toggle_probe_filtering', False):
        cpg_matrix = probe_qc(cpg_matrix, annotation, max_missing_probe)

    # 3. Missing value imputation
    if config.get('toggle_imputation', False):
        cpg_matrix = impute_missing(cpg_matrix, n_neighbours)

    return cpg_matrix


# =====| Preprocessing Helpers |================================================

def sample_qc(cpg_matrix, max_missing):
    """
    Removes samples with too many missing values.

    Parameters
    ----------
    cpg_matrix (DataFrame): a CpG x Sample ID matrix of beta values
    metadata (DataFrame): a metadata matrix corresponding to the beta values
    max_missing (float): maximum allowed fraction of missing CpGs per sample

    Returns
    -------
    (DataFrame): filtered matrix with samples below the threshold removed
    
    """
    sample_missing = cpg_matrix.isnull().mean(axis = 0)
    keep_samples = sample_missing <= max_missing
    return cpg_matrix.loc[:, keep_samples]


def probe_qc(cpg_matrix, annotation, max_missing):
    """
    Removes low-quality probes with too many missing values.

    Parameters
    ----------
    cpg_matrix (DataFrame): a CpG x Sample ID matrix of beta values
    annotations (DataFrame): probe annotations
    max_missing (float): maximum allowed fraction of missing samples per probe

    Returns
    -------
    cpg_matrix (DataFrame): the cleaned CpG matrix 
    """
    
    # Define a filter to remove probes with high missingness
    missing_frac = cpg_matrix.isna().mean(axis = 1)
    keep_missing = missing_frac <= max_missing

    # Remove non-variable probes
    keep_variable = cpg_matrix.var(axi = 1) > 0

    # Remove SNP-affected, cross-reactive, and sex chromosome probes
    keep_snp = annotation["is_snp"] == False
    keep_cross = annotation["is_cross_reactive"] == False
    keep_sex = ~annotation["chr"].isin(["X", "Y", "chrX", "chrY"])

    # Combine all filters and subset the matrix
    keep = keep_missing & keep_variable & keep_snp & keep_cross & keep_sex
    return cpg_matrix.loc[keep]


def impute_missing(cpg_matrix):
    """
    Imputes missing beta values using probe-wise mean imputation as 
    recommende dfor TCGA-style beta value matrices.

    Parameters
    ----------
    cpg_matrix (DataFrame): a CpG x Sample ID matrix of beta values

    Returns
    -------
    cpg_matrix (DataFrame): the imputed CpG matrix
    """
    return cpg_matrix.apply(lambda row: row.fillna(row.mean()), axis = 1)
    
