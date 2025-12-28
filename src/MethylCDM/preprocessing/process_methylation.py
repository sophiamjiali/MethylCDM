# ==============================================================================
# Script:           process_methylation.py
# Purpose:          Performs quality control and preprocesses methylation data
# Author:           Sophia Li
# Affiliation:      CCG Lab, Princess Margaret Cancer Center, UHN, UofT
# Date:             11/25/2025
# ==============================================================================

from pathlib import Path
import pandas as pd
import numpy as np
import anndata as ad
from collections import defaultdict
import os
from MethylCDM.utils.utils import (
    resolve_path,
    load_beta_file,
    load_annotation,
    load_cpg_matrix
)
from MethylCDM.constants import RAW_METHYLATION_DIR

# =====| Preprocessing Wrapper |================================================

def process_methylation(project, metadata, config):
    """
    Preprocesses DNA methylation beta values for a given project, returning a 
    gene-level matrix stored as an AnnData object. 
    
    Raw samples are processed by array type to perform probe-level quality c
    ontrol (QC) and aggregate to gene-level values. The final output is samples x gene matrix with aligned metadata suitable for downstream analyses.

    Metadata from the provided DataFrame is aligned to the surviving samples 
    and stored in the AnnData `obs` attribute. Array type is captured from the
    `platform` column.

    Parameters
    ----------
    project (str): name of a TCGA project with data available on GDC
    metadata (DataFrame): sample metadata including at minimum the `id` column 
        corresponding to the sample identifiers (UUIDs)
    config (dict): Configuration object containing:
        - raw_data_dir (str): path to the output raw data directory

    Returns 
    -------
    (adata): Gene-level methylation matrix with samples as rows and genes as 
        columns. Sample metadata is stored in `adata.obs`.
    """

    # Resolve the project's raw data directory
    raw_data_dir = config.get('download', {}).get('raw_data_dir', '')
    raw_data_dir = resolve_path(raw_data_dir, RAW_METHYLATION_DIR)
    project_data_dir = os.path.join(raw_data_dir, f"{project}")

    # Verify the raw data exists and is not empty
    if not os.path.isdir(project_data_dir):
        raise FileNotFoundError(f"Raw data directory was not found at "
                                f"{project_data_dir}.")
    if not os.listdir(project_data_dir):
        raise FileNotFoundError(f"Raw data directory was empty at "
                                f"{project_data_dir}.")
    
    # Identify and load all nested beta value .txt files
    beta_files = [
        f for f in Path(project_data_dir).glob("*.level3betas.parquet")
    ]

    # Keep the first occurrence of duplicates in the metadata (redundant)
    metadata = metadata.set_index('file_name')
    metadata = metadata[~metadata.index.duplicated(keep = 'first')]

    # Identify the highest coverage array type
    manifests = metadata['platform'].unique()
    annotation, array_type = load_annotation(manifests)

    # Normalize the metadata extensions
    metadata.index = metadata.index.str.removesuffix(".txt")
    valid_stems = set(metadata[metadata['platform'] == array_type].index)

    # Fetch samples that align with the given array type
    beta_files = [f for f in beta_files if f.stem in valid_stems]

    # Preprocess the beta values into a gene-level matrix
    cpg_matrix = load_cpg_matrix(beta_files)
    gene_matrix = process_array_methylation(cpg_matrix, annotation, config)
    gene_matrix = gene_matrix.sort_index()

    # Filter the metadata for surviving samples
    metadata = metadata.loc[gene_matrix.columns]
    metadata = metadata.sort_index()

    # Initialize the gene matrix as an AnnData object
    adata = ad.AnnData(X = gene_matrix.T)
    adata.obs = metadata

    return adata


def process_array_methylation(cpg_matrix, annotation, config):
    """
    Preprocesses a CpG matrix of DNA methylation beta values corresponding to
    a single array-type/platform. Performs sample-level and probe-level quality
    control (QC) and imputation using probe-wise mean if toggled in the 
    associated configuration file.

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

    # 1. Probe quality control
    if preproc_cfg.get('toggle_probe_filtering', False):
        cpg_matrix = probe_qc(cpg_matrix, annotation, max_missing_probe)

    # 2. Sample quality control
    if preproc_cfg.get('toggle_sample_filtering', False):
        cpg_matrix = sample_qc(cpg_matrix, max_missing_sample)

    # 3. Missing value imputation
    if preproc_cfg.get('toggle_imputation', False):
        cpg_matrix = impute_missing(cpg_matrix)

    # 5. Aggregate probes to gene-level
    gene_matrix = aggregate_genes(cpg_matrix, annotation)

    # 6. Clip extreme beta values
    gene_matrix = clip_beta_values(gene_matrix)

    return gene_matrix


# =====| Preprocessing Helpers |================================================

def probe_qc(cpg_matrix, annotation, max_missing):
    """
    Filters low-quality probes with too many missing values, variability,
    SNP/cross-reactive annotations, and sex chromosome probes. Additionally 
    removes non-standard probe IDs that indicate control probes or artifacts.

    Parameters
    ----------
    cpg_matrix (DataFrame): a CpG x Sample ID matrix of beta values
    annotations (DataFrame): probe annotations for the given array type
    max_missing (float): maximum allowed fraction of missing samples per probe

    Returns
    -------
    cpg_matrix (DataFrame): the cleaned CpG matrix 
    """

    # Remove non-standard probes
    standard_probes = annotation[
        annotation['probe_id'].str.startswith('cg')
    ]['probe_id']
    cpg_matrix = cpg_matrix[cpg_matrix.index.isin(standard_probes)]
    
    # Filter by missingness
    missing_frac = cpg_matrix.isna().mean(axis = 1)
    keep_missing = missing_frac <= max_missing

    # Filter by annotation flags
    annotation = annotation.set_index('probe_id')
    annotation = annotation.loc[cpg_matrix.index]

    keep_annotation = (
        ~annotation['is_sex_chr'] &
        ~annotation['has_cpg_snp'] &
        ~annotation['has_sbe_snp'] &
        ~annotation['has_probe_snp'] &
        ~annotation['is_cross_reactive'] &
        ~annotation['is_multi_mapped']
    )

    # Combine filters and subset the matrix
    keep_annotation = keep_annotation.loc[cpg_matrix.index]
    keep = keep_missing & keep_annotation

    return cpg_matrix.loc[keep]

def sample_qc(cpg_matrix, max_missing):
    """ 
    Removes samples based on missingness and global beta distribution checks. 
    
    Parameters 
    ---------- 
    beta_values (DataFrame): a CpG x Sample ID matrix of beta values 
    max_missing (float): maximum allowed fraction of missing CpGs per sample 
    
    Returns 
    ------- 
    (DataFrame): filtered matrix with samples below the threshold removed 
    """ 
    
    # Filter samples for missingness above the provided threshold 
    sample_missing = cpg_matrix.isna().mean(axis = 0) 
    missing_mask = sample_missing[sample_missing > max_missing].index.tolist() 
    
    # Filter for global distribution per sample, flagging beyond three SD 
    sample_mean = cpg_matrix.mean(skipna = True) 
    sample_std = cpg_matrix.std(skipna = True) 
    
    mean_thresh = sample_mean.mean() + np.array([-3, 3]) * sample_mean.std() 
    std_thresh = sample_std.mean() + np.array([-3, 3]) * sample_std.std() 
    
    sd_mask = sample_mean[ 
        (sample_mean < mean_thresh[0]) | (sample_mean > mean_thresh[1]) 
    ].index.tolist() 
    sd_mask += sample_std[ 
        (sample_std < std_thresh[0]) | (sample_std > std_thresh[1]) 
    ].index.tolist() 
    
    # Filter for the combined masks 
    samples_mask = list(set(missing_mask + sd_mask)) 
    return cpg_matrix.drop(columns = samples_mask)


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
    

def aggregate_genes(cpg_matrix, annotation):
    """
    Aggregates a CpG-level beta matrix to the gene level using mean beta values.

    Parameters
    ----------
    cpg_matrix (DataFrame): a CpG x Sample ID matrix of beta values
    annotations (DataFrame): probe annotations for the given array type

    Returns
    -------
    (DataFrame): the gene-level matrix
    """

    # Subset the annotations for the probes present in the matrix
    annotation = annotation.set_index("probe_id")
    annotation = annotation.loc[cpg_matrix.index]

    # Explode multi-gene mappings into long format
    multi_gene_series = annotation["gene_symbol"].str.split(';')
    exploded_mapping = multi_gene_series.explode()

    # Align CpG values to exploded gene mapping
    cpg_long = cpg_matrix.loc[exploded_mapping.index].copy()
    cpg_long['gene'] = exploded_mapping.values

    # Group by gene and aggregate
    gene_matrix = cpg_long.groupby('gene').mean()

    return gene_matrix


def clip_beta_values(gene_matrix):
    """
    Clips gene-level beta values to avoid exact 0 or 1 for numerical 
    stability when training a beta variational autoencoder.

    Parameters
    ----------
    gene_matrix (DataFrame): gene x sample matrix with beta values in [0,1]

    Returns
    -------
    (DataFrame): the gene-level matrix with extreme beta values clipped
    """

    return gene_matrix.clip(lower = 1e-5, upper = 1 - 1e-5)