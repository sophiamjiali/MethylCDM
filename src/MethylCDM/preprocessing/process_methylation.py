# ==============================================================================
# Script:           process_methylation.py
# Purpose:          Performs quality control and preprocesses methylation data
# Author:           Sophia Li
# Affiliation:      CCG Lab, Princess Margaret Cancer Center, UHN, UofT
# Date:             11/25/2025
# ==============================================================================

import pandas as pd
import numpy as np
import os
from sklearn.impute import KNNImputer

from MethylCDM.constants import ANNOTATION_27K, ANNOTATION_450K
from MethylCDM.utils.utils import resolve_path

# =====| Preprocessing Wrapper |================================================

def process_methylation(cpg_matrix, metadata, config):
    """
    Preprocesses and returns DNA methylation data in the form of a CpG x Sample
    ID matrix. Performs sample-level and probe-level quality control (QC) and
    imputation using k-nearest neighbours (KNN) if toggled in the associated 
    configuration file.

    The following preprocessing is applied if toggled:
        - Sample-level QC:  remove samples with > 5% CpGs missing
        - Probe-level QC:   remove CpGs with > 5% values missing, SNP-affected, 
                            cross-reactive, and sex chromosome probes
        - Imputation:       impute missing CpG beta values using KNN imputation

    Parameters
    ----------
    cpg_matrix (DataFrame): a CpG x Sample ID matrix of beta values
    metadata (DataFrame): a metadata matrix corresponding to the beta values
    config (dict): a configuration object containing:
        - toggle_sample_filtering (boolean): perform sample filtering
        - toggle_probe_filtering (boolean): perform probe filtering
        - togge_imputation (boolen): perform imputation

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

    # 1. Sample quality control
    if config.get('toggle_sample_filtering', False):

        cpg_matrix = sample_qc(cpg_matrix, metadata, max_missing_sample)
    
    # 2. Probe quality control
    if config.get('toggle_probe_filtering', False):

        # Load the probe annotation manifests
        annotation_27k = preproc_cfg.get('annotation_27k', '')
        annotation_450k = preproc_cfg.get('annotation_450k', '')
        annotation_27k = resolve_path(annotation_27k, ANNOTATION_27K)
        annotation_450k = resolve_path(annotation_450k, ANNOTATION_450K)
        
        if not os.path.exists(annotation_27k):
            raise FileNotFoundError(f"Illumina 27K probe annotation was not"
                                    f"found at {annotation_27k}.")
        if not os.path.exists(annotation_450k):
            raise FileNotFoundError(f"Illumina 240K probe annotation was not"
                                    f"found at {annotation_27k}.")
        
        annotations = {
            "27k": pd.read_csv(annotation_27k, skiprows = 7),
            "450k": pd.read_csv(annotation_450k, skiprows = 7)
        }

        cpg_matrix = probe_qc(cpg_matrix, metadata, annotations, 
                              max_missing_probe)

    # 3. Missing value imputation
    if config.get('toggle_imputation', False):
        cpg_matrix = impute_missing(cpg_matrix, n_neighbours)

    return cpg_matrix


# =====| Preprocessing Helpers |================================================

def sample_qc(cpg_matrix, metadata, max_missing):
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
    keep_samples = cpg_matrix.columns[cpg_matrix.isnull().mean() <= max_missing]
    return cpg_matrix[keep_samples]


def probe_qc(cpg_matrix, metadata, annotations, max_missing):
    """
    Removes low-quality probes with too many missing values.

    Parameters
    ----------
    cpg_matrix (DataFrame): a CpG x Sample ID matrix of beta values
    metadata (DataFrame):
    annotations (dict):
    max_missing (float): maximum allowed fraction of missing samples per probe


    """
    return

def impute_missing(cpg_matrix, n_neighbours):
    return
    
