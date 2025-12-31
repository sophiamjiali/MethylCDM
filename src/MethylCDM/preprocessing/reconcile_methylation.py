# ==============================================================================
# Script:           reconcile_methylation.py
# Purpose:          Reconciles DNA methylation data across all projects
# Author:           Sophia Li
# Affiliation:      CCG Lab, Princess Margaret Cancer Center, UHN, UofT
# Date:             12/31/2025
# ==============================================================================

import os
import anndata as ad
from pathlib import Path
import scipy.sparse as sp
from sklearn.model_selection import train_test_split

def reconcile_methylation(data_dir):
    """
    Reconciles multiple AnnData objects in a directory into a single AnnData. 
    
    Normalizes and produces a comprehensive AnnData object from all AnnData
    objects in `data_dir`. Converts each into sparse format to reduce memory
    usage.

    Parameters
    ----------
    data_dir (str): path to a processed data directory containing AnnData 
                    objects of multiple project(s)

    Returns
    -------
    (AnnData): the concatenated and normalized AnnData object containing all
               projects.
    """

    # Load each project AnnData object
    adata_files = list(Path(data_dir).glob("*.h5ad"))
    
    adatas = []
    for file in adata_files:

        # Convert to sparse if not already done
        adata = ad.read_h5ad(file)
        if not sp.issparse(adata.X):
            adata.X = sp.csr_matrix(adata.X)
        
        adata.obs['tcga_project'] = file.stem.split('_')[0]
        adatas.append(adata)

    # Concatenate along cells (obs), taking the gene intersection
    cohort_adata = ad.concat(
        adatas, join = "outer", label = "batch", 
        keys = [f.stem.split('_')[0] for f in adata_files]
    )

    return cohort_adata

def split_cohort(cohort_anndata, seed):
    """
    Split a single AnnData object into train, validation, and test sets
    using stratified splits for TCGA project using a 60%-20%-20% ratio.

    Parameters
    ----------
    cohort_anndata (AnnData): Full cohort AnnData object

    Returns
    -------
    train_adata, val_adata, test_adata (AnnData): stratified splits.
    """

    # Split into train+val and test
    train_idx, test_idx = train_test_split(
        cohort_anndata.obs_names,
        test_size = 0.2,
        stratify = cohort_anndata.obs['tcga_project'],
        random_state = seed
    )

    # Split into train and validation
    train_idx, val_idx = train_test_split(
        train_idx,
        test_size = 0.25,
        stratify = cohort_anndata.obs.loc[train_idx, 'tcga_project'],
        random_state = seed
    )

    # Slice the cohort AnnData object
    train_adata = cohort_anndata[train_idx].copy()
    val_adata = cohort_anndata[val_idx].copy()
    test_adata = cohort_anndata[test_idx].copy()

    return train_adata, val_adata, test_adata
