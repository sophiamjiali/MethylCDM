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

def reconcile_methylation(data_dir, verbose):
    """
    Reconciles multiple AnnData objects in a directory into a single AnnData. 
    
    Normalizes and produces a comprehensive AnnData object from all AnnData
    objects in `data_dir`. Converts each into sparse format to reduce memory
    usage.

    Parameters
    ----------
    data_dir (str): path to a processed data directory containing AnnData 
                    objects of multiple project(s)
    verbose (bool): toggle for verbose output

    Returns
    -------
    (AnnData): the concatenated and normalized AnnData object containing all
               projects.
    """

    if verbose:
        print("=" * 50)
        print(f"Beginning to reconcile methylation data")
        print("=" * 50)
        print("\n")

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

        if verbose: print(f"Loaded data for project {file.stem.split('_')[0]}")

    # Concatenate along cells (obs), taking the gene intersection
    cohort_adata = ad.concat(
        adatas, join = "inner", label = "batch", 
        keys = [f.stem.split('_')[0] for f in adata_files]
    )

    if verbose:
        print("=" * 50)
        print(f"Completed methylation data reconciliation")
        print("=" * 50)
        print("\n")

    return cohort_adata

def split_cohort(cohort_anndata, seed, verbose):
    """
    Split a single AnnData object into train, validation, and test sets
    using stratified splits for TCGA project using a 60%-20%-20% ratio.

    Parameters
    ----------
    cohort_anndata (AnnData): Full cohort AnnData object
    seed (int): random state for reproducibility
    verbose (bool): toggle for verbose output

    Returns
    -------
    train_adata, val_adata, test_adata (AnnData): stratified splits.
    """

    if verbose:
        print("=" * 50)
        print(f"Beginning train-val-test splitting")
        print("=" * 50)
        print("\n")

    # Split into train+val and test
    train_idx, test_idx = train_test_split(
        cohort_anndata.obs_names,
        test_size = 0.2,
        stratify = cohort_anndata.obs['tcga_project'],
        random_state = seed,
        shuffle = True
    )

    # Split into train and validation
    train_idx, val_idx = train_test_split(
        train_idx,
        test_size = 0.25,
        stratify = cohort_anndata.obs.loc[train_idx, 'tcga_project'],
        random_state = seed,
        shuffle = True
    )

    # Slice the cohort AnnData object
    train_adata = cohort_anndata[train_idx].copy()
    val_adata = cohort_anndata[val_idx].copy()
    test_adata = cohort_anndata[test_idx].copy()

    if verbose:
        print("=" * 50)
        print(f"Completed train-val-test splitting")
        print("=" * 50)
        print("\n")

    return train_adata, val_adata, test_adata
