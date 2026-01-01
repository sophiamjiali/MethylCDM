#!/usr/bin/env python3
# ==============================================================================
# Script:           prepare_data.py
# Purpose:          Entry-point to split and reconcile DNA methylation data
# Author:           Sophia Li
# Affiliation:      CCG Lab, Princess Margaret Cancer Center, UHN, UofT
# Date:             12/31/2025
#
# Configurations:   methylation_preproc.yaml
# ==============================================================================

import os
import argparse
import anndata
from pathlib import Path

from MethylCDM.utils.utils import init_environment, load_config, resolve_path
from MethylCDM.preprocessing.reconcile_methylation import (
    reconcile_methylation,
    split_cohort
)
from MethylCDM.constants import (
    PROCESSED_METHYLATION_DIR,
    TRAINING_METHYLATION_DIR
)

def main():

    # -----| Environment Initialization |-----

    # Parse the arguments provided to the entry-point script
    parser = argparse.ArgumentParser()
    parser.add_argument("--config_pipeline", type = str, required = True)
    parser.add_argument("--config_betaVAE", type = str, required = True)
    parser.add_argument("--verbose", type = bool, default = False)
    args = parser.parse_args()

    # Load the relevant configuration file
    pipeline_cfg = load_config(args.config_pipeline)
    model_cfg = load_config(args.config_betaVAE)

    # Initialize the environment for reproducible analysis
    init_environment(pipeline_cfg)

    if args.verbose:
        print("=" * 50)
        print(f"~~~~~| Beginning Step 2")
        print("=" * 50)
        print("\n")

    # -----| Data Reconciliation |-----
    data_dir = model_cfg.get('project_data_dir', '')
    data_dir = resolve_path(data_dir, PROCESSED_METHYLATION_DIR)
    cohort_adata = reconcile_methylation(data_dir, args.verbose, 
                                         pipeline_cfg.get('seed', 42))

    # Save the full cohort AnnData Object
    train_data_dir = model_cfg.get('training_data_dir', '')
    train_data_dir = resolve_path(train_data_dir, TRAINING_METHYLATION_DIR)
    Path(train_data_dir).mkdir(parents = True, exist_ok = True)
    cohort_path = os.path.join(train_data_dir, "tcga_cohort_gene_matrix.h5ad")
    anndata.settings.allow_write_nullable_strings = True
    cohort_adata.write_h5ad(cohort_path, compression = "gzip")

    # -----| Train-Validation-Test Split |-----
    train_adata, val_adata, test_adata = split_cohort(
        cohort_adata, pipeline_cfg.get('seed', 42)
    )

    # Save each split individually
    train_path = os.path.join(train_data_dir, "tcga_train_gene_matrix.h5ad")
    val_path = os.path.join(train_data_dir, "tcga_val_gene_matrix.h5ad")
    test_path = os.path.join(train_data_dir, "tcga_test_gene_matrix.h5ad")
    train_adata.write_h5ad(train_path, compression = "gzip")
    val_adata.write_h5ad(val_path, compression = "gzip")
    test_adata.write_h5ad(test_path, compression = "gzip")

    if args.verbose:
        print("=" * 50)
        print(f"~~~~~| Finished Step 2")
        print("=" * 50)

if __name__ == "__main__":
    main()