#!/usr/bin/env python3
# ==============================================================================
# Script:           process_methylation.py
# Purpose:          Entry-point to download and preprocess DNA methylation data
#                   for a specified project
# Author:           Sophia Li
# Affiliation:      CCG Lab, Princess Margaret Cancer Center, UHN, UofT
# Date:             11/18/2025
#
# Configurations:   methylation_preproc.yaml
#
# Notes:            Checks to see if data needs to be downloaded by simply 
#                   checking if the raw data folder is empty, does not support
#                   partial downloads
# ==============================================================================

import pandas as pd
from pathlib import Path
import argparse
import os

from MethylCDM.utils.utils import init_environment, load_config, resolve_path
from MethylCDM.data.load_methylation import download_methylation, merge_cohort
from MethylCDM.preprocessing.process_methylation import (
    process_methylation,
    clean_methylation_data
)
from MethylCDM.constants import (
    RAW_METHYLATION_DIR, 
    INTERMEDIATE_METHYLATION_DIR,
    PROCESSED_METHYLATION_DIR,
    METADATA_METHYLATION_DIR
)

def main():

    # -----| Environment Initialization |-----
    
    # Parse the arguments provided to the entry-point script
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", type = str, required = True)
    parser.add_argument("--config_pipeline", type = str, required = True)
    parser.add_argument("--config_preproc", type = str, required = True)
    parser.add_argument("--verbose", type = bool, default = False)
    args = parser.parse_args()

    # Load the relevant configuration files 
    pipeline_cfg = load_config(args.config_pipeline)
    preproc_cfg = load_config(args.config_preproc)

    # Initialize the environment for reproducible analysis
    init_environment(pipeline_cfg)

    # -----| Data Downloading, Loading, and Cleaning |-----

    # Check to see if DNA methylation data needs to be downloaded
    raw_data_dir = preproc_cfg.get('download', {}).get('raw_data_dir', '')
    raw_data_dir = resolve_path(raw_data_dir, RAW_METHYLATION_DIR)
    project_raw_dir = os.path.join(raw_data_dir, args.project)
    Path(project_raw_dir).mkdir(parents = True, exist_ok = True)

    if not os.listdir(project_raw_dir):
        download_methylation(args.project, preproc_cfg, args.verbose)
        clean_methylation_data(project_raw_dir)

    # Load the methylation data and merge it into a cohort-level matrix
    cpg_matrix = merge_cohort(args.project, preproc_cfg)

    inter_data_dir = (preproc_cfg.get('preprocess', {})
                                 .get('intermediate_data_dir', ''))
    inter_data_dir = resolve_path(inter_data_dir, INTERMEDIATE_METHYLATION_DIR)
    project_inter_dir = os.path.join(inter_data_dir, args.project)
    Path(project_inter_dir).mkdir(parents = True, exist_ok = True)

    inter_file = os.path.join(project_inter_dir, 
                              f"{args.project}_cpg_matrix_raw.parquet")
    cpg_matrix.to_parquet(inter_file) 

    # -----| Data Preprocessing |-----
    metadata = pd.read_csv(METADATA_METHYLATION_DIR)
    cpg_matrix = process_methylation(cpg_matrix, metadata, preproc_cfg)

    proc_data_dir = (preproc_cfg.get('preprocess', {})
                               .get('processed_data_dir', ''))
    proc_data_dir = resolve_path(proc_data_dir, PROCESSED_METHYLATION_DIR)
    proc_file = os.path.join(proc_data_dir, 
                             f"{args.project}_cpg_matrix_processed.parquet")

    cpg_matrix.to_parquet(proc_file)