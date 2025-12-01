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

import argparse
import os

from MethylCDM.utils.utils import init_environment, load_config
from MethylCDM.data.load_methylation import (
    download_methylation, 
    load_raw_methylation
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

    # -----| Data Downloading and Loading |-----

    # Check to see if DNA methylation data needs to be downloaded
    raw_data_dir = preproc_cfg.get('download', {}).get('raw_data_dir', '')
    project_raw_dir = os.path.join(
        raw_data_dir, f"{args.project}_raw_methylation"
    )

    if not os.listdir(project_raw_dir):
        download_methylation(args.project, preproc_cfg, args.verbose)

    # Load the downloaded methylation data as a list of DataFrames (per project)
    data = load_raw_methylation(args.project, preproc_cfg, args.verbose)

    # -----| Data Preprocessing |-----

    # preprocess data

    # -----| Saving Preprocessed Data |-----

    # save preprocessed data