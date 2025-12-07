# ==============================================================================
# Script:           utils.py
# Purpose:          Utility functions for configuration and initialization
# Author:           Sophia Li
# Affiliation:      CCG Lab, Princess Margaret Cancer Center, UHN, UofT
# Date:             11/18/2025
#
# Configurations:   pipeline.yaml
# ==============================================================================

import random
import numpy as np
import pandas as pd
import torch
import yaml
from pathlib import Path

from MethylCDM.constants import (
    CONFIG_DIR,
    ANNOTATION_27K,
    ANNOTATION_450K,
    ANNOTATION_EPIC
)

# =====| File I/O Utilities |===================================================

def resolve_path(path_str, default_path, build_path = False):
    """
    Resolves and returns the path. If a relative path is provided through
    a file or directory name, it is automatically resolved relative to the 
    project root. Else, the absolute path is returned as provided.

    Parameters
    ----------
    path_str (str): path to a YAML configuration file
    default_path (str): default path (constant) to the project root
    build_path (boolean): appends the path to the default if toggled

    Returns
    -------
    path (Path): resolved Path object from pathlib
    """

    p = Path(path_str)

    if p.is_absolute():
        return p.resolve()
    elif build_path:
        return (default_path / path_str).resolve()
    else:
        return default_path.resolve()


def build_meta_fields(fields):
    meta = []
    for f in fields:
        if '.' in f:
            parts = f.split('.')
            if len(parts) == 1:
                meta.append(parts[0])
            else:
                meta.append(parts)
    return meta

def load_annotations(manifests):
    """
    Loads and returns the dominant Illumina methylation manifest
    (EPIC > 450K > 27K) out of the manifests present in the dataset as
    provided by `manifests`, standardizing columns for quality control. 

    Parameters
    ----------
    manifests (list): list of strings of manifests present in the dataset

    Returns
    -------
    (DataFrame): the dominant manifest present in `manifests`

    Raises
    ------
    ValueError: if no viable manifest was provided
    """
    
    # Fetch the dominant manifest present in `manifests`
    if ("Illumina Human Methylation EPIC" in manifests):
        manifest = pd.read_csv(ANNOTATION_EPIC)
    elif ("Illumina Human Methylation 450" in manifests):
        manifest = pd.read_csv(ANNOTATION_450K)
    elif ("Illumina Human Methylation 27" in manifests):
        manifest = pd.read_csv(ANNOTATION_27K)
    else:
        raise ValueError ("No valid manifest provided in `manifests`.")
    
    # Standardize column names, ignoring if they are missing
    col_map = {
        "ID": "probe_id",
        "IlmnID": "probe_id",
        "CHR": "chr",
        "Chromosome": "chr",
        "MAPINFO": "pos",
        "Coordinate_37": "pos",
        "STRAND": "strand",
        "Strand": "strand",
        "Type": "probe_type",
        "Infinium_Design_Type": "probe_type",
        "Relation_to_UCSC_CpG_Island": "island",
        "Relation_to_Island": "island",
        "SNP_ID": "is_snp",
        "Probe_SNPs": "is_snp",
        "MASK_general": "is_cross_reactive",
        "MASK_crosshyb": "is_cross_reactive",
    }
    manifest.rename(columns = {
        k: v for k, v in col_map.items() if k in manifest.columns
    }, inplace = True)

    # Fill any missing values with NA (i.e. if data was not available)
    manifest["is_snp"] = manifest.get(
        "is_snp", pd.Series([False] * len(manifest))
    )
    manifest["is_cross_reactive"] = manifest.get(
        "is_cross_reactive", pd.Series([False] * len(manifest))
    )

    # Reduce the manifest to the standardized set
    keep_cols = ["probe_id", "chr", "pos", "strand", "probe_type", 
                 "island","is_snp", "is_cross_reactive"]
    
    return manifest[keep_cols]
    

# =====| Configuration & Environment |==========================================

def init_environment(config):
    """
    Initializes the current runtime environment for reproducibility.
    
    Parameters
    ----------
    config : a configuration object containing:
        - seed (int): integer value for reproducibility
    """

    # Fetch all relevant values from the configurations object
    seed = config.get('seed', -1)

    # Set the seed for all appropriate packages of the pipeline
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)


def load_config(path_str):
    """
    Loads and returns the configuration file provided by the path. If a relative
    path is provided (filename), automatically resolves it relative to the 
    project root.

    Parameters
    ----------
    path_str (str): path to a YAML configuration file

    Returns
    -------
    config (dict): dictionary of configuration values

    Raises
    ------
    FileNotFoundError: if the file does not exist at the specified path
    ValueError: if the YAML file cannot be parsed into a dictionary
    """

    # Resolve to the project's root if the provided path is not absolute
    path = resolve_path(path_str, CONFIG_DIR, build_path = True)

    if not path.exists():
        raise FileNotFoundError(f"Configuration file not found at {path}.")

    with open(path) as f:
        config = yaml.safe_load(f)

    if not isinstance(config, dict):
        raise ValueError(
            f"Configuration file {path} did not return a dictionary."
        )

    return config
