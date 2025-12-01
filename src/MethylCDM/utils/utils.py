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

from MethylCDM.constants import CONFIG_DIR

# =====| File I/O Utilities |===================================================

def resolve_path(path_str, default_path):
    """
    Resolves and returns the path. If a relative path is provided through
    a file or directory name, it is automatically resolved relative to the 
    project root. Else, the absolute path is returned as provided.

    Parameters
    ----------
    path_str (str): path to a YAML configuration file
    default_path (str): default path (constant) to the project root

    Returns
    -------
    path (Path): resolved Path object from pathlib

    Raises
    ------
    FileNotFoundError: if the file does not exist at the specified path
    """

    # Resolve to the project's root if the provided path is not absolute
    p = Path(path_str)
    return p.resolve() if p.is_absolute() else (default_path / p).resolve()


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
    

# =====| Configuration & Environment |==========================================

def init_environment(config):
    """
    Initializes the current runtime environment for reproducibility.
    
    Parameters
    ----------
    config : Configuration object containing:
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
    path = resolve_path(path_str, CONFIG_DIR)

    if not path.exists():
        raise FileNotFoundError(f"Configuration file not found at {path}.")

    with open(path) as f:
        config = yaml.safe_load(f)

    if not isinstance(config, dict):
        raise ValueError(
            f"Configuration file {path} did not return a dictionary."
        )

    return config