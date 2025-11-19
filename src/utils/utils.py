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
import torch
import yaml
import os

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


def load_config(path):
    """
    Loads and returns the configuration file provided by the path.

    Parameters
    ----------
    path (str): path to a YAML configuration file

    Returns
    -------
    config (dict): dictionary of configuration values

    Raises
    ------
    FileNotFoundError: if the file does not exist at the specified path
    ValueError: if the YAML file cannot be parsed into a dictionary
    """

    if not os.path.isfile(path):
        raise FileNotFoundError(f"Configuration file not found at {path}.")

    with open(path) as f:
        config = yaml.safe_load(f)

    if not isinstance(config, dict):
        raise ValueError(
            f"Configuration file {path} did not return a dictionary."
        )

    return config