#!/usr/bin/env python3
# ==============================================================================
# Script:           train_betaVAE.py
# Purpose:          Entry-point to train a betaVAE using a pancancer dataset
# Author:           Sophia Li
# Affiliation:      CCG Lab, Princess Margaret Cancer Center, UHN, UofT
# Date:             12/31/2025
#
# Configurations:   betaVAE.yaml
#
# Notes:            Begins an Optuna hyperparameter sweep
# ==============================================================================

import os
import argparse
import torch
import numpy as np
import pandas as pd
from torch.optim import Adam
from torch.utils.data import DataLoader
from tensorboardX import SummaryWriter
from sklearn.model_selection import train_test_split
from warmup_scheduler import GradualWarmupScheduler

from MethylCDM.utils.utils import init_environment, load_config, resolve_path

def main():

    # -----| Environment Initialization |-----

    # Parse the arguments provided to the entry-point script
    parser = argparse.ArgumentParser()
    parser.add_argument("--config_pipeline", type = str, required = True)
    parser.add_argument("--config_train", type = str, required = True)
    parser.add_argument("--verbose", type = bool, default = False)
    args = parser.parse_args()


if __name__ == "__main__":
    main()