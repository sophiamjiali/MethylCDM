# ==============================================================================
# Script:           betaVAE_objective.py
# Purpose:          Defines the Optuna objective function for training.
# Author:           Sophia Li
# Affiliation:      CCG Lab, Princess Margaret Cancer Center, UHN, UofT
# Date:             12/31/2025
#
# Configurations:   betaVAE.yaml
#
# Notes:            Compatible with an Optuna hyperparameter sweep
# ==============================================================================

import wandb
import lightning.pytorch as pl
from pathlib import Path

from MethylCDM.utils.utils import resolve_path
from MethylCDM.models.betaVAE import BetaVAE
from MethylCDM.data.methylation_datamodule import MethylDataModule
from MethylCDM.constants import BETAVAE_CHECKPOINT_DIR
from MethylCDM.utils.training_utils import (
    configure_callbacks, 
    configure_loggers
)

def objective(trial, config):
    """
    Objective function for an an Optuna hyperparameter optimization trial.

    Parameters
    ----------
    trial (Trial): Optuna Trial object to sample hyperparameter values
    config (Dict): Dictionary containing hyperparameter sweep value ranges

    Returns
    -------
    (float): validation loss to be minimized by Optuna
    """

    # Suggest hyperparameters
    latent_dim = trial.suggest_int("latent_dim", *config['latent_dim'])
    beta = trial.suggest_float("beta", *config['beta'], log = True)
    lr = trial.suggest_float("lr", *config['lr'], log = True)
    batch_size = trial.suggest_int("batch_size", *config['batch_size'])

    # Initialize the model
    model = BetaVAE(
        input_dim = config['input_dim'],
        latent_dim = latent_dim,
        encoder_dims = config['encoder_dims'],
        decoder_dims = config['decoder_dims'],
        beta = beta,
        lr = lr
    )

    # Initialize the DataModule
    datamodule = MethylDataModule(
        train_adata_path = config['train_adata_path'],
        val_adata_path = config["val_adata_path"],
        test_adata_path = config["test_adata_path"],
        batch_size = batch_size,
        num_workers = config['num_workers']
    )

    # Initialize callbacks and loggers
    checkpoint_dir = config.get('checkpoint_dir', '')
    checkpoint_dir = resolve_path(checkpoint_dir, BETAVAE_CHECKPOINT_DIR)
    Path(checkpoint_dir).mkdir(parents = True, exist_ok = True)

    callbacks = configure_callbacks(trial, checkpoint_dir)
    logger = configure_loggers(trial)[0]

    # Initialize the Trainer and train the model
    trainer = pl.Trainer(
        max_epochs = config['max_epochs'],
        callbacks = callbacks,
        logger = logger,
        accelerator = "auto",
        devices = "auto"
    )
    trainer.fit(model, datamodule)

    # Fetch validation metrics
    val_loss = trainer.callback_metrics['val_total_loss'].item()
    wandb.finish()

    return val_loss

