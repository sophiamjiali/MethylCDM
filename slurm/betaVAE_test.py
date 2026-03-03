# ==============================================================================
# Script:           train_betaVAE.py
# Purpose:          Single training run for BetaVAE on TCGA methylation data.
#                   Intended for sanity checking before a full Optuna sweep.
# Author:           Sophia Li
# Affiliation:      CCG Lab, Princess Margaret Cancer Center, UHN, UofT
# ==============================================================================

import os
import gc
import argparse
import yaml
import torch
import wandb
import pytorch_lightning as pl
from pathlib import Path
import sys

project_root = Path('/ddn_exa/campbell/sli/methylcdm-project')
sys.path.append(str(project_root / "src"))

from MethylCDM.models.betaVAE import BetaVAE
from MethylCDM.data.methylation_datamodule import MethylDataModule
from MethylCDM.utils.training_utils import configure_callbacks, configure_loggers
from MethylCDM.utils.utils import resolve_path
from MethylCDM.constants import BETAVAE_CHECKPOINT_DIR


# ==============================================================================
# Sensible single-run hyperparameters
# ==============================================================================
# These are chosen to be representative of the mid-range of each sweep
# dimension — not the best configuration, but a reasonable one that should
# train stably and expose any issues before the full sweep.
#
# Architecture rationale:
#   encoder_dims [1024, 256, 128]: middle candidate; ~3× compression per
#       layer from 211580 → 1024 → 256 → 128 → latent_dim.
#   latent_dim 128: mid-range; sufficient for pan-cancer tissue/subtype signal.
#   decoder_dims [256, 1024]: automatic mirror of encoder (minus last dim).
#
# Regularisation rationale:
#   beta 0.005: mid-range; conservative enough to avoid collapse on first run.
#   input_dropout 0.1: lower bound; preserves co-methylation structure.
#   num_cycles 4: more cycles → better latent utilisation at 200 epochs.
#
# Optimisation rationale:
#   lr 1e-3: slightly below RNA-CDM default (3e-3); MSE on M-values produces
#       larger raw gradients than normalised RNA counts, so a conservative LR
#       reduces the risk of instability on the first run.
#   batch_size 128: mid-range; stable BatchNorm statistics at n=1214.

SINGLE_RUN_CONFIG = {
    # Architecture
    "input_dim":     211580,
    "latent_dim":    128,
    "encoder_dims":  [1024, 256, 128],
    "decoder_dims":  [256, 1024],
    # Regularisation
    "beta":          0.005,
    "input_dropout": 0.1,
    "num_cycles":    4,
    # Optimisation
    "lr":            1e-3,
    "batch_size":    128,
    # Training
    "max_epochs":    200,
}


def parse_args():
    parser = argparse.ArgumentParser(
        description="Single BetaVAE training run for TCGA methylation data."
    )
    parser.add_argument(
        "--config", type=str, required=True,
        help="Path to beta_vae.yaml configuration file."
    )
    parser.add_argument(
        "--run_name", type=str, default="single_run_baseline",
        help="W&B run name and checkpoint subdirectory name."
    )
    # Optional overrides for quick iteration without editing the config
    parser.add_argument("--latent_dim",    type=int,   default=None)
    parser.add_argument("--beta",          type=float, default=None)
    parser.add_argument("--lr",            type=float, default=None)
    parser.add_argument("--batch_size",    type=int,   default=None)
    parser.add_argument("--max_epochs",    type=int,   default=None)
    parser.add_argument("--num_cycles",    type=int,   default=None)
    parser.add_argument("--input_dropout", type=float, default=None)
    # Sanity check modes
    parser.add_argument(
        "--fast_dev_run", action="store_true",
        help="Run one train/val/test batch only (Lightning fast_dev_run). "
             "Use to verify the full pipeline before training."
    )
    parser.add_argument(
        "--overfit_batches", type=int, default=0,
        help="Number of batches to overfit on. Set to 1 to verify the model "
             "can memorise a single batch. 0 = disabled (default)."
    )
    return parser.parse_args()


def load_config(config_path: str) -> dict:
    with open(config_path, "r") as f:
        return yaml.safe_load(f)


def build_run_config(yaml_config: dict, args) -> dict:
    """
    Merges SINGLE_RUN_CONFIG defaults with yaml_config path/infra settings,
    then applies any CLI overrides. Priority: CLI > SINGLE_RUN_CONFIG > yaml.
    """
    config = {**SINGLE_RUN_CONFIG}

    # Pull infrastructure settings from yaml (paths, num_workers, etc.)
    infra_keys = [
        "train_adata_path", "val_adata_path", "test_adata_path",
        "num_workers", "experiment_dir", "checkpoint_dir"
    ]
    for key in infra_keys:
        if key in yaml_config:
            config[key] = yaml_config[key]

    # Apply CLI overrides for quick iteration
    cli_overrides = {
        "latent_dim":    args.latent_dim,
        "beta":          args.beta,
        "lr":            args.lr,
        "batch_size":    args.batch_size,
        "max_epochs":    args.max_epochs,
        "num_cycles":    args.num_cycles,
        "input_dropout": args.input_dropout,
    }
    for key, val in cli_overrides.items():
        if val is not None:
            config[key] = val

    return config


def main():
    args = parse_args()

    # ------------------------------------------------------------------
    # Config
    # ------------------------------------------------------------------
    yaml_config = load_config(args.config)
    config      = build_run_config(yaml_config, args)

    pl.seed_everything(42, workers=True)

    print("\n========== Run Configuration ==========")
    for k, v in config.items():
        print(f"  {k}: {v}")
    print("=======================================\n")

    # ------------------------------------------------------------------
    # Checkpoint directory
    # ------------------------------------------------------------------
    checkpoint_dir = resolve_path(
        config.get("checkpoint_dir", ""), BETAVAE_CHECKPOINT_DIR
    )
    checkpoint_dir = os.path.join(checkpoint_dir, args.run_name)
    Path(checkpoint_dir).mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------
    # Model
    # ------------------------------------------------------------------
    model = BetaVAE(
        input_dim     = config["input_dim"],
        latent_dim    = config["latent_dim"],
        encoder_dims  = config["encoder_dims"],
        decoder_dims  = config["decoder_dims"],
        beta          = config["beta"],
        input_dropout = config["input_dropout"],
        num_cycles    = config["num_cycles"],
        lr            = config["lr"],
    )

    # ------------------------------------------------------------------
    # DataModule
    # ------------------------------------------------------------------
    datamodule = MethylDataModule(
        train_adata_path = config["train_adata_path"],
        val_adata_path   = config["val_adata_path"],
        test_adata_path  = config["test_adata_path"],
        batch_size       = config["batch_size"],
        num_workers      = config["num_workers"],
    )

    # ------------------------------------------------------------------
    # Callbacks and logger
    # ------------------------------------------------------------------
    # Pass trial=None — this is a single run, not an Optuna trial,
    # so no pruning callback will be added
    callbacks = configure_callbacks(trial=None, checkpoint_dir=checkpoint_dir)
    loggers   = configure_loggers(trial=None, study_name=args.run_name)

    # ------------------------------------------------------------------
    # Trainer
    # ------------------------------------------------------------------
    trainer_kwargs = dict(
        max_epochs        = config["max_epochs"],
        callbacks         = callbacks,
        logger            = loggers,
        accelerator       = "auto",
        devices           = "auto",
        deterministic     = True,
        log_every_n_steps = 1,
    )

    # Sanity check modes — mutually exclusive; fast_dev_run takes priority
    if args.fast_dev_run:
        print("\n[INFO] Running fast_dev_run — one train/val/test batch only.\n")
        trainer_kwargs["fast_dev_run"] = True
        # fast_dev_run overrides max_epochs internally; no need to set it
    elif args.overfit_batches > 0:
        print(f"\n[INFO] Overfitting on {args.overfit_batches} batch(es). "
              f"train_recon should collapse to near zero.\n")
        trainer_kwargs["overfit_batches"] = args.overfit_batches

    trainer = pl.Trainer(**trainer_kwargs)

    # ------------------------------------------------------------------
    # Train
    # ------------------------------------------------------------------
    try:
        trainer.fit(model, datamodule)

        # Skip test and metric reporting in sanity check modes
        if not args.fast_dev_run and args.overfit_batches == 0:
            trainer.test(model, datamodule)

            val_loss  = trainer.callback_metrics.get("val_loss")
            test_loss = trainer.callback_metrics.get("test_loss")
            print("\n========== Final Metrics ==========")
            print(f"  val_loss:  {val_loss:.5f}"  if val_loss  else "  val_loss:  N/A")
            print(f"  test_loss: {test_loss:.5f}" if test_loss else "  test_loss: N/A")
            print("===================================\n")

    finally:
        wandb.finish()
        del model, datamodule, trainer
        torch.cuda.empty_cache()
        gc.collect()


if __name__ == "__main__":
    main()