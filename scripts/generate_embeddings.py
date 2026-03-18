import anndata as ad
import numpy as np
from pathlib import Path
import torch
from tqdm import tqdm

import sys

project_root = Path.cwd().parent
sys.path.append(str(project_root / "src"))

from MethylCDM.models.betaVAE import BetaVAE

DEVICE = "cpu"
CHECKPOINT_PATH = "/ddn_exa/campbell/sli/methylcdm-project/models/beta_vae/betaVAE_sweep_20260306_133219/trial_72/best-epoch=151-val_loss=1.3515.ckpt"
DATA_PATH = "/ddn_exa/campbell/sli/methylcdm-project/data/training/methylation/pancancer_cohort_adata.h5ad"
OUTPUT_BASE = Path("/ddn_exa/campbell/sli/methylcdm-project/data/embeddings")


def main():

    print("~~~~~| Beginning to generate embeddings |~~~~~")

    adata = ad.read_h5ad(DATA_PATH)
    print(f"Loaded AnnData: {adata.shape}")

    model = BetaVAE.load_from_checkpoint(CHECKPOINT_PATH, map_location=DEVICE)
    model.eval()
    model.to(DEVICE)
    print("Model loaded on CPU")

    # Load the AnnData cohort object
    for i, (idx, row) in enumerate(tqdm(adata.obs.iterrows(), total = adata.n_obs)):
        barcode = row['barcode']
        file_id = idx.upper()
        project = row['project_id']
    
        project_dir = OUTPUT_BASE / project
        project_dir.mkdir(exist_ok = True, parents = True)
    
        sample_data = adata[i, :].X
        if not isinstance(sample_data, np.ndarray):
            sample_data = sample_data.toarray() if hasattr(sample_data, "toarray") else np.array(sample_data)
        sample_tensor = torch.tensor(sample_data, dtype=torch.float32, device=DEVICE)
    
        embedding = generate_embeddings(model, sample_tensor)
        filename = f"{barcode}.{file_id}.npy"
        np.save(project_dir / filename, embedding)

        if i < 5:  # print only first few to avoid flooding
            print(f"Saved embedding: {project_dir / filename}, shape: {embedding.shape}")
    

    print("~~~~~| Finished generating embeddings |~~~~~")


def generate_embeddings(model, sample_tensor):
    """
    Computes beta-VAE latente embeddings on a single sample tensor.
    """
    with torch.no_grad():
        z_mu, _, _ = model.encode(sample_tensor)
    return z_mu.squeeze(0).cpu().numpy()

if __name__ == "__main__":
    main()

# [END]