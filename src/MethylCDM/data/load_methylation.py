# ==============================================================================
# Script:           load_methylation.py
# Purpose:          Loads and downloads DNA methylation data for model input
# Author:           Sophia Li
# Affiliation:      CCG Lab, Princess Margaret Cancer Center, UHN, UofT
# Date:             11/13/2025
# ==============================================================================

from pathlib import Path
import pandas as pd
import subprocess
import time
import requests
import json
import os
import shutil
import tempfile

from MethylCDM.utils.utils import resolve_path
from MethylCDM.constants import (
    GDC_CLIENT_PATH, 
    RAW_METHYLATION_DIR,
    METADATA_METHYLATION_DIR,
    METADATA_METHYLATION,
    CHUNK_SIZE,
    MAX_RETRIES,
    RETRY_SLEEP
)

# =====| Download Methylation Data |============================================

def download_methylation(project, config, verbose = False):
    """
    Downloads DNA methylation data (beta values, .txt) from TCGA from the 
    specified projects, saving them per project. A manifest of file IDs is
    created to use the `gdc-client` API, preserving metadata.

    A retry loop will re-run the download command if it fails. The manifest is temporarily split into chunks to break up long API requests (which is the 
    source of main downloading issues).

    Parameters
    ----------
    project (str): name of a TCGA project with data available on GDC

    config (dict): a configuration object containing:
        - projects (list): list of TCGA project names (string)
        - raw_data_dir (str): path to the output raw data directory
        - metadata_dir (str): path to the output metadata directory

    verbose (boolean): toggle for verbose processing 
    """

    # Fetch all relevant values from the configurations object
    download_cfg = config.get('download', {})
    raw_data_dir = download_cfg.get('raw_data_dir', '')
    metadata_dir = download_cfg.get('metadata_dir', '')

    # Resolve directory paths relative to the project root
    raw_data_dir = resolve_path(raw_data_dir, RAW_METHYLATION_DIR)
    metadata_dir = resolve_path(metadata_dir, METADATA_METHYLATION_DIR)

    # Initialize the output directories if necessary
    project_data_dir = os.path.join(raw_data_dir, project)
    project_metadata_dir = os.path.join(metadata_dir, project)

    if raw_data_dir: os.makedirs(project_data_dir, exist_ok = True)
    if metadata_dir: os.makedirs(project_metadata_dir, exist_ok = True)

    if verbose:
        print("=" * 50)
        print(f"Beginning to download TCGA project {project}")
        print("=" * 50)


    # -----| Query DNA Methylation Files |-----

    # Query the GDC API for DNA methylation beta values from the projects
    query_url = "https://api.gdc.cancer.gov/files"

    filters = {
        "op": "and",
        "content": [
            {"op": "in", "content": {"field": "cases.project.project_id", 
                                     "value": [project]}},
            {"op": "in", "content": {"field": "files.experimental_strategy", 
                                     "value": "Methylation Array"}},
            {"op": "in", "content": {"field": "files.data_category", 
                                     "value": ["dna methylation"]}},
            {"op": "in", "content": {"field": "files.data_type", 
                                     "value": ["Methylation Beta Value"]}},
            {"op": "in", "content": {"field": "files.access", 
                                     "value": ["open"]}},
        ]
    }
    data_query = {
        "filters": json.dumps(filters),
        "fields": "file_id,file_name",
        "format": "JSON",
        "size": 10000
    }

    # Query the GDC API for methylation files, filtered by configurations
    response = requests.post(query_url, json = data_query)
    data = response.json()
    files = data['data']['hits']


    # -----| Generate DNA Methylation Data Manifest |-----
    manifest_file = os.path.join(project_metadata_dir, 
                                 f"{project}_manifest.txt")
    
    with open(manifest_file, "w") as f:
        f.write("id\tfilename\n")
        for file in files:
            f.write(f"{file['file_id']}\t{file['file_name']}\n")
    
    if verbose: print(f"Generated manifest with {len(files)} files.")


    # -----| Download Data from GDC Client with Manifest |-----

    # Require the GDC-client tool to be downloaded 
    if not GDC_CLIENT_PATH.exists():
        raise FileNotFoundError("`gdc-client` was not found. Please download"
                                " and place it in the `tools/` directory.")
    
    # Break the manifest into subprocesses to reduce long API requests
    with Path(manifest_file).open() as f:
        lines = f.readlines()

    header = lines[0]
    entries = lines[1:]

    # Create a temporary directory for chunked manifests
    with tempfile.TemporaryDirectory() as tmpdir:

        # Split into chunks
        tmpdir = Path(tmpdir)
        manifest_chunks = []
        for i in range(0, len(entries), CHUNK_SIZE):
            chunk_path = tmpdir / f"manifest_chunk_{i // CHUNK_SIZE}.txt"
            with chunk_path.open("w") as f:
                f.write(header)
                f.writelines(entries[i : i + CHUNK_SIZE])
            manifest_chunks.append(chunk_path)

        if verbose: print(f"Split manifest into {len(manifest_chunks)} chunks.")

        # Download each chunk with retries to handle API failures
        chunk_idx = 0
        for chunk in manifest_chunks:
            chunk_idx += 1
            for attempt in range(1, MAX_RETRIES + 1):
                try:
                    # Spawn a sub-process to call the `gdc-client`
                    subprocess.run([GDC_CLIENT_PATH, "download", "-m", 
                                    str(chunk), "-d", project_data_dir], 
                                    check = True,
                                    capture_output = True,
                                    text = True)
                    if verbose: print(f"Downloaded {CHUNK_SIZE} files.\n")
                    break

                except subprocess.CalledProcessError as e:
                    if verbose:
                        print(f"Attempt {attempt} failed.")
                        print("STDOUT:\n", e.stdout)
                        print("STDERR:\n", e.stderr)
                    
                    # Re-attempt if within five tries, else throw final error
                    if attempt < MAX_RETRIES:
                        print(f"Retrying in {RETRY_SLEEP} seconds...")
                        time.sleep(15)
                    else:
                        raise RuntimeError(
                            f"GDC-client failed after 5 attempts"
                        ) from e
                    
            if verbose: print(f"Processed chunk {chunk_idx} out of "
                              f"{len(manifest_chunks)}")

    # -----| Download Metadata |-----

    # Query for the metadata fields provided in the configurations
    metadata_query = {
        "filters": json.dumps(filters),
        "fields": ",".join(METADATA_METHYLATION),
        "format": "JSON",
        "size": 10000
    }

    response = requests.post(query_url, json = metadata_query)
    data = response.json()
    metadata = pd.json_normalize(data['data']['hits'])


    # -----| Normalize Metadata |-----
    metadata = metadata.explode('cases', ignore_index = True)
    cases_metadata = pd.json_normalize(metadata['cases'])

    metadata = pd.concat([
        metadata.drop(columns = ["cases"]).reset_index(drop = True),
        cases_metadata.reset_index(drop = True)
    ], axis = 1)

    # Expand the diagnoses and samples fields nested in the cases
    metadata = metadata.explode('diagnoses').explode('samples')
    diagnoses_metadata = pd.json_normalize(metadata['diagnoses'])
    samples_metadata = pd.json_normalize(metadata['samples'])

    metadata = pd.concat([
        metadata.drop(columns=['diagnoses', 'samples']).reset_index(drop=True),
        diagnoses_metadata.reset_index(drop = True),
        samples_metadata.reset_index(drop = True)
    ], axis = 1)

    # Clean any remaining prefixed column names
    metadata.columns = [col.split('.')[-1] for col in metadata.columns]

    # Save the metadata to the designated data folder
    metadata_file = os.path.join(project_metadata_dir,
                                 f"{project}_metadata.csv")
    metadata.to_csv(metadata_file, index = False)

    if verbose: print(f"Downloaded project metadata.")

    if verbose:
        print("=" * 50)
        print(f"Completed downloading TCGA project {project}.")
        print("=" * 50)

# =====| Clean Methylation Data |===============================================

def clean_methylation_data(dir_path, verbose = False):
    """
    Cleans all beta value files downloaded using the `gdc-client`, converting
    *.txt to *.parquet and removing accessory files.

    Parameters
    ----------
    dir_path (str): path to the directory containing nested beta value files
    verbose (boolean): toggle for verbose processing 
    """

    if verbose:
        print("=" * 50)
        print(f"Beginning to clean methylation data")
        print("=" * 50)

    # Recursively find all beta value *.txt files
    txt_files = list(Path(dir_path).rglob("*.level3betas.txt"))

    if not txt_files:
        print("All files are already cleaned.")
        
    else:
        for txt_path in txt_files:

            # Build the parquet filename to match the original
            parquet_path = dir_path / Path(txt_path.stem + ".parquet")

            # Convert the .txt to .parquet
            txt = pd.read_csv(txt_path, sep = "\t", index_col = 0)
            txt.index.name = "probe_id"
            txt.columns = ["beta_value"]
            txt = txt.astype("float32")
            txt.to_parquet(parquet_path, index = True)

            # Delete the old .txt file
            txt_path.unlink()

            print(f"Successfully cleaned file {txt_path.stem}")


        # Delete all nested directories and files
        for subdir in Path(dir_path).iterdir():
            if subdir.is_dir(): shutil.rmtree(subdir)

    if verbose:
        print("=" * 50)
        print(f"Completed cleaning methylation data")
        print("=" * 50)


# =====| Load Methylation Data |================================================

# def merge_cohort(project, config):
#     """
#     Loads individual DNA methylation data (beta values) and merges them into
#     a cohort-level Pandas DataFrame, returning the matrix. The matrix is 
#     filtered to only contain data from the common set of probes between
#     all samples (i.e. common set between Illumina manifests).

#     Parameters
#     ----------
#     project (str): name of a TCGA project with data available on GDC
#     config (dict): Configuration object containing:
#         - raw_data_dir (str): path to the output raw data directory

#     Returns
#     -------
#     beta_values (DataFrame): raw cohort-level CpG x Sample ID matrix of
#                              beta values of a given project

#     Raises
#     ------
#     FileNotFoundError: if the raw data directory does not exist or is empty at 
#                        the specified path
#     """

#     # Resolve the project's raw data directory
#     raw_data_dir = config.get('download', {}).get('raw_data_dir', '')
#     raw_data_dir = resolve_path(raw_data_dir, RAW_METHYLATION_DIR)
#     project_data_dir = os.path.join(raw_data_dir, f"{project}")

#     # Verify the raw data exists and is not empty
#     if not os.path.isdir(project_data_dir):
#         raise FileNotFoundError(f"Raw data directory was not found at "
#                                 f"{project_data_dir}.")
#     if not os.listdir(project_data_dir):
#         raise FileNotFoundError(f"Raw data directory was empty at "
#                                 f"{project_data_dir}.")
    
#     # Identify and load all nested beta value .txt files
#     beta_files = [f for f in Path(project_data_dir).glob("*/*.level3betas.txt")]
#     beta_values = [load_beta_file(f) for f in beta_files]

#     # Identify the common probes between present manifests and filter
#     common_probes = beta_values[0].index
#     for df in beta_values[1:]:
#         common_probes = common_probes.intersection(df.index)
#     beta_values = [df.loc[common_probes].astype("float32") 
#                    for df in beta_values]

#     # Merge beta values into a single matrix, keeping only common probes
#     beta_values = pd.concat(beta_values, axis = 1)
    
#     # Sort the CpGs
#     beta_values = beta_values.sort_index()

#     return beta_values


# def load_methylation():
#     """
#     Loads processed DNA methylation data (beta values) as a Pandas DataFrame.

#     Parameters
#     ----------
#     ...

#     Returns
#     -------
#     ...
#     """

#     # Prior to gene group definition, just preprocessing of beta values per
#     # project