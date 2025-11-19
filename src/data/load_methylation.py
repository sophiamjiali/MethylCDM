# ==============================================================================
# Script:           load_methylation.py
# Purpose:          Loads and downloads DNA methylation data for model input
# Author:           Sophia Li
# Affiliation:      CCG Lab, Princess Margaret Cancer Center, UHN, UofT
# Date:             11/13/2025
# ==============================================================================

import pandas as pd
import subprocess
import requests
import json
import os

# =====| Download Methylation Data |============================================

def download_methylation(project, config, verbose = False):
    """
    Downloads DNA methylation data (beta values, .txt) from TCGA from the 
    specified projects, saving them per project. A manifest of file IDs is
    created to use the `gdc-client` API, preserving metadata.

    Parameters
    ----------
    project (str): name of a TCGA project with data available on GDC

    config : Configuration object containing:
        - projects (list): list of TCGA project names (string)
        - raw_data_dir (str): path to the output raw data directory
        - metadata_dir (str): path to the output metadata directory
        - metadata_fields (list): list of metadata columns to fetch from GDC

    verbose (boolean): toggle for verbose processing 
    """

    # Fetch all relevant values from the configurations object
    download_cfg = config.get('download', {})
    raw_data_dir = download_cfg.get('raw_data_dir', '')
    metadata_dir = download_cfg.get('metadata_dir', '')
    metadata_fields = download_cfg.get('metadata_fields', [])

    # Initialize the output directories if necessary
    project_data_dir = os.path.join(
        raw_data_dir, f"{project}_raw_methylation"
    )
    if raw_data_dir: os.makedirs(project_data_dir, exist_ok = True)
    if metadata_dir: os.makedirs(metadata_dir, exist_ok = True)

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
            {"op": "in", "content": {"field": "files.data_category", 
                                        "value": ["DNA Methylation"]}},
            {"op": "in", "content": {"field": "files.data_type", 
                                        "value": ["Methylation Beta Value"]}},
            {"op": "in", "content": {"field": "files.access", 
                                        "value": ["open"]}},
        ]
    }
    data_query = {
        "filters": json.dumps(filters),
        "fields": "file_id,file_name",
        "format": "JSON"
    }

    # Query the GDC API for methylation files, filtered by configurations
    response = requests.post(query_url, json = data_query)
    data = response.json()
    files = data['data']['hits']


    # -----| Generate DNA Methylation Data Manifest |-----
    manifest_file = os.path.join(
        metadata_dir, f"{project}_methylation_manifest.txt"
    )
    
    with open(manifest_file, "w") as f:
        f.write("id\tfilename\n")
        for file in files:
            f.write(f"{file['file_id']}\t{file['file_name']}\n")
    
    if verbose: print(f"Generated manifest with {len(files)} files.")


    # -----| Download Data from GDC Client with Manifest |-----

    # Spawn a sub-process to call the `gdc-client` with the manifest
    subprocess.run(["gdc-client", "download", "-m", 
                    manifest_file, "-d", project_data_dir], 
                    check = True)
    
    if verbose: print(f"Downloaded {len(files)} files.\n")


    # -----| Download Metadata |-----

    # Query for the metadata fields provided in the configurations
    metadata_query = {
        "filters": json.dumps(filters),
        "fields": ",".join(metadata_fields),
        "format": "JSON"
    }

    response = requests.post(query_url, json = metadata_query)
    data = response.json()
    metadata = pd.json_normalize(data['data']['hits'])

    # Normalize into a table and save it
    metadata_file = os.path.join(
        metadata_dir, f"{project}_methylation_metadata.csv"
    )
    metadata.to_csv(metadata_file, index = False)

    if verbose: print(f"Downloaded project metadata.")


    if verbose:
        print("=" * 50)
        print(f"Completed downloading TCGA project {project}.")
        print("=" * 50), 


# =====| Load Methylation Data |============================================

def load_raw_methylation(project, config, verbose = False):
    """
    Loads DNA methylation data (beta values) as a Pandas DataFrame after being 
    downloaded using `download_methylation()`.

    Parameters
    ----------
    project (str): name of a TCGA project with data available on GDC

    config : Configuration object containing:
        - 

    verbose (boolean): toggle for verbose processing 

    Returns
    -------
    methylation (list): list of (probes x samples) DataFrames per project
    """

    # load each file contained in the project raw data sub-directory

    return


def load_methylation():
    """
    Loads processed DNA methylation data (beta values) as a Pandas DataFrame.

    Parameters
    ----------
    ...

    Returns
    -------
    ...
    """

    # Prior to gene group definition, just preprocessing of beta values per
    # project