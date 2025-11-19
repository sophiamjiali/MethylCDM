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

def download_methylation(config, verbose = False):
    """
    Downloads DNA methylation data (beta values, .txt) from TCGA from the 
    specified projects, saving them per project. A manifest of file IDs is
    created to use the `gdc-client` API, preserving metadata.

    Parameters
    ----------
    config : Configuration object containing:
        - projects (list): list of TCGA project names (string)
        - data_dir (string): absolute path to the output data directory
        - metadata_dir (string): absolute path to the output metadata directory
        - metadata_fields (list): list of metadata columns to fetch from GDC

    verbose (boolean): toggle for verbose processing 
    """

    # Fetch all relevant values from the configurations object
    projects = config.get('projects', [])
    data_dir = config.get('data_dir', '')
    metadata_dir = config.get('metadata_dir', '')
    metadata_fields = config.get('metadata_fields', [])

    # Initialize the output directories if necessary
    if data_dir: os.makedirs(data_dir, exist_ok = True)
    if metadata_dir: os.makedirs(metadata_dir, exist_ok = True)

    if verbose:
        print("=" * 50)
        print(f"Beginning to download {len(projects)} projects...")
        print("=" * 50)


    # Query the GDC API for DNA methylation beta values from the projects
    query_url = "https://api.gdc.cancer.gov/files"

    # Request and download data per project, saving data in individual folders
    for project in projects:

        if verbose: print(f"Beginning to download project {project}:")

        # -----| Query DNA Methylation Files |-----
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
        manifest_file = os.path.join(metadata_dir, 
                                     f"{project}_methylation_manifest.txt")
        
        with open(manifest_file, "w") as f:
            f.write("id\tfilename\n")
            for file in files:
                f.write(f"{file['file_id']}\t{file['file_name']}\n")
        
        if verbose: print(f"Generated manifest with {len(files)} files.")


        # -----| Download Data from GDC Client with Manifest |-----
        project_data_dir = os.path.join(data_dir, f"{project}_methylation_data")
        os.makedirs(project_data_dir, exist_ok = True)

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
        project_metadata_dir = os.path.join(
            metadata_dir, f"{project}_methylation_metadata"
        )
        os.makedirs(project_metadata_dir, exist_ok = True)

        metadata_file = os.path.join(
            project_metadata_dir, f"{project}_methylation_metadata.csv"
        )
        metadata.to_csv(metadata_file, index = False)

        if verbose: print(f"Downloaded project metadata.")


    if verbose:
        print("=" * 50)
        print(f"Completed downloading {len(projects)} projects.")
        print("=" * 50), 
