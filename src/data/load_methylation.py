# ==============================================================================
# Script:           load_methylation.py
# Purpose:          Loads and downloads DNA methylation data for model input
# Author:           Sophia Li
# Affiliation:      CCG Lab, Princess Margaret Cancer Center, UHN, UofT
# Date:             11/13/2025
# ==============================================================================

import requests
from io import StringIO
import pandas as pd
import json

def download_methylation(config):
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
    """

    # Request and download data per project, saving data in individual folders
    for project in config.get('projects', []):

        # -----| Download DNA Methylation Data |-----

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
                                         "value": ["Methylation Beta Values"]}},
                {"op": "in", "content": {"field": "files.access", 
                                         "value": ["open"]}},
            ]
        }
        params = {
            "filters": json.dumps(filters),
            "fields": "file_id",
            "format": "JSON"
        }

        response = requests.get(query_url, params = params)
        content = StringIO(str(response.content, 'utf-8'))

        file_ids = pd.read_csv(content)



        # -----| Download Metadata |-----
