# ==============================================================================
# Script:           constants.py
# Purpose:          Defines global constants used across the workflow
# Author:           Sophia Li
# Affiliation:      CCG Lab, Princess Margaret Cancer Center, UHN, UofT
# Date:             11/22/2025
# ==============================================================================

from pathlib import Path

# =====| Paths |================================================================

# Compute project path relative to this file
PROJECT_ROOT = Path(__file__).parent.parent.parent

# Compute default paths
CONFIG_DIR = PROJECT_ROOT / "config"
DATA_DIR = PROJECT_ROOT / "data"

# Compute the GDC-Client path
GDC_CLIENT_PATH = PROJECT_ROOT / "tools/gdc-client"


# =====| GDC Metadata |=========================================================

# Define DNA Methylation metadata fields
METHYLATION_METADATA = [

    # File-Level Metadata
    "project.project_id",
    "file_id",
    "file_name",
    "data_category",
    "data_type",
    "experimental_strategy",
    "platform",
    "state",
    "data_format",
    "cases.case_id",
    "cases.submitter_id",
    "cases.samples.sample_id",
    "cases.samples.sample_type",

    # Clinical Metadata
    "cases.diagnoses.age_at_diagnosis",
    "cases.demographic.gender",
    "cases.demographic.race",
    "cases.demographic.ethnicity",
    "cases.diagnoses.primary_diagnosis",
    "cases.diagnoses.morphology",
    "cases.diagnoses.tumor_stage",
    "cases.diagnoses.tumor_grade",
    "cases.diagnoses.days_to_death",
    "cases.diagnoses.days_to_last_follow_up",
    "cases.diagnoses.vital_status",
    "cases.prior_malignancy",
    "cases.sites_of_resection_or_biopsy"
]