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

# Compute the GDC-Client path
GDC_CLIENT_PATH = PROJECT_ROOT / "tools/gdc-client"

# =====| Data Directory |=======================================================

DATA_DIR = PROJECT_ROOT / "data"

RAW_METHYLATION_DIR = DATA_DIR / "raw" / "methylation"
INTERMEDIATE_METHYLATION_DIR = DATA_DIR / "intermediate" / "methylation"
PROCESSED_METHYLATION_DIR = DATA_DIR / "processed" / "methylation"
METADATA_METHYLATION_DIR = DATA_DIR / "metadata" / "methylation"

RAW_WSI_DIR = DATA_DIR / "raw" / "wsi"
INTERMEDIATE_WSI_DIR = DATA_DIR / "intermediate" / "wsi"
PROCESSED_WSI_DIR = DATA_DIR / "processed" / "wsi"
METADATA_WSI_DIR = DATA_DIR / "metadata" / "wsi"

# =====| Probe Annotations |====================================================

RESOURCES_DIR = PROJECT_ROOT / "resources"
ANNOTATION_27K = RESOURCES_DIR / "IlluminaHumanMethylation450.csv"
ANNOTATION_450K = RESOURCES_DIR / "IlluminaHumanMethylation27.csv"
ANNOTATION_EPIC = None # empty for now

# =====| GDC Metadata |=========================================================

# Define DNA Methylation metadata fields
METADATA_METHYLATION = [

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