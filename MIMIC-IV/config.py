# config.py
from pathlib import Path
from multiprocessing import cpu_count
# from scratch import OUTPUT_DIR

# Versioning
CLEANED_VERSION = 0

# Directories
DATA_DIR = Path('/ewsc/jpickard/physionet/mimiciv/physionet.org/files/mimiciv/3.1')
CLEANED_DIR = DATA_DIR / f"patient_trajectories_cleaned_v{CLEANED_VERSION}"
ANALYSIS_DIR = Path("/ewsc/jpickard/EHR-DE/MIMIC-IV/EDA")
EDA_DIR = ANALYSIS_DIR / f"eda_analysis_v{CLEANED_VERSION}"

# Performance
NUM_WORKERS = min(64, cpu_count() - 1)
SAMPLE_SIZE = 10000 # None  # Set to integer to analyze subset
