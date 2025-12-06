import os
import json
import time
import pandas as pd
import numpy as np
from pathlib import Path
from datetime import datetime, timedelta
from tqdm import tqdm
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict, Counter
from multiprocessing import Pool, cpu_count
import logging
import warnings
from scipy import stats
warnings.filterwarnings('ignore')

from clinical_constants import *

def load_patient_files(cleaned_dir: Path, sample_size: int = None, logger=None) -> list:
    """
    Load all patient trajectory JSON files.
    Optionally subsample for speed.
    """
    patient_files = sorted(cleaned_dir.glob("*.json"))
    
    if sample_size is not None and sample_size < len(patient_files):
        np.random.seed(123)
        patient_files = list(np.random.choice(patient_files, size=sample_size, replace=False))
    if logger is not None:
        logger.info(f"Loading {len(patient_files):,} patient files")
    return patient_files

def identify_sepsis_bool(patient_data: dict) -> bool:
    """Quick sepsis identification."""
    records = patient_data.get('records', [])
    if not isinstance(records, list):
        return False
    
    for record in records:
        if not isinstance(record, dict):
            continue
        if record.get('data_type') == 'diagnoses':
            icd_code = record.get('icd_code', '')
            icd_version = record.get('icd_version', 10)
            
            if icd_version == 10:
                for sepsis_code in SEPSIS_ICD10_CODES:
                    if icd_code.startswith(sepsis_code):
                        return True
            elif icd_version == 9:
                for sepsis_code in SEPSIS_ICD9_CODES:
                    if icd_code.startswith(sepsis_code):
                        return True
    
    return False

def identify_sepsis_detailed(patient_data: dict) -> dict:
    """
    Identify sepsis in patient record using multiple criteria.
    
    Returns:
        Dictionary with sepsis flags and details:
        - has_sepsis: Boolean
        - sepsis_type: 'explicit', 'suspected', or None
        - sepsis_admissions: List of hadm_ids with sepsis
        - first_sepsis_time: Earliest sepsis diagnosis
        - sepsis_codes: ICD codes found
    """
    result = {
        'has_sepsis': False,
        'sepsis_type': None,
        'sepsis_admissions': [],
        'first_sepsis_time': None,
        'sepsis_codes': [],
        'organ_dysfunctions': [],
    }
    
    # Check diagnoses for explicit sepsis codes
    diagnoses = patient_data.get('records', [])
    if not isinstance(diagnoses, list):
        diagnoses = []
    
    sepsis_diagnoses = []
    sepsis_hadm_ids = set()  # Use set to avoid duplicates
    
    for record in diagnoses:
        if not isinstance(record, dict):
            continue
            
        if record.get('data_type') == 'diagnoses':
            icd_code = record.get('icd_code', '')
            if not icd_code:
                continue
                
            icd_version = record.get('icd_version', 10)
            
            # Check ICD-10
            if icd_version == 10:
                for sepsis_code in SEPSIS_ICD10_CODES:
                    if icd_code.startswith(sepsis_code):
                        sepsis_diagnoses.append({
                            'code': icd_code,
                            'hadm_id': record.get('hadm_id'),
                            'version': 10
                        })
                        result['sepsis_codes'].append(icd_code)
                        if record.get('hadm_id'):
                            sepsis_hadm_ids.add(record['hadm_id'])
                        break
            
            # Check ICD-9
            elif icd_version == 9:
                for sepsis_code in SEPSIS_ICD9_CODES:
                    if icd_code.startswith(sepsis_code):
                        sepsis_diagnoses.append({
                            'code': icd_code,
                            'hadm_id': record.get('hadm_id'),
                            'version': 9
                        })
                        result['sepsis_codes'].append(icd_code)
                        if record.get('hadm_id'):
                            sepsis_hadm_ids.add(record['hadm_id'])
                        break
            
            # Check for organ dysfunction
            for organ, codes in ORGAN_DYSFUNCTION_ICD10.items():
                for dysfunction_code in codes:
                    if icd_code.startswith(dysfunction_code):
                        result['organ_dysfunctions'].append(organ)
                        break
    
    if sepsis_diagnoses:
        result['has_sepsis'] = True
        result['sepsis_type'] = 'explicit'
        result['sepsis_admissions'] = list(sepsis_hadm_ids)
    
    # Check for suspected sepsis (infection + organ dysfunction)
    if not result['has_sepsis'] and len(set(result['organ_dysfunctions'])) >= 2:
        # Check for infection indicators
        has_infection = False
        for record in diagnoses:
            if not isinstance(record, dict):
                continue
            if record.get('data_type') == 'diagnoses':
                icd_code = record.get('icd_code', '')
                # Infection codes (simplified)
                if any(icd_code.startswith(prefix) for prefix in ['A', 'B', 'J']):
                    has_infection = True
                    break
        
        if has_infection:
            result['has_sepsis'] = True
            result['sepsis_type'] = 'suspected'
    
    result['sepsis_codes'] = list(set(result['sepsis_codes']))
    result['organ_dysfunctions'] = list(set(result['organ_dysfunctions']))
    
    return result