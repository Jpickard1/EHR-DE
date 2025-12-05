#!/usr/bin/env python3
"""
MIMIC-IV Patient Trajectory Data Cleaner
=========================================
Converts raw JSONL patient trajectories to cleaned, time-sorted JSON format
suitable for conversion to pandas DataFrames.
"""

import os
import json
import pandas as pd
from pathlib import Path
from datetime import datetime
from tqdm import tqdm
import time
from multiprocessing import Pool, cpu_count
from collections import defaultdict
import logging

# ============================================================================
# CONFIGURATION
# ============================================================================

# Version control
TRAJECTORY_VERSION = 0
CLEANED_VERSION = 0

# Paths
DATA_DIR = Path('/ewsc/jpickard/physionet/mimiciv/physionet.org/files/mimiciv/3.1')
INPUT_DIR = DATA_DIR / f"patient_trajectories_v{TRAJECTORY_VERSION}"
OUTPUT_DIR = DATA_DIR / f"patient_trajectories_cleaned_v{CLEANED_VERSION}"
PROGRESS_DIR = OUTPUT_DIR / 'progress'
LOG_DIR = OUTPUT_DIR / 'logs'

# Performance settings
MAX_CPUS = min(20, cpu_count() - 1)
NUM_WORKERS = max(1, MAX_CPUS)
BATCH_SIZE = 100  # Number of files to process before saving progress

# Create directories
OUTPUT_DIR.mkdir(exist_ok=True)
PROGRESS_DIR.mkdir(exist_ok=True)
LOG_DIR.mkdir(exist_ok=True)

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler(LOG_DIR / f'cleaning_{datetime.now().strftime("%Y%m%d_%H%M%S")}.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# ============================================================================
# TIME EXTRACTION UTILITIES
# ============================================================================

# Define time column mappings for different data types
TIME_COLUMN_MAP = {
    'admissions': 'admittime',
    'icu_stays': 'intime',
    'events': 'charttime',  # Primary time field
    'prescriptions': 'starttime',
    'diagnoses': None,  # Untimed, links to admission
    'procedures': None,  # Untimed, links to admission
    'microbiology': 'charttime',
}

# Alternative time columns (fallback if primary is missing)
ALTERNATIVE_TIME_COLUMNS = [
    'charttime', 'starttime', 'endtime', 'storetime', 
    'admittime', 'dischtime', 'intime', 'outtime',
    'edregtime', 'edouttime'
]

def extract_timestamp(record: dict, data_type: str) -> pd.Timestamp:
    """
    Extract timestamp from a record with intelligent fallback.
    
    Args:
        record: Dictionary containing event data
        data_type: Type of data (e.g., 'events', 'admissions')
        
    Returns:
        pandas Timestamp or None if no valid timestamp found
    """
    # Try primary time column for this data type
    primary_col = TIME_COLUMN_MAP.get(data_type)
    if primary_col and primary_col in record:
        try:
            return pd.to_datetime(record[primary_col])
        except:
            pass
    
    # Try alternative time columns
    for col in ALTERNATIVE_TIME_COLUMNS:
        if col in record and record[col] is not None:
            try:
                return pd.to_datetime(record[col])
            except:
                continue
    
    return None

def extract_hadm_id(record: dict) -> int:
    """Extract hospital admission ID with None as default."""
    hadm_id = record.get('hadm_id')
    return int(hadm_id) if hadm_id is not None and pd.notna(hadm_id) else None

def extract_icustay_id(record: dict) -> int:
    """Extract ICU stay ID with None as default."""
    # Try different possible column names
    for col in ['stay_id', 'icustay_id']:
        if col in record:
            stay_id = record.get(col)
            return int(stay_id) if stay_id is not None and pd.notna(stay_id) else None
    return None

# ============================================================================
# PROGRESS TRACKING
# ============================================================================

class CleaningProgressTracker:
    """Track cleaning progress with checkpoint capability."""
    
    def __init__(self, progress_file: Path):
        self.progress_file = progress_file
        self.completed_patients = set()
        self.failed_patients = {}
        self.stats = {
            'total_processed': 0,
            'successful': 0,
            'failed': 0,
            'start_time': None,
            'last_update': None
        }
        self.load()
    
    def load(self):
        """Load progress from checkpoint file."""
        if self.progress_file.exists():
            try:
                with open(self.progress_file, 'r') as f:
                    data = json.load(f)
                    self.completed_patients = set(data.get('completed_patients', []))
                    self.failed_patients = data.get('failed_patients', {})
                    self.stats = data.get('stats', self.stats)
                logger.info(f"Loaded progress: {len(self.completed_patients)} patients already processed")
            except Exception as e:
                logger.warning(f"Could not load progress file: {e}")
    
    def save(self):
        """Save current progress to checkpoint file."""
        try:
            with open(self.progress_file, 'w') as f:
                json.dump({
                    'completed_patients': list(self.completed_patients),
                    'failed_patients': self.failed_patients,
                    'stats': self.stats
                }, f, indent=2)
        except Exception as e:
            logger.error(f"Failed to save progress: {e}")
    
    def mark_completed(self, patient_id: int):
        """Mark a patient as successfully processed."""
        self.completed_patients.add(patient_id)
        self.stats['successful'] += 1
        self.stats['total_processed'] += 1
        self.stats['last_update'] = datetime.now().isoformat()
    
    def mark_failed(self, patient_id: int, error: str):
        """Mark a patient as failed with error message."""
        self.failed_patients[str(patient_id)] = {
            'error': str(error),
            'timestamp': datetime.now().isoformat()
        }
        self.stats['failed'] += 1
        self.stats['total_processed'] += 1
        self.stats['last_update'] = datetime.now().isoformat()
    
    def is_completed(self, patient_id: int) -> bool:
        """Check if patient has been successfully processed."""
        return patient_id in self.completed_patients
    
    def get_summary(self) -> dict:
        """Get summary statistics."""
        return {
            **self.stats,
            'completion_rate': (self.stats['successful'] / self.stats['total_processed'] * 100 
                              if self.stats['total_processed'] > 0 else 0)
        }

# ============================================================================
# DATA CLEANING FUNCTIONS
# ============================================================================

def flatten_record(record: dict, data_type: str) -> dict:
    """
    Flatten a single record into a format ready for DataFrame conversion.
    
    Args:
        record: Dictionary containing the record data
        data_type: Type of data (demographics, events, etc.)
        
    Returns:
        Flattened dictionary with standardized keys
    """
    flattened = {
        'data_type': data_type,
        'time': None,
        'hadm_id': None,
        'stay_id': None,
    }
    
    # Extract structured fields
    timestamp = extract_timestamp(record, data_type)
    if timestamp is not None:
        flattened['time'] = timestamp.isoformat()
    
    flattened['hadm_id'] = extract_hadm_id(record)
    flattened['stay_id'] = extract_icustay_id(record)
    
    # Add all other fields from the record
    for key, value in record.items():
        if key not in flattened:
            # Convert non-serializable types
            if pd.isna(value):
                flattened[key] = None
            elif isinstance(value, (pd.Timestamp, datetime)):
                flattened[key] = value.isoformat()
            else:
                flattened[key] = value
    
    return flattened

def clean_patient_file(patient_id: int, input_dir: Path, output_dir: Path) -> dict:
    """
    Clean a single patient's JSONL file and save as sorted JSON.
    
    Args:
        patient_id: Patient subject_id
        input_dir: Directory containing raw JSONL files
        output_dir: Directory for cleaned JSON output
        
    Returns:
        Dictionary with processing statistics
    """
    jsonl_file = input_dir / f'patient_{patient_id}.jsonl'
    json_file = output_dir / f'patient_{patient_id}.json'
    
    if not jsonl_file.exists():
        raise FileNotFoundError(f"Input file not found: {jsonl_file}")
    
    # Collect all records
    timed_records = []
    untimed_records = []
    demographics = {}
    
    stats = {
        'total_records': 0,
        'timed_records': 0,
        'untimed_records': 0,
        'data_types': defaultdict(int)
    }
    
    # Read and process JSONL
    with open(jsonl_file, 'r') as f:
        for line_num, line in enumerate(f, 1):
            if not line.strip():
                continue
                
            try:
                record = json.loads(line)
                data_type = record.get('data_type', 'unknown')
                data = record.get('data', {})
                
                stats['total_records'] += 1
                stats['data_types'][data_type] += 1
                
                # Handle demographics separately (single record)
                if data_type == 'demographics':
                    demographics = data
                    continue
                
                # Flatten and categorize record
                flattened = flatten_record(data, data_type)
                
                if flattened['time'] is not None:
                    timed_records.append(flattened)
                    stats['timed_records'] += 1
                else:
                    untimed_records.append(flattened)
                    stats['untimed_records'] += 1
                    
            except json.JSONDecodeError as e:
                logger.warning(f"Patient {patient_id}: JSON decode error on line {line_num}: {e}")
                continue
            except Exception as e:
                logger.warning(f"Patient {patient_id}: Error processing line {line_num}: {e}")
                continue
    
    # Sort timed records by timestamp
    timed_records.sort(key=lambda x: x['time'])
    
    # Combine: untimed first, then timed
    all_records = untimed_records + timed_records
    
    # Create final structure
    cleaned_data = {
        'subject_id': int(patient_id),
        'demographics': demographics,
        'records': all_records,
        'metadata': {
            'total_records': len(all_records),
            'timed_records': len(timed_records),
            'untimed_records': len(untimed_records),
            'data_type_counts': dict(stats['data_types']),
            'time_range': {
                'earliest': timed_records[0]['time'] if timed_records else None,
                'latest': timed_records[-1]['time'] if timed_records else None
            },
            'processed_at': datetime.now().isoformat(),
            'source_file': str(jsonl_file)
        }
    }
    
    # Write cleaned JSON
    with open(json_file, 'w') as f:
        json.dump(cleaned_data, f, indent=2, default=str)
    
    return stats

def clean_patient_wrapper(args):
    """Wrapper for parallel processing."""
    patient_id, input_dir, output_dir, progress_file = args
    
    try:
        stats = clean_patient_file(patient_id, input_dir, output_dir)
        return {
            'patient_id': patient_id,
            'status': 'success',
            'stats': stats
        }
    except Exception as e:
        logger.error(f"Failed to process patient {patient_id}: {e}")
        return {
            'patient_id': patient_id,
            'status': 'failed',
            'error': str(e)
        }

# ============================================================================
# MAIN PROCESSING PIPELINE
# ============================================================================

def main():
    """Main execution function."""
    
    print("=" * 80)
    print("MIMIC-IV PATIENT TRAJECTORY CLEANER")
    print("=" * 80)
    print(f"Start time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"\nInput directory:  {INPUT_DIR}")
    print(f"Output directory: {OUTPUT_DIR}")
    print(f"Parallel workers: {NUM_WORKERS}")
    print("=" * 80)
    print()
    
    start_time = time.time()
    
    # Initialize progress tracker
    progress_file = PROGRESS_DIR / 'cleaning_progress.json'
    tracker = CleaningProgressTracker(progress_file)
    tracker.stats['start_time'] = datetime.now().isoformat()
    
    # Find all patient files
    print("Scanning for patient files...")
    all_patient_files = sorted(INPUT_DIR.glob('patient_*.jsonl'))
    patient_ids = []
    
    for f in all_patient_files:
        try:
            patient_id = int(f.stem.split('_')[1])
            patient_ids.append(patient_id)
        except (ValueError, IndexError):
            logger.warning(f"Could not parse patient ID from {f.name}")
            continue
    
    print(f"Found {len(patient_ids):,} patient files")
    
    # Filter to unprocessed patients
    patients_to_process = [
        pid for pid in patient_ids 
        if not tracker.is_completed(pid)
    ]
    
    print(f"Already processed: {len(tracker.completed_patients):,}")
    print(f"Remaining: {len(patients_to_process):,}")
    print()
    
    if not patients_to_process:
        print("✓ All patients already processed!")
        return
    
    # Process patients
    print(f"Processing {len(patients_to_process):,} patients with {NUM_WORKERS} workers...")
    print()
    
    if NUM_WORKERS > 1:
        # Parallel processing
        args_list = [
            (pid, INPUT_DIR, OUTPUT_DIR, progress_file)
            for pid in patients_to_process
        ]
        
        with Pool(NUM_WORKERS) as pool:
            results_iter = pool.imap_unordered(clean_patient_wrapper, args_list, chunksize=10)
            
            batch_results = []
            for result in tqdm(results_iter, total=len(patients_to_process), desc="Cleaning"):
                batch_results.append(result)
                
                # Update progress
                if result['status'] == 'success':
                    tracker.mark_completed(result['patient_id'])
                else:
                    tracker.mark_failed(result['patient_id'], result.get('error', 'Unknown error'))
                
                # Save progress every BATCH_SIZE patients
                if len(batch_results) >= BATCH_SIZE:
                    tracker.save()
                    batch_results = []
    else:
        # Sequential processing (for debugging)
        for patient_id in tqdm(patients_to_process, desc="Cleaning"):
            result = clean_patient_wrapper((patient_id, INPUT_DIR, OUTPUT_DIR, progress_file))
            
            if result['status'] == 'success':
                tracker.mark_completed(result['patient_id'])
            else:
                tracker.mark_failed(result['patient_id'], result.get('error', 'Unknown error'))
            
            if tracker.stats['total_processed'] % BATCH_SIZE == 0:
                tracker.save()
    
    # Final save
    tracker.save()
    
    # Summary
    elapsed = time.time() - start_time
    summary = tracker.get_summary()
    
    print()
    print("=" * 80)
    print("CLEANING COMPLETE")
    print("=" * 80)
    print(f"Total processed:    {summary['total_processed']:,}")
    print(f"Successful:         {summary['successful']:,}")
    print(f"Failed:             {summary['failed']:,}")
    print(f"Completion rate:    {summary['completion_rate']:.2f}%")
    print(f"Elapsed time:       {elapsed/60:.1f} minutes")
    print(f"Processing rate:    {summary['total_processed']/(elapsed/60):.1f} patients/min")
    print()
    print(f"Output saved to: {OUTPUT_DIR}")
    print(f"Logs saved to:   {LOG_DIR}")
    print("=" * 80)
    
    # Report failures
    if tracker.failed_patients:
        print()
        print(f"WARNING: {len(tracker.failed_patients)} patients failed to process.")
        print(f"See log file for details.")
        failures_file = LOG_DIR / 'failed_patients.json'
        with open(failures_file, 'w') as f:
            json.dump(tracker.failed_patients, f, indent=2)
        print(f"Failed patient details saved to: {failures_file}")

# ============================================================================
# UTILITY: Convert cleaned JSON to DataFrame
# ============================================================================

def load_patient_dataframe(patient_id: int, cleaned_dir: Path = OUTPUT_DIR) -> pd.DataFrame:
    """
    Load a cleaned patient JSON file as a pandas DataFrame.
    
    Args:
        patient_id: Patient subject_id
        cleaned_dir: Directory containing cleaned JSON files
        
    Returns:
        pandas DataFrame with guaranteed columns: time, hadm_id, stay_id
        
    Example:
        >>> df = load_patient_dataframe(10000032)
        >>> print(df.columns)
        Index(['time', 'hadm_id', 'stay_id', 'data_type', ...])
    """
    json_file = cleaned_dir / f'patient_{patient_id}.json'
    
    with open(json_file, 'r') as f:
        data = json.load(f)
    
    # Convert records to DataFrame
    df = pd.DataFrame(data['records'])
    
    # Ensure guaranteed columns exist
    for col in ['time', 'hadm_id', 'stay_id']:
        if col not in df.columns:
            df[col] = None
    
    # Convert time column to datetime
    if 'time' in df.columns:
        df['time'] = pd.to_datetime(df['time'], errors='coerce')
    
    # Add subject_id
    df['subject_id'] = data['subject_id']
    
    # Reorder columns: guaranteed columns first
    cols = ['subject_id', 'time', 'hadm_id', 'stay_id', 'data_type']
    other_cols = [c for c in df.columns if c not in cols]
    df = df[cols + other_cols]
    
    return df

# ============================================================================
# ENTRY POINT
# ============================================================================

if __name__ == '__main__':
    main()