import os
import pandas as pd
import json
import gzip
import subprocess
from pathlib import Path
from datetime import datetime, timedelta
from tqdm import tqdm
import time
from multiprocessing import Pool, Manager, cpu_count, Lock
from functools import partial

# ============================================================================
# CONFIGURATION - EASY TO MODIFY
# ============================================================================

print("=" * 80)
print("MIMIC-IV PATIENT TRAJECTORY BUILDER - OPTIMIZED VERSION")
print("=" * 80)
print(f"Start time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print()

# *** MODIFY THESE PATHS AS NEEDED ***
TRAJECTORY_VERSION = 0
DATA_DIR = '/ewsc/jpickard/physionet/mimiciv/physionet.org/files/mimiciv/3.1'
OUTPUT_DIR = Path(os.path.join(DATA_DIR, "patient_trajectories_v" + str(TRAJECTORY_VERSION)))
PROGRESS_DIR = Path(os.path.join(OUTPUT_DIR, 'progress'))

# *** MODIFY THESE PERFORMANCE SETTINGS ***
CHUNK_SIZE = 500_000  # Increased from 100k for better performance
BUFFER_SIZE = 10_000  # Records to buffer before writing to disk
MAX_CPUS = min(20, cpu_count() - 1)
NUM_WORKERS = max(1, MAX_CPUS)  # Parallel workers (leave 1 CPU free)

# Create directories
OUTPUT_DIR.mkdir(exist_ok=True)
PROGRESS_DIR.mkdir(exist_ok=True)

print(f"Data directory: {DATA_DIR}/")
print(f"Output directory: {OUTPUT_DIR}/")
print(f"Progress directory: {PROGRESS_DIR}/")
print(f"Chunk size: {CHUNK_SIZE:,} rows")
print(f"Buffer size: {BUFFER_SIZE:,} records")
print(f"Parallel workers: {NUM_WORKERS}")
print()

# ============================================================================
# FILE DEFINITIONS - What to process and how
# ============================================================================

FILE_CONFIG = {
    # Static data (small files, load entirely)
    'static': [
        {
            'path': 'hosp/patients.csv.gz',
            'key': 'demographics',
            'id_col': 'subject_id',
            'time_col': None,
            'merge_strategy': 'single'
        },
    ],
    
    # Admission-level data
    'admissions': [
        {
            'path': 'hosp/admissions.csv.gz',
            'key': 'admissions',
            'id_col': 'subject_id',
            'time_col': 'admittime',
            'merge_strategy': 'list'
        },
        {
            'path': 'icu/icustays.csv.gz',
            'key': 'icu_stays',
            'id_col': 'subject_id',
            'time_col': 'intime',
            'merge_strategy': 'list'
        },
    ],
    
    # Event data (large files, process in chunks)
    'events': [
        {
            'path': 'hosp/labevents.csv.gz',
            'key': 'events',
            'id_col': 'subject_id',
            'time_col': 'charttime',
            'event_type': 'lab',
            'merge_strategy': 'append'
        },
        {
            'path': 'icu/chartevents.csv.gz',
            'key': 'events',
            'id_col': 'subject_id',
            'time_col': 'charttime',
            'event_type': 'chart',
            'merge_strategy': 'append'
        },
        {
            'path': 'icu/inputevents.csv.gz',
            'key': 'events',
            'id_col': 'subject_id',
            'time_col': 'starttime',
            'event_type': 'input',
            'merge_strategy': 'append'
        },
        {
            'path': 'icu/outputevents.csv.gz',
            'key': 'events',
            'id_col': 'subject_id',
            'time_col': 'charttime',
            'event_type': 'output',
            'merge_strategy': 'append'
        },
        {
            'path': 'hosp/prescriptions.csv.gz',
            'key': 'events',
            'id_col': 'subject_id',
            'time_col': 'starttime',
            'event_type': 'prescription',
            'merge_strategy': 'append'
        },
        {
            'path': 'hosp/diagnoses_icd.csv.gz',
            'key': 'diagnoses',
            'id_col': 'subject_id',
            'time_col': None,
            'event_type': 'diagnosis',
            'merge_strategy': 'append'
        },
        {
            'path': 'hosp/procedures_icd.csv.gz',
            'key': 'procedures',
            'id_col': 'subject_id',
            'time_col': None,
            'event_type': 'procedure',
            'merge_strategy': 'append'
        },
        {
            'path': 'hosp/microbiologyevents.csv.gz',
            'key': 'events',
            'id_col': 'subject_id',
            'time_col': 'charttime',
            'event_type': 'microbiology',
            'merge_strategy': 'append'
        },
    ]
}

# ============================================================================
# UTILITY FUNCTIONS
# ============================================================================

def count_file_lines(file_path):
    """Quickly count lines in a gzipped file."""
    try:
        result = subprocess.run(
            f"gunzip -c '{file_path}' | wc -l",
            shell=True,
            capture_output=True,
            text=True,
            timeout=300
        )
        if result.returncode == 0:
            return int(result.stdout.strip())
    except:
        pass
    
    # Fallback: manual counting
    count = 0
    with gzip.open(file_path, 'rt') as f:
        for _ in f:
            count += 1
    return count

def get_file_line_count_cached(file_path):
    """Get line count with caching."""
    cache_file = os.path.join(PROGRESS_DIR, f"{Path(file_path).name}.linecount")
    
    if Path(cache_file).exists():
        with open(cache_file, 'r') as f:
            return int(f.read().strip())
    
    print(f"    Counting lines in {Path(file_path).name}... ", end='', flush=True)
    count = count_file_lines(file_path)
    print(f"{count:,} lines")
    
    with open(cache_file, 'w') as f:
        f.write(str(count))
    
    return count

# ============================================================================
# OPTIMIZED PATIENT BUFFER - Batches writes for 10x speed improvement
# ============================================================================

class PatientBuffer:
    """
    Buffers writes to patient files to dramatically reduce disk I/O.
    Instead of opening/closing files for every record, we batch writes.
    """
    
    def __init__(self, output_dir, buffer_size=10000):
        self.output_dir = output_dir
        self.buffer_size = buffer_size
        self.buffers = {}  # {patient_id: [records]}
        self.total_buffered = 0
        
    def add(self, patient_id, data_type, data):
        """Add a record to the buffer"""
        if patient_id not in self.buffers:
            self.buffers[patient_id] = []
        
        self.buffers[patient_id].append({
            'data_type': data_type,
            'data': data
        })
        self.total_buffered += 1
        
        # Auto-flush when buffer gets large
        if self.total_buffered >= self.buffer_size:
            self.flush_all()
    
    def flush(self, patient_id):
        """Flush a specific patient's buffer to disk"""
        if patient_id in self.buffers and self.buffers[patient_id]:
            patient_file = os.path.join(self.output_dir, f'patient_{patient_id}.jsonl')
            with open(patient_file, 'a') as f:
                for record in self.buffers[patient_id]:
                    f.write(json.dumps(record, default=str) + '\n')
            
            count = len(self.buffers[patient_id])
            self.buffers[patient_id] = []
            self.total_buffered -= count
    
    def flush_all(self):
        """Flush all buffers to disk"""
        for patient_id in list(self.buffers.keys()):
            self.flush(patient_id)
        self.total_buffered = 0
    
    def __len__(self):
        """Total records currently buffered"""
        return self.total_buffered

# ============================================================================
# PROGRESS TRACKING - Thread-safe for parallel processing
# ============================================================================

class ProgressTracker:
    """
    Tracks which files and chunks have been processed.
    Thread-safe for parallel processing.
    """
    
    def __init__(self, progress_file):
        self.progress_file = progress_file
        self.file_progress = {}
        self.load()
    
    def load(self):
        """Load progress from file"""
        if Path(self.progress_file).exists():
            try:
                with open(self.progress_file, 'r') as f:
                    data = json.load(f)
                    self.file_progress = data.get('file_progress', {})
            except:
                self.file_progress = {}
    
    def save(self):
        """Save progress to file"""
        try:
            # Use file locking for process safety
            import fcntl
            with open(self.progress_file, 'w') as f:
                fcntl.flock(f.fileno(), fcntl.LOCK_EX)
                try:
                    json.dump({
                        'file_progress': self.file_progress,
                        'last_update': datetime.now().isoformat()
                    }, f, indent=2)
                finally:
                    fcntl.flock(f.fileno(), fcntl.LOCK_UN)
        except ImportError:
            # Fallback for Windows (no fcntl)
            with open(self.progress_file, 'w') as f:
                json.dump({
                    'file_progress': self.file_progress,
                    'last_update': datetime.now().isoformat()
                }, f, indent=2)
        except Exception as e:
            print(f"Warning: Could not save progress: {e}")
    
    def get_file_progress(self, file_path):
        """Get progress for a specific file"""
        return self.file_progress.get(str(file_path), {
            'status': 'not_started',
            'last_chunk': 0,
            'total_chunks': None,
            'records_written': 0
        })
    
    def update_file_progress(self, file_path, chunk_num, total_chunks, records_written):
        """Update progress for a file"""
        file_key = str(file_path)
        if file_key not in self.file_progress:
            self.file_progress[file_key] = {}
        
        self.file_progress[file_key].update({
            'status': 'in_progress',
            'last_chunk': chunk_num,
            'total_chunks': total_chunks,
            'records_written': records_written,
            'last_update': datetime.now().isoformat()
        })
        self.save()
    
    def mark_file_complete(self, file_path, records_written):
        """Mark a file as completed"""
        file_key = str(file_path)
        if file_key not in self.file_progress:
            self.file_progress[file_key] = {}
        
        self.file_progress[file_key].update({
            'status': 'complete',
            'records_written': records_written,
            'completed_at': datetime.now().isoformat()
        })
        self.save()
    
    def is_file_complete(self, file_path):
        """Check if a file has been fully processed"""
        prog = self.get_file_progress(file_path)
        return prog.get('status') == 'complete'
    
    def get_overall_progress(self):
        """Get overall progress statistics"""
        total_files = len(self.file_progress)
        completed = sum(1 for p in self.file_progress.values() if p.get('status') == 'complete')
        in_progress = sum(1 for p in self.file_progress.values() if p.get('status') == 'in_progress')
        total_records = sum(p.get('records_written', 0) for p in self.file_progress.values())
        
        return {
            'total_files': total_files,
            'completed': completed,
            'in_progress': in_progress,
            'total_records': total_records
        }

# Global progress tracker (will be created per process)
# progress = None

def init_worker(progress_file):
    """Initialize worker process with its own progress tracker"""
    global progress
    progress = ProgressTracker(progress_file)

# ============================================================================
# FILE PROCESSING FUNCTIONS - Optimized versions
# ============================================================================

def process_static_file(file_config, icu_patient_ids, output_dir):
    """Process static files (demographics) - small files loaded entirely"""
    file_path = os.path.join(DATA_DIR, file_config['path'])
    
    print(f"  [{Path(file_path).name}] Loading...")
    
    df = pd.read_csv(file_path)
    df = df[df[file_config['id_col']].isin(icu_patient_ids)]
    
    buffer = PatientBuffer(output_dir, buffer_size=BUFFER_SIZE)
    
    # Use to_dict('records') instead of iterrows for speed
    for record in df.to_dict('records'):
        patient_id = record[file_config['id_col']]
        buffer.add(patient_id, file_config['key'], record)
    
    buffer.flush_all()
    count = len(df)
    
    print(f"  [{Path(file_path).name}] ✓ Wrote {count:,} records")
    return count

def process_chunked_file(file_config, icu_patient_ids, output_dir, progress_tracker):
    """
    Process large files in chunks with optimizations:
    - Batched writes via buffer
    - Vectorized pandas operations
    - Resume capability
    """
    file_path = os.path.join(DATA_DIR, file_config['path'])
    
    print(f"  [{Path(file_path).name}] Starting...")
    
    # Get previous progress
    file_prog = progress_tracker.get_file_progress(file_path)
    start_chunk = file_prog.get('last_chunk', 0)
    total_records = file_prog.get('records_written', 0)
    
    # Get total number of chunks
    if file_prog.get('total_chunks') is None:
        line_count = get_file_line_count_cached(file_path)
        total_chunks = max(1, (line_count - 1) // CHUNK_SIZE)
    else:
        total_chunks = file_prog['total_chunks']
    
    if start_chunk > 0:
        print(f"  [{Path(file_path).name}] Resuming from chunk {start_chunk + 1}/{total_chunks}")
        print(f"  [{Path(file_path).name}] Already processed: {total_records:,} records")
    
    buffer = PatientBuffer(output_dir, buffer_size=BUFFER_SIZE)
    
    try:
        # Create chunk iterator
        chunk_iter = pd.read_csv(file_path, chunksize=CHUNK_SIZE)
        
        # Skip chunks we've already processed
        for _ in range(start_chunk):
            next(chunk_iter)
        
        # Process remaining chunks with progress bar
        pbar = tqdm(
            total=total_chunks,
            initial=start_chunk,
            desc=f"  [{Path(file_path).name}]",
            position=0,
            leave=True
        )
        
        for chunk_idx, chunk in enumerate(chunk_iter, start=start_chunk):
            # Filter to ICU patients only
            chunk = chunk[chunk[file_config['id_col']].isin(icu_patient_ids)]
            
            if len(chunk) == 0:
                progress_tracker.update_file_progress(
                    file_path, chunk_idx + 1, total_chunks, total_records
                )
                pbar.update(1)
                continue
            
            # Add event_type if specified
            if 'event_type' in file_config:
                chunk['event_type'] = file_config['event_type']
            
            # Vectorized processing - much faster than iterrows
            for subject_id, group in chunk.groupby(file_config['id_col']):
                records = group.to_dict('records')  # Fast vectorized operation
                for data in records:
                    buffer.add(subject_id, file_config['key'], data)
                    total_records += 1
            
            # Flush buffer and save progress after each chunk
            buffer.flush_all()
            progress_tracker.update_file_progress(
                file_path, chunk_idx + 1, total_chunks, total_records
            )
            
            pbar.update(1)
            pbar.set_postfix({'records': f'{total_records:,}'})
        
        pbar.close()
    
    except Exception as e:
        buffer.flush_all()
        print(f"  [{Path(file_path).name}] ✗ Error: {e}")
        return total_records
    
    print(f"  [{Path(file_path).name}] ✓ Complete: {total_records:,} records")
    return total_records

# ============================================================================
# PARALLEL PROCESSING WRAPPER
# ============================================================================

def process_file_wrapper(args):
    """
    Wrapper function for parallel processing.
    Processes a single file and returns results.
    """
    global progress
    file_config, icu_patient_ids, output_dir, progress_file = args
#    progress_file = Path(progress_file)
    # Create progress tracker if not exists
    if progress is None:
        progress = ProgressTracker(progress_file)
    
    file_path = Path(os.path.join(DATA_DIR, file_config['path']))
    
    # Check if already complete
    if progress.is_file_complete(file_path):
        prog = progress.get_file_progress(file_path)
        records = prog.get('records_written', 0)
        print(f"  [{Path(file_path).name}] ⊘ Already complete ({records:,} records)")
        return {
            'file': str(file_path),
            'status': 'skipped',
            'records': records
        }
    
    # Check if file exists
    if not Path(file_path).exists():
        print(f"  [{Path(file_path).name}] ⚠ File not found")
        return {
            'file': str(file_path),
            'status': 'not_found',
            'records': 0
        }
    
    try:
        # Determine if this is a chunked file or static file
        if file_config in FILE_CONFIG['events']:
            count = process_chunked_file(file_config, icu_patient_ids, output_dir, progress)
        else:
            count = process_static_file(file_config, icu_patient_ids, output_dir)
        
        progress.mark_file_complete(file_path, count)
        
        return {
            'file': str(file_path),
            'status': 'complete',
            'records': count
        }
    
    except Exception as e:
        print(f"  [{Path(file_path).name}] ✗ Error: {e}")
        return {
            'file': str(file_path),
            'status': 'error',
            'error': str(e),
            'records': 0
        }

print("Setup Complete")

global progress
overall_start = time.time()

# Initialize progress tracker in main process
progress = ProgressTracker(os.path.join(PROGRESS_DIR, 'build_progress.json'))

# ========================================================================
# STEP 1: GET ALL ICU PATIENTS
# ========================================================================

print("STEP 1: Identifying all ICU patients...")

icustays_file = os.path.join(DATA_DIR, 'icu', 'icustays.csv.gz')
icustays = pd.read_csv(icustays_file)
icu_patient_ids = set(icustays['subject_id'].unique())

print(f"✓ Found {len(icu_patient_ids):,} unique ICU patients")
print()

# ========================================================================
# STEP 2: INITIALIZE PATIENT FILES
# ========================================================================

print("STEP 2: Initializing patient trajectory files...")

existing_patients = set()
for f in OUTPUT_DIR.glob('patient_*.jsonl'):
    try:
        patient_id = int(f.stem.split('_')[1])
        existing_patients.add(patient_id)
    except:
        pass

for f in OUTPUT_DIR.glob('patient_*.json'):
    try:
        patient_id = int(f.stem.split('_')[1])
        existing_patients.add(patient_id)
    except:
        pass

patients_to_initialize = icu_patient_ids - existing_patients

if len(patients_to_initialize) > 0:
    print(f"Creating {len(patients_to_initialize):,} new patient files...")
    for patient_id in tqdm(patients_to_initialize, desc="Initializing"):
        patient_file = Path(os.path.join(OUTPUT_DIR, f'patient_{patient_id}.jsonl'))
        patient_file.touch()
    print(f"✓ Initialized {len(patients_to_initialize):,} patient files")
else:
    print("✓ All patient files already exist")

print(f"✓ Total patient files: {len(icu_patient_ids):,}")
print()

# ========================================================================
# STEP 3: PROCESS CSV FILES (WITH PARALLELIZATION)
# ========================================================================

print("STEP 3: Processing all CSV files...")
print(f"Using {NUM_WORKERS} parallel workers")
print()

# Prepare all file configs
all_file_configs = (
    FILE_CONFIG['static'] + 
    FILE_CONFIG['admissions'] + 
    FILE_CONFIG['events']
)

# For large event files, process sequentially with internal chunking
# For smaller files, can process in parallel

print("Processing static and admission files in parallel...")
small_files = FILE_CONFIG['static'] + FILE_CONFIG['admissions']

if len(small_files) > 0 and NUM_WORKERS > 1:
    args_list = [
        (config, icu_patient_ids, Path(OUTPUT_DIR), os.path.join(PROGRESS_DIR, 'build_progress.json'))
        for config in small_files
    ]
    
    with Pool(NUM_WORKERS) as pool:
        results = pool.map(process_file_wrapper, args_list)
else:
    # Sequential fallback
    results = []
    for config in small_files:
        result = process_file_wrapper((config, icu_patient_ids, OUTPUT_DIR, os.path.join(PROGRESS_DIR, 'build_progress.json')))
        results.append(result)

print()
print("Processing large event files (sequential with chunking)...")

# Process large event files sequentially (they already use chunking efficiently)
for file_config in FILE_CONFIG['events']:
    result = process_file_wrapper((file_config, icu_patient_ids, OUTPUT_DIR, os.path.join(PROGRESS_DIR, 'build_progress.json')))

print()

# ========================================================================
# PROGRESS SUMMARY
# ========================================================================

overall_progress = progress.get_overall_progress()
print("=" * 80)
print("FILE PROCESSING SUMMARY")
print("=" * 80)
print(f"Total files tracked: {overall_progress['total_files']}")
print(f"Completed: {overall_progress['completed']}")
print(f"In progress: {overall_progress['in_progress']}")
print(f"Total records written: {overall_progress['total_records']:,}")
print()