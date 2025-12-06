"""
MIMIC-IV Steroid Analysis
==========================
Comprehensive analysis of steroid administration
in sepsis vs non-sepsis ICU patients.
"""

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

from clinical_constants import STEROIDS
from utils import identify_sepsis_bool

# ============================================================================
# CONFIGURATION
# ============================================================================
from config import *

OUTPUT_DIR = EDA_DIR / f"steroids"
FIGURES_DIR = OUTPUT_DIR / 'figures'
CACHE_DIR = OUTPUT_DIR / 'cache'

# Create directories if they don't exist
for d in [OUTPUT_DIR, FIGURES_DIR, CACHE_DIR]:
    d.mkdir(parents=True, exist_ok=True)

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler(OUTPUT_DIR / f'steroid_analysis_{datetime.now().strftime("%Y%m%d_%H%M%S")}.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# Plotting style
sns.set_style("whitegrid")
plt.rcParams['figure.figsize'] = (14, 8)
plt.rcParams['font.size'] = 10


def classify_steroid(drug_name: str) -> tuple:
    """
    Classify a drug as a steroid.
    Returns: (is_steroid: bool, steroid_type: str or None)
    """
    if not drug_name or not isinstance(drug_name, str):
        return False, None
    
    drug_lower = drug_name.lower().strip()
    
    for steroid_type, keywords in STEROIDS.items():
        for keyword in keywords:
            if keyword in drug_lower:
                return True, steroid_type
    
    return False, None


def extract_patient_steroids(patient_file: Path) -> dict:
    """
    Extract steroid medication data from a patient file.
    Focus on steroids during ICU stays.
    """
    try:
        with open(patient_file, 'r') as f:
            data = json.load(f)
        
        subject_id = data.get('subject_id')
        has_sepsis = identify_sepsis_bool(data)
        
        records = data.get('records', [])
        if not isinstance(records, list):
            records = []
        
        # Get ICU stay IDs
        icu_stays = set()
        for record in records:
            if isinstance(record, dict) and record.get('data_type') == 'icu_stays':
                stay_id = record.get('stay_id')
                if stay_id:
                    icu_stays.add(stay_id)
        
        # Extract steroid events
        steroids = []
        
        for record in records:
            if not isinstance(record, dict):
                continue
            
            # Look at prescription and input events
            event_type = record.get('event_type')
            if event_type not in ['prescription', 'input']:
                continue
            
            stay_id = record.get('stay_id')
            drug_name = record.get('drug') or record.get('ordercategorydescription', '')
            
            # Classify medication
            is_steroid, steroid_type = classify_steroid(drug_name)
            
            if is_steroid:
                steroids.append({
                    'subject_id': subject_id,
                    'stay_id': stay_id,
                    'drug': drug_name,
                    'steroid_type': steroid_type,
                    'starttime': record.get('starttime'),
                    'stoptime': record.get('stoptime'),
                    'dose_val_rx': record.get('dose_val_rx'),
                    'dose_unit_rx': record.get('dose_unit_rx'),
                    'rate': record.get('rate'),
                    'rate_uom': record.get('rate_uom'),
                    'amount': record.get('amount'),
                    'amount_uom': record.get('amountuom'),
                    'event_type': event_type,
                    'route': record.get('route'),
                })
        
        return {
            'subject_id': subject_id,
            'has_sepsis': has_sepsis,
            'num_icu_stays': len(icu_stays),
            'steroids': steroids,
        }
        
    except Exception as e:
        logger.error(f"Error extracting steroids from {patient_file.name}: {e}")
        return {
            'subject_id': None,
            'error': str(e)
        }


def extract_steroids_wrapper(patient_file: Path) -> dict:
    """Wrapper for parallel processing."""
    return extract_patient_steroids(patient_file)


# ============================================================================
# ANALYSIS FUNCTIONS
# ============================================================================

class SteroidAnalyzer:
    """Comprehensive steroid analysis."""
    
    def __init__(self, steroid_data: list):
        """Initialize with list of patient steroid dictionaries."""
        self.data = steroid_data
        
        # Create DataFrame
        all_steroids = []
        
        for patient in steroid_data:
            if patient.get('subject_id') is None:
                continue
            
            has_sepsis = patient.get('has_sepsis', False)
            subject_id = patient['subject_id']
            
            for steroid in patient.get('steroids', []):
                steroid['has_sepsis'] = has_sepsis
                all_steroids.append(steroid)
        
        self.steroid_df = pd.DataFrame(all_steroids)
        
        # Patient-level summary
        self.patient_summary = self._create_patient_summary()
        
        logger.info(f"Loaded {len(self.steroid_df):,} steroid records")
        logger.info(f"Patient summary: {len(self.patient_summary):,} patients")
    
    def _create_patient_summary(self) -> pd.DataFrame:
        """Create patient-level steroid summary."""
        patient_data = []
        
        for patient in self.data:
            if patient.get('subject_id') is None:
                continue
            
            subject_id = patient['subject_id']
            has_sepsis = patient.get('has_sepsis', False)
            
            # Count steroids
            steroids = patient.get('steroids', [])
            steroid_types = Counter(s['steroid_type'] for s in steroids if s.get('steroid_type'))
            
            # Get routes of administration
            routes = Counter(s['route'] for s in steroids if s.get('route'))
            
            patient_data.append({
                'subject_id': subject_id,
                'has_sepsis': has_sepsis,
                'num_icu_stays': patient.get('num_icu_stays', 0),
                'received_steroids': len(steroids) > 0,
                'num_steroid_administrations': len(steroids),
                'num_unique_steroids': len(steroid_types),
                'steroid_types': list(steroid_types.keys()),
                'routes': list(routes.keys()),
            })
        
        return pd.DataFrame(patient_data)
    
    def generate_summary_statistics(self) -> dict:
        """Generate comprehensive summary statistics."""
        
        sepsis_patients = self.patient_summary[self.patient_summary['has_sepsis'] == True]
        non_sepsis_patients = self.patient_summary[self.patient_summary['has_sepsis'] == False]
        
        summary = {
            'cohort_overview': {
                'total_patients': len(self.patient_summary),
                'sepsis_patients': len(sepsis_patients),
                'non_sepsis_patients': len(non_sepsis_patients),
            },
            
            'steroid_usage': {
                'overall': {
                    'patients_received': int(self.patient_summary['received_steroids'].sum()),
                    'percentage': float(self.patient_summary['received_steroids'].mean() * 100),
                    'total_administrations': len(self.steroid_df),
                    'avg_per_patient': float(self.patient_summary['num_steroid_administrations'].mean()),
                },
                'sepsis': {
                    'patients_received': int(sepsis_patients['received_steroids'].sum()) if len(sepsis_patients) > 0 else 0,
                    'percentage': float(sepsis_patients['received_steroids'].mean() * 100) if len(sepsis_patients) > 0 else 0,
                    'avg_administrations': float(sepsis_patients['num_steroid_administrations'].mean()) if len(sepsis_patients) > 0 else 0,
                },
                'non_sepsis': {
                    'patients_received': int(non_sepsis_patients['received_steroids'].sum()) if len(non_sepsis_patients) > 0 else 0,
                    'percentage': float(non_sepsis_patients['received_steroids'].mean() * 100) if len(non_sepsis_patients) > 0 else 0,
                    'avg_administrations': float(non_sepsis_patients['num_steroid_administrations'].mean()) if len(non_sepsis_patients) > 0 else 0,
                },
            },
            
            'steroid_types': self._analyze_steroid_types(),
            'administration_routes': self._analyze_routes(),
            'timing_analysis': self._analyze_timing(),
        }
        
        return summary
    
    def _analyze_steroid_types(self) -> dict:
        """Analyze distribution of steroid types."""
        result = {}
        
        if len(self.steroid_df) == 0:
            return result
        
        # Overall distribution
        overall_counts = self.steroid_df['steroid_type'].value_counts().to_dict()
        result['overall_distribution'] = overall_counts
        
        # By sepsis status
        sepsis_steroid = self.steroid_df[self.steroid_df['has_sepsis'] == True]
        non_sepsis_steroid = self.steroid_df[self.steroid_df['has_sepsis'] == False]
        
        result['sepsis_distribution'] = sepsis_steroid['steroid_type'].value_counts().to_dict()
        result['non_sepsis_distribution'] = non_sepsis_steroid['steroid_type'].value_counts().to_dict()
        
        return result
    
    def _analyze_routes(self) -> dict:
        """Analyze routes of administration."""
        result = {}
        
        if len(self.steroid_df) == 0:
            return result
        
        routes_df = self.steroid_df[self.steroid_df['route'].notna()]
        
        if len(routes_df) == 0:
            return result
        
        result['overall_distribution'] = routes_df['route'].value_counts().to_dict()
        
        sepsis_routes = routes_df[routes_df['has_sepsis'] == True]
        non_sepsis_routes = routes_df[routes_df['has_sepsis'] == False]
        
        result['sepsis_distribution'] = sepsis_routes['route'].value_counts().to_dict()
        result['non_sepsis_distribution'] = non_sepsis_routes['route'].value_counts().to_dict()
        
        return result
    
    def _analyze_timing(self) -> dict:
        """Analyze timing of steroid administration."""
        result = {}
        
        if len(self.steroid_df) == 0:
            return result
        
        # Calculate duration for records with both start and stop times
        timed_df = self.steroid_df[
            self.steroid_df['starttime'].notna() & 
            self.steroid_df['stoptime'].notna()
        ].copy()
        
        if len(timed_df) > 0:
            timed_df['starttime'] = pd.to_datetime(timed_df['starttime'], errors='coerce')
            timed_df['stoptime'] = pd.to_datetime(timed_df['stoptime'], errors='coerce')
            timed_df['duration_hours'] = (timed_df['stoptime'] - timed_df['starttime']).dt.total_seconds() / 3600
            
            # Remove negative or extreme durations
            timed_df = timed_df[(timed_df['duration_hours'] > 0) & (timed_df['duration_hours'] < 720)]  # < 30 days
            
            if len(timed_df) > 0:
                result['overall_duration'] = {
                    'mean_hours': float(timed_df['duration_hours'].mean()),
                    'median_hours': float(timed_df['duration_hours'].median()),
                    'std_hours': float(timed_df['duration_hours'].std()),
                }
                
                sepsis_timed = timed_df[timed_df['has_sepsis'] == True]
                non_sepsis_timed = timed_df[timed_df['has_sepsis'] == False]
                
                if len(sepsis_timed) > 0:
                    result['sepsis_duration'] = {
                        'mean_hours': float(sepsis_timed['duration_hours'].mean()),
                        'median_hours': float(sepsis_timed['duration_hours'].median()),
                        'std_hours': float(sepsis_timed['duration_hours'].std()),
                    }
                
                if len(non_sepsis_timed) > 0:
                    result['non_sepsis_duration'] = {
                        'mean_hours': float(non_sepsis_timed['duration_hours'].mean()),
                        'median_hours': float(non_sepsis_timed['duration_hours'].median()),
                        'std_hours': float(non_sepsis_timed['duration_hours'].std()),
                    }
        
        return result
    
    def plot_steroid_distribution(self, output_dir: Path):
        """Plot steroid usage distribution."""
        if len(self.steroid_df) == 0:
            logger.warning("No steroid data to plot")
            return
        
        fig, axes = plt.subplots(2, 2, figsize=(16, 12))
        
        # 1. Overall steroid type distribution
        steroid_counts = self.steroid_df['steroid_type'].value_counts()
        axes[0, 0].barh(steroid_counts.index, steroid_counts.values, color='steelblue')
        axes[0, 0].set_xlabel('Number of Administrations')
        axes[0, 0].set_title('Steroid Type Distribution (All Patients)')
        axes[0, 0].grid(alpha=0.3, axis='x')
        
        # 2. Steroid usage: Sepsis vs Non-Sepsis
        sepsis_steroid = self.steroid_df[self.steroid_df['has_sepsis'] == True]
        non_sepsis_steroid = self.steroid_df[self.steroid_df['has_sepsis'] == False]
        
        sepsis_counts = sepsis_steroid['steroid_type'].value_counts()
        non_sepsis_counts = non_sepsis_steroid['steroid_type'].value_counts()
        
        all_types = sorted(set(list(sepsis_counts.index) + list(non_sepsis_counts.index)))
        x = np.arange(len(all_types))
        width = 0.35
        
        sepsis_vals = [sepsis_counts.get(t, 0) for t in all_types]
        non_sepsis_vals = [non_sepsis_counts.get(t, 0) for t in all_types]
        
        axes[0, 1].bar(x - width/2, sepsis_vals, width, label='Sepsis', color='coral')
        axes[0, 1].bar(x + width/2, non_sepsis_vals, width, label='Non-Sepsis', color='lightblue')
        axes[0, 1].set_xlabel('Steroid Type')
        axes[0, 1].set_ylabel('Number of Administrations')
        axes[0, 1].set_title('Steroid Usage by Sepsis Status')
        axes[0, 1].set_xticks(x)
        axes[0, 1].set_xticklabels(all_types, rotation=45, ha='right')
        axes[0, 1].legend()
        axes[0, 1].grid(alpha=0.3, axis='y')
        
        # 3. Patient-level: Percentage receiving steroids
        sepsis_pts = self.patient_summary[self.patient_summary['has_sepsis'] == True]
        non_sepsis_pts = self.patient_summary[self.patient_summary['has_sepsis'] == False]
        
        pct_data = pd.DataFrame({
            'Received': [
                sepsis_pts['received_steroids'].sum(),
                non_sepsis_pts['received_steroids'].sum()
            ],
            'Not Received': [
                len(sepsis_pts) - sepsis_pts['received_steroids'].sum(),
                len(non_sepsis_pts) - non_sepsis_pts['received_steroids'].sum()
            ]
        }, index=['Sepsis', 'Non-Sepsis'])
        
        pct_data.plot(kind='bar', stacked=True, ax=axes[1, 0], color=['coral', 'lightgray'])
        axes[1, 0].set_ylabel('Number of Patients')
        axes[1, 0].set_title('Steroid Receipt by Cohort')
        axes[1, 0].set_xticklabels(['Sepsis', 'Non-Sepsis'], rotation=0)
        axes[1, 0].legend(title='Steroid')
        axes[1, 0].grid(alpha=0.3, axis='y')
        
        # 4. Number of steroid administrations per patient
        sepsis_num = sepsis_pts[sepsis_pts['received_steroids']]['num_steroid_administrations']
        non_sepsis_num = non_sepsis_pts[non_sepsis_pts['received_steroids']]['num_steroid_administrations']
        
        if len(sepsis_num) > 0 or len(non_sepsis_num) > 0:
            axes[1, 1].hist([non_sepsis_num, sepsis_num], bins=30, label=['Non-Sepsis', 'Sepsis'], 
                           color=['lightblue', 'coral'], alpha=0.7, edgecolor='black')
            axes[1, 1].set_xlabel('Number of Steroid Administrations')
            axes[1, 1].set_ylabel('Number of Patients')
            axes[1, 1].set_title('Steroid Administrations per Patient')
            axes[1, 1].legend()
            axes[1, 1].grid(alpha=0.3, axis='y')
            axes[1, 1].set_xlim(0, min(50, max(sepsis_num.max() if len(sepsis_num) > 0 else 0, 
                                              non_sepsis_num.max() if len(non_sepsis_num) > 0 else 0)))
        
        plt.tight_layout()
        plt.savefig(output_dir / 'steroid_distribution.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        logger.info("Saved steroid distribution plot")
    
    def plot_route_analysis(self, output_dir: Path):
        """Plot administration route analysis."""
        if len(self.steroid_df) == 0:
            logger.warning("No steroid data to plot routes")
            return
        
        routes_df = self.steroid_df[self.steroid_df['route'].notna()]
        
        if len(routes_df) == 0:
            logger.warning("No route data available")
            return
        
        fig, axes = plt.subplots(1, 2, figsize=(16, 6))
        
        # 1. Overall route distribution
        route_counts = routes_df['route'].value_counts().head(10)
        axes[0].barh(route_counts.index, route_counts.values, color='steelblue')
        axes[0].set_xlabel('Number of Administrations')
        axes[0].set_title('Top 10 Administration Routes (All Patients)')
        axes[0].grid(alpha=0.3, axis='x')
        
        # 2. Routes by sepsis status
        sepsis_routes = routes_df[routes_df['has_sepsis'] == True]
        non_sepsis_routes = routes_df[routes_df['has_sepsis'] == False]
        
        sepsis_route_counts = sepsis_routes['route'].value_counts().head(10)
        non_sepsis_route_counts = non_sepsis_routes['route'].value_counts().head(10)
        
        all_routes = sorted(set(list(sepsis_route_counts.index) + list(non_sepsis_route_counts.index)))
        x = np.arange(len(all_routes))
        width = 0.35
        
        sepsis_vals = [sepsis_route_counts.get(r, 0) for r in all_routes]
        non_sepsis_vals = [non_sepsis_route_counts.get(r, 0) for r in all_routes]
        
        axes[1].bar(x - width/2, sepsis_vals, width, label='Sepsis', color='coral')
        axes[1].bar(x + width/2, non_sepsis_vals, width, label='Non-Sepsis', color='lightblue')
        axes[1].set_xlabel('Administration Route')
        axes[1].set_ylabel('Number of Administrations')
        axes[1].set_title('Administration Routes by Sepsis Status')
        axes[1].set_xticks(x)
        axes[1].set_xticklabels(all_routes, rotation=45, ha='right')
        axes[1].legend()
        axes[1].grid(alpha=0.3, axis='y')
        
        plt.tight_layout()
        plt.savefig(output_dir / 'steroid_routes.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        logger.info("Saved steroid routes plot")
    
    def plot_duration_analysis(self, output_dir: Path):
        """Plot steroid duration analysis."""
        timed_df = self.steroid_df[
            self.steroid_df['starttime'].notna() & 
            self.steroid_df['stoptime'].notna()
        ].copy()
        
        if len(timed_df) == 0:
            logger.warning("No timing data available")
            return
        
        timed_df['starttime'] = pd.to_datetime(timed_df['starttime'], errors='coerce')
        timed_df['stoptime'] = pd.to_datetime(timed_df['stoptime'], errors='coerce')
        timed_df['duration_hours'] = (timed_df['stoptime'] - timed_df['starttime']).dt.total_seconds() / 3600
        
        # Remove negative or extreme durations
        timed_df = timed_df[(timed_df['duration_hours'] > 0) & (timed_df['duration_hours'] < 720)]
        
        if len(timed_df) == 0:
            logger.warning("No valid duration data after filtering")
            return
        
        fig, axes = plt.subplots(1, 2, figsize=(16, 6))
        
        # 1. Duration distribution histogram
        sepsis_dur = timed_df[timed_df['has_sepsis'] == True]['duration_hours']
        non_sepsis_dur = timed_df[timed_df['has_sepsis'] == False]['duration_hours']
        
        axes[0].hist([non_sepsis_dur, sepsis_dur], bins=50, label=['Non-Sepsis', 'Sepsis'],
                     color=['lightblue', 'coral'], alpha=0.7, edgecolor='black')
        axes[0].set_xlabel('Duration (hours)')
        axes[0].set_ylabel('Number of Administrations')
        axes[0].set_title('Steroid Administration Duration')
        axes[0].legend()
        axes[0].grid(alpha=0.3, axis='y')
        
        # 2. Duration boxplot by sepsis status
        combined = timed_df[['has_sepsis', 'duration_hours']].copy()
        combined['has_sepsis'] = combined['has_sepsis'].map({True: 'Sepsis', False: 'Non-Sepsis'})
        
        sns.boxplot(x='has_sepsis', y='duration_hours', data=combined, ax=axes[1])
        axes[1].set_xlabel('Cohort')
        axes[1].set_ylabel('Duration (hours)')
        axes[1].set_title('Steroid Duration by Cohort')
        axes[1].grid(alpha=0.3, axis='y')
        
        plt.tight_layout()
        plt.savefig(output_dir / 'steroid_duration.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        logger.info("Saved steroid duration plot")


# ============================================================================
# MAIN SCRIPT
# ============================================================================

def load_patient_files(cleaned_dir: Path, sample_size: int = None) -> list:
    """
    Load all patient trajectory JSON files.
    Optionally subsample for speed.
    """
    patient_files = sorted(cleaned_dir.glob("*.json"))
    
    if sample_size is not None and sample_size < len(patient_files):
        np.random.seed(123)
        patient_files = list(np.random.choice(patient_files, size=sample_size, replace=False))
    
    logger.info(f"Loading {len(patient_files):,} patient files")
    return patient_files


def to_serializable(obj):
    """Convert numpy types to native Python types for JSON serialization."""
    if isinstance(obj, (np.integer, np.int64)):
        return int(obj)
    if isinstance(obj, (np.floating, np.float64)):
        return float(obj)
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    return obj


def run_parallel_extraction(patient_files: list, num_workers: int = 4) -> list:
    """
    Run steroid extraction in parallel using multiprocessing.
    """
    logger.info(f"Starting parallel extraction with {num_workers} workers...")
    
    start = time.time()
    with Pool(num_workers) as pool:
        results = list(tqdm(pool.imap(extract_steroids_wrapper, patient_files),
                            total=len(patient_files)))
    end = time.time()
    
    logger.info(f"Completed extraction in {end - start:.2f} seconds")
    return results


def save_raw_extraction_results(results: list, cache_dir: Path):
    """
    Save raw extraction results to disk for reuse.
    """
    cache_file = cache_dir / "steroid_extraction_results.json"
    
    with open(cache_file, "w") as f:
        json.dump(results, f)
    
    logger.info(f"Saved steroid extraction results to {cache_file}")


def load_cached_results(cache_dir: Path) -> list | None:
    """
    Load previously saved extraction results if available.
    """
    cache_file = cache_dir / "steroid_extraction_results.json"
    
    if cache_file.exists():
        logger.info(f"Loading cached extraction results from {cache_file}")
        with open(cache_file, "r") as f:
            return json.load(f)
    else:
        return None


def main_steroids():
    logger.info("====================================================")
    logger.info("      MIMIC-IV Steroid Utilization Analysis        ")
    logger.info("====================================================")
    
    # ------------------------------------------------------------
    # LOAD PATIENT FILES
    # ------------------------------------------------------------
    patient_files = load_patient_files(CLEANED_DIR, sample_size=SAMPLE_SIZE)
    
    # ------------------------------------------------------------
    # LOAD OR CREATE EXTRACTION CACHE
    # ------------------------------------------------------------
    cached = load_cached_results(CACHE_DIR)
    
    if cached is not None and False:
        steroid_data = cached
    else:
        steroid_data = run_parallel_extraction(patient_files, NUM_WORKERS)
        save_raw_extraction_results(steroid_data, CACHE_DIR)
    
    # ------------------------------------------------------------
    # ANALYZE
    # ------------------------------------------------------------
    logger.info("Initializing SteroidAnalyzer...")
    analyzer = SteroidAnalyzer(steroid_data)
    
    logger.info("Generating summary statistics...")
    summary_stats = analyzer.generate_summary_statistics()
    
    # Save summary statistics
    summary_file = OUTPUT_DIR / "summary_statistics.json"
    with open(summary_file, "w") as f:
        json.dump(summary_stats, f, indent=4, default=to_serializable)
    
    logger.info(f"Saved summary statistics to {summary_file}")
    
    # ------------------------------------------------------------
    # GENERATE PLOTS
    # ------------------------------------------------------------
    logger.info("Generating plots...")
    
    try:
        analyzer.plot_steroid_distribution(FIGURES_DIR)
    except Exception as e:
        logger.error(f"Error plotting steroid distribution: {e}")
    
    try:
        analyzer.plot_route_analysis(FIGURES_DIR)
    except Exception as e:
        logger.error(f"Error plotting route analysis: {e}")
    
    try:
        analyzer.plot_duration_analysis(FIGURES_DIR)
    except Exception as e:
        logger.error(f"Error plotting duration analysis: {e}")
    
    logger.info("All tasks completed successfully.")


if __name__ == "__main__":
    main_steroids()