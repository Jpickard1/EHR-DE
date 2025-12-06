#!/usr/bin/env python3
"""
MIMIC-IV Temporal Analysis for RL Feature Engineering
=====================================================
Deep temporal analysis of patient trajectories including:
- Event frequency and timing patterns
- Medication patterns and dosing frequencies
- Measurement density and regularity
- Clinical state transitions
- Time-to-event analysis

Author: Research Software Engineering Team
Date: 2025-12-05
Version: 1.0.0
"""

import os
import json
import pandas as pd
import numpy as np
from pathlib import Path
from datetime import datetime, timedelta
from tqdm import tqdm
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict, Counter
from multiprocessing import Pool, cpu_count
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

from utils import identify_sepsis_bool

# ============================================================================
# CONFIGURATION
# ============================================================================
from config import *

OUTPUT_DIR = EDA_DIR / f"time_analysis"
FIGURES_DIR = OUTPUT_DIR / 'figures'
CACHE_DIR = OUTPUT_DIR / 'cache'

# Create directories if they don't exist
for d in [OUTPUT_DIR, FIGURES_DIR, CACHE_DIR]:
    d.mkdir(parents=True, exist_ok=True)

# Plotting
sns.set_style("whitegrid")
plt.rcParams['figure.figsize'] = (14, 6)


# ============================================================================
# TEMPORAL ANALYSIS FUNCTIONS
# ============================================================================

def extract_patient_timeline(patient_file: Path) -> dict:
    """
    Comprehensive temporal analysis of a single patient's trajectory.
    
    Returns detailed timing statistics suitable for RL feature engineering.
    """
    try:
        with open(patient_file, 'r') as f:
            data = json.load(f)
        
        subject_id = data.get('subject_id')
        has_sepsis = identify_sepsis_bool(data)
        demographics = data.get('demographics', {})
        records = data.get('records', [])
        
        # Initialize result structure
        result = {
            'subject_id': subject_id,
            'has_sepsis': has_sepsis,
            'anchor_age': demographics.get('anchor_age'),
            'gender': demographics.get('gender'),
            'died': demographics.get('dod') is not None,
        }
        
        # ====================================================================
        # ADMISSION-LEVEL TEMPORAL ANALYSIS
        # ====================================================================
        
        admissions = [r for r in records if r.get('data_type') == 'admissions']
        icu_stays = [r for r in records if r.get('data_type') == 'icu_stays']
        
        result['num_admissions'] = len(admissions)
        result['num_icu_stays'] = len(icu_stays)
        
        # Hospital lengths of stay
        hospital_los = []
        admission_times = []
        
        for adm in admissions:
            admit_time = pd.to_datetime(adm.get('admittime'), errors='coerce')
            disch_time = pd.to_datetime(adm.get('dischtime'), errors='coerce')
            
            if pd.notna(admit_time) and pd.notna(disch_time):
                los_hours = (disch_time - admit_time).total_seconds() / 3600
                hospital_los.append(los_hours)
                admission_times.append(admit_time)
        
        result['hospital_los_hours_mean'] = np.mean(hospital_los) if hospital_los else None
        result['hospital_los_hours_median'] = np.median(hospital_los) if hospital_los else None
        result['hospital_los_hours_std'] = np.std(hospital_los) if hospital_los else None
        result['hospital_los_hours_max'] = np.max(hospital_los) if hospital_los else None
        result['total_hospital_hours'] = np.sum(hospital_los) if hospital_los else 0
        
        # ICU lengths of stay
        icu_los = []
        icu_admission_times = []
        
        for icu in icu_stays:
            in_time = pd.to_datetime(icu.get('intime'), errors='coerce')
            out_time = pd.to_datetime(icu.get('outtime'), errors='coerce')
            
            if pd.notna(in_time) and pd.notna(out_time):
                los_hours = (out_time - in_time).total_seconds() / 3600
                icu_los.append(los_hours)
                icu_admission_times.append(in_time)
        
        result['icu_los_hours_mean'] = np.mean(icu_los) if icu_los else None
        result['icu_los_hours_median'] = np.median(icu_los) if icu_los else None
        result['icu_los_hours_std'] = np.std(icu_los) if icu_los else None
        result['icu_los_hours_max'] = np.max(icu_los) if icu_los else None
        result['total_icu_hours'] = np.sum(icu_los) if icu_los else 0
        
        # Time between admissions (readmission patterns)
        if len(admission_times) > 1:
            admission_times_sorted = sorted(admission_times)
            gaps = [(admission_times_sorted[i+1] - admission_times_sorted[i]).total_seconds() / 86400 
                    for i in range(len(admission_times_sorted) - 1)]
            result['readmission_gap_days_mean'] = np.mean(gaps)
            result['readmission_gap_days_min'] = np.min(gaps)
            result['has_30day_readmission'] = any(g < 30 for g in gaps)
        else:
            result['readmission_gap_days_mean'] = None
            result['readmission_gap_days_min'] = None
            result['has_30day_readmission'] = False
        
        # ====================================================================
        # EVENT FREQUENCY ANALYSIS
        # ====================================================================
        
        # Get all timed events
        timed_events = [r for r in records if r.get('time')]
        
        if not timed_events:
            # Fill with None/0 for patients with no timed events
            result.update({
                'total_timed_events': 0,
                'timeline_hours': None,
                'events_per_hour': None,
            })
            return result
        
        # Convert times to datetime
        event_times = pd.to_datetime([e['time'] for e in timed_events], errors='coerce')
        event_times = event_times.dropna().sort_values()
        
        if len(event_times) == 0:
            result['total_timed_events'] = 0
            result['timeline_hours'] = None
            result['events_per_hour'] = None
            return result
        
        first_time = event_times.min()
        last_time = event_times.max()
        timeline_hours = (last_time - first_time).total_seconds() / 3600
        
        result['first_event_time'] = first_time.isoformat()
        result['last_event_time'] = last_time.isoformat()
        result['timeline_hours'] = timeline_hours
        result['total_timed_events'] = len(timed_events)
        result['events_per_hour'] = len(timed_events) / max(timeline_hours, 1)
        
        # ====================================================================
        # EVENT TYPE-SPECIFIC FREQUENCIES
        # ====================================================================
        
        event_type_times = defaultdict(list)
        
        for event in timed_events:
            event_type = event.get('event_type', event.get('data_type', 'unknown'))
            event_time = pd.to_datetime(event.get('time'), errors='coerce')
            if pd.notna(event_time):
                event_type_times[event_type].append(event_time)
        
        # Analyze each event type
        for event_type, times in event_type_times.items():
            times_sorted = sorted(times)
            count = len(times_sorted)
            
            # Frequency metrics
            result[f'{event_type}_count'] = count
            result[f'{event_type}_per_hour'] = count / max(timeline_hours, 1)
            
            # Inter-event intervals (sampling regularity)
            if len(times_sorted) > 1:
                intervals = [(times_sorted[i+1] - times_sorted[i]).total_seconds() / 60 
                            for i in range(len(times_sorted) - 1)]
                
                result[f'{event_type}_interval_minutes_mean'] = np.mean(intervals)
                result[f'{event_type}_interval_minutes_median'] = np.median(intervals)
                result[f'{event_type}_interval_minutes_std'] = np.std(intervals)
                result[f'{event_type}_interval_minutes_min'] = np.min(intervals)
                result[f'{event_type}_interval_minutes_max'] = np.max(intervals)
                
                # Coefficient of variation (regularity measure)
                cv = np.std(intervals) / max(np.mean(intervals), 0.001)
                result[f'{event_type}_interval_cv'] = cv
                
                # Percentage of regular measurements (within 2x median interval)
                median_interval = np.median(intervals)
                regular_count = sum(1 for i in intervals if i < 2 * median_interval)
                result[f'{event_type}_regularity_pct'] = regular_count / len(intervals) * 100
            else:
                result[f'{event_type}_interval_minutes_mean'] = None
                result[f'{event_type}_interval_minutes_median'] = None
                result[f'{event_type}_interval_cv'] = None
                result[f'{event_type}_regularity_pct'] = None
        
        # ====================================================================
        # PRESCRIPTION PATTERN ANALYSIS
        # ====================================================================
        
        prescriptions = [e for e in timed_events if e.get('event_type') == 'prescription']
        
        result['total_prescriptions'] = len(prescriptions)
        
        if prescriptions:
            # Unique medications
            unique_drugs = set()
            drug_frequencies = Counter()
            drug_durations = defaultdict(list)
            
            for rx in prescriptions:
                drug = rx.get('drug', '').lower().strip()
                if drug:
                    unique_drugs.add(drug)
                    drug_frequencies[drug] += 1
                    
                    # Calculate duration if available
                    start = pd.to_datetime(rx.get('starttime'), errors='coerce')
                    stop = pd.to_datetime(rx.get('stoptime'), errors='coerce')
                    if pd.notna(start) and pd.notna(stop):
                        duration_hours = (stop - start).total_seconds() / 3600
                        drug_durations[drug].append(duration_hours)
            
            result['unique_medications'] = len(unique_drugs)
            result['prescriptions_per_hour'] = len(prescriptions) / max(timeline_hours, 1)
            
            # Most common medications
            if drug_frequencies:
                top_3_drugs = drug_frequencies.most_common(3)
                result['most_common_medication'] = top_3_drugs[0][0] if top_3_drugs else None
                result['most_common_medication_count'] = top_3_drugs[0][1] if top_3_drugs else 0
            
            # Medication duration statistics
            all_durations = [d for durations in drug_durations.values() for d in durations]
            if all_durations:
                result['medication_duration_hours_mean'] = np.mean(all_durations)
                result['medication_duration_hours_median'] = np.median(all_durations)
                result['medication_duration_hours_max'] = np.max(all_durations)
            
            # Polypharmacy indicator (5+ concurrent medications)
            prescription_times = pd.to_datetime([p.get('time') for p in prescriptions], errors='coerce').dropna()
            if len(prescription_times) > 0:
                # Sample time windows and count concurrent meds
                max_concurrent = 0
                for t in prescription_times[::max(len(prescription_times)//10, 1)]:  # Sample 10 timepoints
                    concurrent = sum(1 for p in prescriptions 
                                   if pd.to_datetime(p.get('starttime'), errors='coerce') <= t 
                                   and (pd.isna(pd.to_datetime(p.get('stoptime'), errors='coerce')) 
                                        or pd.to_datetime(p.get('stoptime'), errors='coerce') >= t))
                    max_concurrent = max(max_concurrent, concurrent)
                result['max_concurrent_medications'] = max_concurrent
                result['has_polypharmacy'] = max_concurrent >= 5
            else:
                result['max_concurrent_medications'] = 0
                result['has_polypharmacy'] = False
        else:
            result['unique_medications'] = 0
            result['prescriptions_per_hour'] = 0
            result['max_concurrent_medications'] = 0
            result['has_polypharmacy'] = False
        
        # ====================================================================
        # LAB EVENT TEMPORAL PATTERNS
        # ====================================================================
        
        lab_events = [e for e in timed_events if e.get('event_type') == 'lab']
        
        if lab_events:
            # Unique lab tests
            unique_labs = set(e.get('itemid') for e in lab_events if e.get('itemid'))
            result['unique_lab_tests'] = len(unique_labs)
            result['labs_per_hour'] = len(lab_events) / max(timeline_hours, 1)
            
            # Panel ordering patterns (multiple labs at same time)
            lab_times = [pd.to_datetime(e.get('time'), errors='coerce') for e in lab_events]
            lab_times = [t for t in lab_times if pd.notna(t)]
            
            if lab_times:
                time_counts = Counter([t.floor('H') for t in lab_times])  # Count by hour
                result['max_labs_per_hour'] = max(time_counts.values())
                result['avg_labs_per_panel'] = np.mean(list(time_counts.values()))
        else:
            result['unique_lab_tests'] = 0
            result['labs_per_hour'] = 0
        
        # ====================================================================
        # VITAL SIGNS (CHART EVENTS) TEMPORAL PATTERNS
        # ====================================================================
        
        chart_events = [e for e in timed_events if e.get('event_type') == 'chart']
        
        if chart_events:
            result['chart_events_count'] = len(chart_events)
            result['charts_per_hour'] = len(chart_events) / max(timeline_hours, 1)
            
            # Identify vital sign patterns (high frequency measurements)
            chart_times = [pd.to_datetime(e.get('time'), errors='coerce') for e in chart_events]
            chart_times = sorted([t for t in chart_times if pd.notna(t)])
            
            if len(chart_times) > 1:
                # Detect monitoring intensity periods
                intervals = [(chart_times[i+1] - chart_times[i]).total_seconds() / 60 
                            for i in range(len(chart_times) - 1)]
                
                # High-frequency monitoring (< 30 min intervals)
                intensive_monitoring_pct = sum(1 for i in intervals if i < 30) / len(intervals) * 100
                result['intensive_monitoring_pct'] = intensive_monitoring_pct
        else:
            result['chart_events_count'] = 0
            result['charts_per_hour'] = 0
            result['intensive_monitoring_pct'] = 0
        
        # ====================================================================
        # INPUT/OUTPUT EVENT PATTERNS
        # ====================================================================
        
        input_events = [e for e in timed_events if e.get('event_type') == 'input']
        output_events = [e for e in timed_events if e.get('event_type') == 'output']
        
        result['input_events_count'] = len(input_events)
        result['output_events_count'] = len(output_events)
        result['inputs_per_hour'] = len(input_events) / max(timeline_hours, 1) if input_events else 0
        result['outputs_per_hour'] = len(output_events) / max(timeline_hours, 1) if output_events else 0
        
        # Fluid balance monitoring intensity
        result['fluid_monitoring_ratio'] = (len(input_events) + len(output_events)) / max(timeline_hours, 1)
        
        # ====================================================================
        # MICROBIOLOGY EVENT PATTERNS
        # ====================================================================
        
        micro_events = [e for e in timed_events if e.get('event_type') == 'microbiology']
        
        result['microbiology_events_count'] = len(micro_events)
        result['micro_per_admission'] = len(micro_events) / max(len(admissions), 1)
        
        # ====================================================================
        # TEMPORAL CLUSTERING (Early vs Late Events)
        # ====================================================================
        
        if len(event_times) > 0:
            # Calculate relative time from first event
            relative_times = [(t - first_time).total_seconds() / 3600 for t in event_times]
            
            # Divide timeline into quartiles
            q1, q2, q3 = np.percentile(relative_times, [25, 50, 75])
            
            events_q1 = sum(1 for t in relative_times if t <= q1)
            events_q2 = sum(1 for t in relative_times if q1 < t <= q2)
            events_q3 = sum(1 for t in relative_times if q2 < t <= q3)
            events_q4 = sum(1 for t in relative_times if t > q3)
            
            # Early loading (more events in first quartile)
            result['events_first_quartile_pct'] = events_q1 / len(relative_times) * 100
            result['events_last_quartile_pct'] = events_q4 / len(relative_times) * 100
            result['temporal_skew'] = (events_q1 - events_q4) / len(relative_times)  # Positive = front-loaded
        
        # ====================================================================
        # DATA COMPLETENESS AND QUALITY METRICS
        # ====================================================================
        
        # Percentage of events with timestamps
        total_records = len(records)
        result['temporal_coverage_pct'] = len(timed_events) / max(total_records, 1) * 100
        
        # Missing data patterns
        missing_count = sum(1 for r in records if r.get('time') is None)
        result['missing_timestamps_count'] = missing_count
        result['missing_timestamps_pct'] = missing_count / max(total_records, 1) * 100
        
        return result
        
    except Exception as e:
        return {
            'subject_id': None,
            'error': str(e)
        }

def extract_patient_timeline_wrapper(patient_file: Path) -> dict:
    """Wrapper for parallel processing."""
    return extract_patient_timeline(patient_file)

# ============================================================================
# TEMPORAL ANALYSIS VISUALIZATIONS
# ============================================================================

class TemporalAnalyzer:
    """Analyze and visualize temporal patterns in the cohort."""
    
    def __init__(self, temporal_df: pd.DataFrame):
        self.df = temporal_df
        self.sepsis_df = self.df[self.df['has_sepsis'] == True]
        self.non_sepsis_df = self.df[self.df['has_sepsis'] == False]
    
    def plot_event_frequency_distributions(self, output_dir: Path):
        """Plot event frequency distributions by type."""
        
        fig, axes = plt.subplots(2, 3, figsize=(18, 10))
        
        event_types = [
            ('lab', 'Lab Events'),
            ('chart', 'Chart Events'),
            ('prescription', 'Prescriptions'),
            ('input', 'Input Events'),
            ('output', 'Output Events'),
            ('microbiology', 'Microbiology')
        ]
        
        for idx, (event_type, title) in enumerate(event_types):
            row, col = idx // 3, idx % 3
            ax = axes[row, col]
            
            col_name = f'{event_type}_per_hour'
            if col_name in self.df.columns:
                data_all = self.df[col_name].dropna()
                data_sepsis = self.sepsis_df[col_name].dropna()
                data_non_sepsis = self.non_sepsis_df[col_name].dropna()
                
                if len(data_all) > 0:
                    # Plot distributions
                    ax.hist(data_non_sepsis, bins=50, alpha=0.6, label='Non-Sepsis',
                           color='blue', edgecolor='black', range=(0, np.percentile(data_all, 95)))
                    ax.hist(data_sepsis, bins=50, alpha=0.6, label='Sepsis',
                           color='red', edgecolor='black', range=(0, np.percentile(data_all, 95)))
                    
                    ax.set_xlabel('Events per Hour')
                    ax.set_ylabel('Number of Patients')
                    ax.set_title(f'{title} Frequency')
                    ax.legend()
                    ax.grid(alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(output_dir / 'event_frequency_distributions.png', dpi=150, bbox_inches='tight')
        plt.close()
    
    def plot_temporal_patterns(self, output_dir: Path):
        """Plot various temporal patterns."""
        
        fig, axes = plt.subplots(2, 2, figsize=(14, 10))
        
        # 1. Length of stay comparison
        ax = axes[0, 0]
        los_data = [
            self.non_sepsis_df['hospital_los_hours_mean'].dropna() / 24,
            self.sepsis_df['hospital_los_hours_mean'].dropna() / 24
        ]
        bp = ax.boxplot(los_data, labels=['Non-Sepsis', 'Sepsis'], patch_artist=True)
        bp['boxes'][0].set_facecolor('lightblue')
        bp['boxes'][1].set_facecolor('lightcoral')
        ax.set_ylabel('Hospital Length of Stay (days)')
        ax.set_title('Hospital LOS: Sepsis vs Non-Sepsis')
        ax.set_ylim(0, 30)
        ax.grid(alpha=0.3)
        
        # 2. ICU length of stay
        ax = axes[0, 1]
        icu_los_data = [
            self.non_sepsis_df['icu_los_hours_mean'].dropna(),
            self.sepsis_df['icu_los_hours_mean'].dropna()
        ]
        bp = ax.boxplot(icu_los_data, labels=['Non-Sepsis', 'Sepsis'], patch_artist=True)
        bp['boxes'][0].set_facecolor('lightblue')
        bp['boxes'][1].set_facecolor('lightcoral')
        ax.set_ylabel('ICU Length of Stay (hours)')
        ax.set_title('ICU LOS: Sepsis vs Non-Sepsis')
        ax.set_ylim(0, 500)
        ax.grid(alpha=0.3)
        
        # 3. Medication frequency
        ax = axes[1, 0]
        med_data = [
            self.non_sepsis_df['unique_medications'].dropna(),
            self.sepsis_df['unique_medications'].dropna()
        ]
        bp = ax.boxplot(med_data, labels=['Non-Sepsis', 'Sepsis'], patch_artist=True)
        bp['boxes'][0].set_facecolor('lightblue')
        bp['boxes'][1].set_facecolor('lightcoral')
        ax.set_ylabel('Number of Unique Medications')
        ax.set_title('Medication Complexity')
        ax.grid(alpha=0.3)
        
        # 4. Data density (events per hour)
        ax = axes[1, 1]
        density_data = [
            self.non_sepsis_df['events_per_hour'].dropna(),
            self.sepsis_df['events_per_hour'].dropna()
        ]
        bp = ax.boxplot(density_data, labels=['Non-Sepsis', 'Sepsis'], patch_artist=True)
        bp['boxes'][0].set_facecolor('lightblue')
        bp['boxes'][1].set_facecolor('lightcoral')
        ax.set_ylabel('Events per Hour')
        ax.set_title('Data Density')
        ax.grid(alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(output_dir / 'temporal_patterns_comparison.png', dpi=150, bbox_inches='tight')
        plt.close()
    
    def plot_measurement_regularity(self, output_dir: Path):
        """Plot measurement interval regularity."""
        
        fig, axes = plt.subplots(2, 2, figsize=(14, 10))
        
        event_types = ['lab', 'chart', 'prescription', 'input']
        titles = ['Lab Tests', 'Vital Signs', 'Prescriptions', 'Fluid Inputs']
        
        for idx, (event_type, title) in enumerate(zip(event_types, titles)):
            row, col = idx // 2, idx % 2
            ax = axes[row, col]
            
            col_name = f'{event_type}_interval_minutes_median'
            if col_name in self.df.columns:
                data_non_sepsis = self.non_sepsis_df[col_name].dropna()
                data_sepsis = self.sepsis_df[col_name].dropna()
                
                if len(data_non_sepsis) > 0 or len(data_sepsis) > 0:
                    bins = np.logspace(0, 4, 50)  # Log scale bins
                    
                    ax.hist(data_non_sepsis, bins=bins, alpha=0.6, label='Non-Sepsis',
                           color='blue', edgecolor='black')
                    ax.hist(data_sepsis, bins=bins, alpha=0.6, label='Sepsis',
                           color='red', edgecolor='black')
                    
                    ax.set_xscale('log')
                    ax.set_xlabel('Median Interval Between Measurements (minutes)')
                    ax.set_ylabel('Number of Patients')
                    ax.set_title(f'{title} Measurement Intervals')
                    ax.legend()
                    ax.grid(alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(output_dir / 'measurement_regularity.png', dpi=150, bbox_inches='tight')
        plt.close()
    
    def generate_summary_table(self) -> pd.DataFrame:
        """Generate summary statistics table."""
        
        metrics = {
            'Hospital LOS (days)': ('hospital_los_hours_mean', lambda x: x/24),
            'ICU LOS (hours)': ('icu_los_hours_mean', lambda x: x),
            'Unique Medications': ('unique_medications', lambda x: x),
            'Events per Hour': ('events_per_hour', lambda x: x),
            'Labs per Hour': ('lab_per_hour', lambda x: x),
            'Charts per Hour': ('chart_per_hour', lambda x: x),
            'Timeline (hours)': ('timeline_hours', lambda x: x),
            'Prescriptions per Hour': ('prescriptions_per_hour', lambda x: x),
        }
        
        summary_data = []
        
        for metric_name, (col_name, transform) in metrics.items():
            if col_name in self.df.columns:
                all_data = self.df[col_name].dropna().apply(transform)
                sepsis_data = self.sepsis_df[col_name].dropna().apply(transform)
                non_sepsis_data = self.non_sepsis_df[col_name].dropna().apply(transform)
                
                summary_data.append({
                    'Metric': metric_name,
                    'All Mean': f"{all_data.mean():.2f}",
                    'All Median': f"{all_data.median():.2f}",
                    'Sepsis Mean': f"{sepsis_data.mean():.2f}",
                    'Non-Sepsis Mean': f"{non_sepsis_data.mean():.2f}",
                    'Ratio (S/NS)': f"{sepsis_data.mean() / max(non_sepsis_data.mean(), 0.001):.2f}x"
                })
        
        return pd.DataFrame(summary_data)

# ============================================================================
# MAIN EXECUTION
# ============================================================================

def main_time():
    """Main temporal analysis pipeline."""
    
    print("=" * 80)
    print("MIMIC-IV TEMPORAL PATTERN ANALYSIS")
    print("=" * 80)
    print(f"Start time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Input directory: {CLEANED_DIR}")
    print(f"Output directory: {OUTPUT_DIR}")
    print(f"Workers: {NUM_WORKERS}")
    print("=" * 80)
    print()
    
    import time
    overall_start = time.time()
    
    # ========================================================================
    # STEP 1: EXTRACT TEMPORAL PATTERNS
    # ========================================================================
    
    cache_file = CACHE_DIR / 'temporal_statistics.parquet'
    
    if cache_file.exists():
        print("Loading cached temporal statistics...")
        temporal_df = pd.read_parquet(cache_file)
        print(f"✓ Loaded {len(temporal_df):,} patient records from cache")
    else:
        print("STEP 1: Extracting temporal patterns from patient files...")
        
        # Find all patient files
        patient_files = sorted(CLEANED_DIR.glob('patient_*.json'))
        print(f"Found {len(patient_files):,} patient files")
        
        # Sample if specified
        if SAMPLE_SIZE and SAMPLE_SIZE < len(patient_files):
            patient_files = patient_files[:SAMPLE_SIZE]
            print(f"Analyzing sample of {len(patient_files):,} patients")
        
        # Parallel analysis
        print(f"Analyzing with {NUM_WORKERS} workers...")
        with Pool(NUM_WORKERS) as pool:
            results = list(tqdm(
                pool.imap(extract_patient_timeline_wrapper, patient_files, chunksize=10),
                total=len(patient_files),
                desc="Extracting timelines"
            ))
        
        # Convert to DataFrame
        temporal_df = pd.DataFrame([r for r in results if r.get('subject_id') is not None])
        
        # Save cache
        temporal_df.to_parquet(cache_file, index=False)
        print(f"✓ Cached temporal statistics to {cache_file}")
    
    print(f"\n✓ Analyzed {len(temporal_df):,} patients")
    print()
    
    # ========================================================================
    # STEP 2: GENERATE VISUALIZATIONS
    # ========================================================================
    
    print("STEP 2: Generating temporal visualizations...")
    
    analyzer = TemporalAnalyzer(temporal_df)
    
    analyzer.plot_event_frequency_distributions(FIGURES_DIR)
    print("  ✓ Event frequency distributions")
    
    analyzer.plot_temporal_patterns(FIGURES_DIR)
    print("  ✓ Temporal patterns comparison")
    
    analyzer.plot_measurement_regularity(FIGURES_DIR)
    print("  ✓ Measurement regularity analysis")
    
    print(f"\n✓ Saved figures to {FIGURES_DIR}")
    
    # ========================================================================
    # STEP 3: GENERATE SUMMARY STATISTICS
    # ========================================================================
    
    print("\nSTEP 3: Generating summary statistics...")
    
    summary_table = analyzer.generate_summary_table()
    summary_file = OUTPUT_DIR / 'temporal_summary_table.csv'
    summary_table.to_csv(summary_file, index=False)
    print(f"✓ Saved summary table to {summary_file}")
    
    # Save full temporal dataset
    full_output = OUTPUT_DIR / 'temporal_features_full.parquet'
    temporal_df.to_parquet(full_output, index=False)
    print(f"✓ Saved full temporal dataset to {full_output}")
    
    # ========================================================================
    # STEP 4: PRINT SUMMARY REPORT
    # ========================================================================
    
    print("\n" + "=" * 80)
    print("TEMPORAL ANALYSIS SUMMARY")
    print("=" * 80)
    
    sepsis_df = temporal_df[temporal_df['has_sepsis'] == True]
    non_sepsis_df = temporal_df[temporal_df['has_sepsis'] == False]
    
    print(f"\n📊 COHORT SIZE")
    print(f"  Total patients:           {len(temporal_df):,}")
    print(f"  Sepsis patients:          {len(sepsis_df):,}")
    print(f"  Non-sepsis patients:      {len(non_sepsis_df):,}")
    
    print(f"\n⏱️  LENGTH OF STAY")
    print(f"  Hospital LOS (Sepsis):    {sepsis_df['hospital_los_hours_mean'].mean()/24:.1f} days")
    print(f"  Hospital LOS (Non-Sepsis):{non_sepsis_df['hospital_los_hours_mean'].mean()/24:.1f} days")
    print(f"  ICU LOS (Sepsis):         {sepsis_df['icu_los_hours_mean'].mean():.1f} hours")
    print(f"  ICU LOS (Non-Sepsis):     {non_sepsis_df['icu_los_hours_mean'].mean():.1f} hours")
    
    print(f"\n📈 EVENT FREQUENCIES")
    print(f"  Events/hour (Sepsis):     {sepsis_df['events_per_hour'].mean():.1f}")
    print(f"  Events/hour (Non-Sepsis): {non_sepsis_df['events_per_hour'].mean():.1f}")
    print(f"  Labs/hour (Sepsis):       {sepsis_df['lab_per_hour'].mean():.2f}")
    print(f"  Labs/hour (Non-Sepsis):   {non_sepsis_df['lab_per_hour'].mean():.2f}")
    
    print(f"\n💊 MEDICATION PATTERNS")
    print(f"  Unique meds (Sepsis):     {sepsis_df['unique_medications'].mean():.1f}")
    print(f"  Unique meds (Non-Sepsis): {non_sepsis_df['unique_medications'].mean():.1f}")
    print(f"  Polypharmacy rate (Sepsis):{sepsis_df['has_polypharmacy'].mean()*100:.1f}%")
    print(f"  Polypharmacy (Non-Sepsis):{non_sepsis_df['has_polypharmacy'].mean()*100:.1f}%")
    
    print(f"\n🔬 MEASUREMENT PATTERNS")
    if 'lab_interval_minutes_median' in temporal_df.columns:
        print(f"  Lab interval (Sepsis):    {sepsis_df['lab_interval_minutes_median'].median():.0f} min")
        print(f"  Lab interval (Non-Sepsis):{non_sepsis_df['lab_interval_minutes_median'].median():.0f} min")
    if 'chart_interval_minutes_median' in temporal_df.columns:
        print(f"  Vital sign interval (S):  {sepsis_df['chart_interval_minutes_median'].median():.0f} min")
        print(f"  Vital sign interval (NS): {non_sepsis_df['chart_interval_minutes_median'].median():.0f} min")
    
    print(f"\n📊 DATA QUALITY")
    print(f"  Temporal coverage:        {temporal_df['temporal_coverage_pct'].mean():.1f}%")
    print(f"  Timeline duration:        {temporal_df['timeline_hours'].mean():.1f} hours")
    print(f"  Intensive monitoring (S): {sepsis_df['intensive_monitoring_pct'].mean():.1f}%")
    print(f"  Intensive monitoring (NS):{non_sepsis_df['intensive_monitoring_pct'].mean():.1f}%")
    
    print("\n" + "=" * 80)
    print("\n📋 SUMMARY TABLE")
    print("=" * 80)
    print(summary_table.to_string(index=False))
    
    print("\n" + "=" * 80)
    elapsed = time.time() - overall_start
    print(f"Analysis completed in {elapsed/60:.1f} minutes")
    print(f"Results saved to: {OUTPUT_DIR}")
    print("=" * 80)
    print()

if __name__ == '__main__':
    main_time()
