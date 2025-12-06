#!/usr/bin/env python3
"""
MIMIC-IV Sepsis-Focused Exploratory Data Analysis
==================================================
Comprehensive analysis of patient trajectories with sepsis identification
and dataset statistics for RL preparation.
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
warnings.filterwarnings('ignore')

from clinical_constants import SEPSIS_ICD10_CODES, SEPSIS_ICD9_CODES, ORGAN_DYSFUNCTION_ICD10
from utils import identify_sepsis_detailed

# ============================================================================
# CONFIGURATION
# ============================================================================
from config import *

OUTPUT_DIR = EDA_DIR / f"sepsis"
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
        logging.FileHandler(OUTPUT_DIR / f'eda_{datetime.now().strftime("%Y%m%d_%H%M%S")}.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# Plotting style
sns.set_style("whitegrid")
plt.rcParams['figure.figsize'] = (12, 6)
plt.rcParams['font.size'] = 10

def is_valid_date(date_value):
    """Check if a date value is valid (not None, not 'nan', not empty)."""
    if date_value is None:
        return False
    if isinstance(date_value, str):
        date_str = date_value.strip().lower()
        if date_str in ['', 'nan', 'none', 'null', 'nat']:
            return False
    if pd.isna(date_value):
        return False
    return True

# ============================================================================
# PATIENT-LEVEL ANALYSIS
# ============================================================================

def analyze_patient(patient_file: Path) -> dict:
    """
    Analyze a single patient file and extract comprehensive statistics.
    
    Returns:
        Dictionary with all patient-level metrics
    """
    try:
        with open(patient_file, 'r') as f:
            data = json.load(f)
        
        subject_id = data.get('subject_id')
        demographics = data.get('demographics', {})
        if not isinstance(demographics, dict):
            demographics = {}
            
        records = data.get('records', [])
        if not isinstance(records, list):
            records = []
            
        metadata = data.get('metadata', {})
        if not isinstance(metadata, dict):
            metadata = {}
        
        # Basic demographics
        stats = {
            'subject_id': subject_id,
            'gender': demographics.get('gender'),
            'anchor_age': demographics.get('anchor_age'),
            'dod': demographics.get('dod'),
        }
        
        # FIX: Properly check if patient died
        dod = demographics.get('dod')
        stats['died'] = is_valid_date(dod)
        
        # Sepsis identification
        sepsis_info = identify_sepsis_detailed(data)
        stats.update({
            'has_sepsis': sepsis_info['has_sepsis'],
            'sepsis_type': sepsis_info['sepsis_type'],
            'num_sepsis_admissions': len(sepsis_info['sepsis_admissions']),
            'sepsis_codes': sepsis_info['sepsis_codes'],
            'organ_dysfunctions': sepsis_info['organ_dysfunctions'],
        })
        
        # Categorize records by type
        data_types = defaultdict(list)
        for record in records:
            if not isinstance(record, dict):
                continue
            dtype = record.get('data_type', 'unknown')
            data_types[dtype].append(record)
        
        # Admission-level statistics
        admissions = data_types.get('admissions', [])
        icu_stays = data_types.get('icu_stays', [])
        
        stats['num_admissions'] = len(admissions)
        stats['num_icu_stays'] = len(icu_stays)
        
        # Hospital lengths of stay
        los_list = []
        for adm in admissions:
            if not isinstance(adm, dict):
                continue
            admit_time = pd.to_datetime(adm.get('admittime'), errors='coerce')
            disch_time = pd.to_datetime(adm.get('dischtime'), errors='coerce')
            if pd.notna(admit_time) and pd.notna(disch_time):
                los_days = (disch_time - admit_time).total_seconds() / 86400
                if los_days >= 0:  # Sanity check
                    los_list.append(los_days)
        
        stats['mean_los_days'] = np.mean(los_list) if los_list else None
        stats['total_hospital_days'] = sum(los_list) if los_list else 0
        
        # ICU lengths of stay
        icu_los_list = []
        for icu in icu_stays:
            if not isinstance(icu, dict):
                continue
            in_time = pd.to_datetime(icu.get('intime'), errors='coerce')
            out_time = pd.to_datetime(icu.get('outtime'), errors='coerce')
            if pd.notna(in_time) and pd.notna(out_time):
                los_hours = (out_time - in_time).total_seconds() / 3600
                if los_hours >= 0:  # Sanity check
                    icu_los_list.append(los_hours)
        
        stats['mean_icu_los_hours'] = np.mean(icu_los_list) if icu_los_list else None
        stats['total_icu_hours'] = sum(icu_los_list) if icu_los_list else 0
        
        # Diagnoses
        diagnoses = data_types.get('diagnoses', [])
        stats['num_diagnoses'] = len(diagnoses)
        unique_diag_codes = set()
        for d in diagnoses:
            if isinstance(d, dict) and d.get('icd_code'):
                unique_diag_codes.add(d['icd_code'])
        stats['unique_diagnoses'] = len(unique_diag_codes)
        
        # Procedures
        procedures = data_types.get('procedures', [])
        stats['num_procedures'] = len(procedures)
        unique_proc_codes = set()
        for p in procedures:
            if isinstance(p, dict) and p.get('icd_code'):
                unique_proc_codes.add(p['icd_code'])
        stats['unique_procedures'] = len(unique_proc_codes)
        
        # Events (labs, vitals, medications, etc.)
        events = data_types.get('events', [])
        stats['num_events'] = len(events)
        
        # Break down by event type
        event_types = Counter()
        for e in events:
            if isinstance(e, dict) and e.get('event_type'):
                event_types[e['event_type']] += 1
        
        stats['num_lab_events'] = event_types.get('lab', 0)
        stats['num_chart_events'] = event_types.get('chart', 0)
        stats['num_input_events'] = event_types.get('input', 0)
        stats['num_output_events'] = event_types.get('output', 0)
        stats['num_prescription_events'] = event_types.get('prescription', 0)
        stats['num_microbiology_events'] = event_types.get('microbiology', 0)
        
        # Medications (from prescriptions)
        prescriptions = [e for e in events if isinstance(e, dict) and e.get('event_type') == 'prescription']
        unique_drugs = set()
        for rx in prescriptions:
            drug = rx.get('drug')
            if drug and isinstance(drug, str):
                unique_drugs.add(drug.lower().strip())
        stats['num_unique_medications'] = len(unique_drugs)
        
        # Timeline statistics
        timed_records = [r for r in records if isinstance(r, dict) and r.get('time')]
        if timed_records:
            times = pd.to_datetime([r['time'] for r in timed_records], errors='coerce')
            times = times.dropna()
            if len(times) > 0:
                stats['first_record_time'] = times.min().isoformat()
                stats['last_record_time'] = times.max().isoformat()
                timeline_seconds = (times.max() - times.min()).total_seconds()
                stats['timeline_days'] = timeline_seconds / 86400 if timeline_seconds > 0 else 0
            else:
                stats['timeline_days'] = None
        else:
            stats['timeline_days'] = None
        
        # Data density (events per day)
        if stats['timeline_days'] and stats['timeline_days'] > 0:
            stats['events_per_day'] = stats['num_events'] / stats['timeline_days']
        else:
            stats['events_per_day'] = None
        
        return stats
        
    except Exception as e:
        logger.error(f"Error analyzing {patient_file.name}: {e}")
        return {'subject_id': None, 'error': str(e)}

def analyze_patient_wrapper(patient_file: Path) -> dict:
    """Wrapper for parallel processing."""
    return analyze_patient(patient_file)

# ============================================================================
# COHORT-LEVEL ANALYSIS
# ============================================================================

class CohortAnalyzer:
    """Comprehensive cohort-level analysis."""
    
    def __init__(self, patient_stats_df: pd.DataFrame):
        self.df = patient_stats_df.copy()
        # Ensure died is boolean
        self.df['died'] = self.df['died'].fillna(False).astype(bool)
        self.df['has_sepsis'] = self.df['has_sepsis'].fillna(False).astype(bool)
        
        self.sepsis_df = self.df[self.df['has_sepsis'] == True]
        self.non_sepsis_df = self.df[self.df['has_sepsis'] == False]
        
    def generate_summary_statistics(self) -> dict:
        """Generate comprehensive summary statistics."""
        
        summary = {
            'cohort_overview': {
                'total_patients': len(self.df),
                'sepsis_patients': len(self.sepsis_df),
                'non_sepsis_patients': len(self.non_sepsis_df),
                'sepsis_prevalence': len(self.sepsis_df) / len(self.df) * 100 if len(self.df) > 0 else 0,
            },
            
            'demographics': {
                'age_mean': float(self.df['anchor_age'].mean()) if 'anchor_age' in self.df else None,
                'age_std': float(self.df['anchor_age'].std()) if 'anchor_age' in self.df else None,
                'age_median': float(self.df['anchor_age'].median()) if 'anchor_age' in self.df else None,
                'gender_distribution': self.df['gender'].value_counts().to_dict() if 'gender' in self.df else {},
                'mortality_rate': float(self.df['died'].sum() / len(self.df) * 100) if len(self.df) > 0 else 0,
            },
            
            'hospitalizations': {
                'admissions_per_patient_mean': float(self.df['num_admissions'].mean()),
                'admissions_per_patient_median': float(self.df['num_admissions'].median()),
                'icu_stays_per_patient_mean': float(self.df['num_icu_stays'].mean()),
                'icu_stays_per_patient_median': float(self.df['num_icu_stays'].median()),
                'mean_los_days': float(self.df['mean_los_days'].mean()),
                'mean_icu_los_hours': float(self.df['mean_icu_los_hours'].mean()),
            },
            
            'clinical_complexity': {
                'diagnoses_per_patient_mean': float(self.df['num_diagnoses'].mean()),
                'diagnoses_per_patient_median': float(self.df['num_diagnoses'].median()),
                'procedures_per_patient_mean': float(self.df['num_procedures'].mean()),
                'procedures_per_patient_median': float(self.df['num_procedures'].median()),
                'medications_per_patient_mean': float(self.df['num_unique_medications'].mean()),
                'medications_per_patient_median': float(self.df['num_unique_medications'].median()),
            },
            
            'data_characteristics': {
                'events_per_patient_mean': float(self.df['num_events'].mean()),
                'events_per_patient_median': float(self.df['num_events'].median()),
                'lab_events_per_patient': float(self.df['num_lab_events'].mean()),
                'chart_events_per_patient': float(self.df['num_chart_events'].mean()),
                'timeline_days_mean': float(self.df['timeline_days'].mean()),
                'events_per_day_mean': float(self.df['events_per_day'].mean()),
            },
            
            'sepsis_comparison': self._compare_sepsis_vs_non_sepsis(),
        }
        
        return summary
    
    def _compare_sepsis_vs_non_sepsis(self) -> dict:
        """Compare sepsis vs non-sepsis patients."""
        
        metrics = [
            'anchor_age', 'num_admissions', 'num_icu_stays',
            'mean_los_days', 'mean_icu_los_hours', 'num_diagnoses',
            'num_procedures', 'num_unique_medications', 'num_events'
        ]
        
        comparison = {}
        
        # Handle died separately since it's boolean
        if len(self.sepsis_df) > 0 and len(self.non_sepsis_df) > 0:
            comparison['died'] = {
                'sepsis_mean': float(self.sepsis_df['died'].mean()),
                'non_sepsis_mean': float(self.non_sepsis_df['died'].mean()),
                'difference': float(self.sepsis_df['died'].mean() - self.non_sepsis_df['died'].mean()),
                'ratio': float(self.sepsis_df['died'].mean() / max(self.non_sepsis_df['died'].mean(), 0.001)),
            }
        
        for metric in metrics:
            if metric in self.df.columns and len(self.sepsis_df) > 0 and len(self.non_sepsis_df) > 0:
                sepsis_mean = self.sepsis_df[metric].mean()
                non_sepsis_mean = self.non_sepsis_df[metric].mean()
                
                comparison[metric] = {
                    'sepsis_mean': float(sepsis_mean) if pd.notna(sepsis_mean) else None,
                    'non_sepsis_mean': float(non_sepsis_mean) if pd.notna(non_sepsis_mean) else None,
                    'difference': float(sepsis_mean - non_sepsis_mean) if pd.notna(sepsis_mean) and pd.notna(non_sepsis_mean) else None,
                    'ratio': float(sepsis_mean / max(non_sepsis_mean, 0.001)) if pd.notna(sepsis_mean) and pd.notna(non_sepsis_mean) else None,
                }
        
        return comparison
    
    def plot_distributions(self, output_dir: Path):
        """Generate comprehensive distribution plots."""
        
        metrics_to_plot = [
            ('anchor_age', 'Age Distribution', 'Age (years)', (0, 100)),
            ('num_admissions', 'Hospital Admissions per Patient', 'Number of Admissions', (0, 20)),
            ('num_icu_stays', 'ICU Stays per Patient', 'Number of ICU Stays', (0, 10)),
            ('num_diagnoses', 'Diagnoses per Patient', 'Number of Diagnoses', (0, 100)),
            ('num_procedures', 'Procedures per Patient', 'Number of Procedures', (0, 50)),
            ('num_unique_medications', 'Unique Medications per Patient', 'Number of Medications', (0, 100)),
            ('mean_los_days', 'Hospital Length of Stay', 'Days', (0, 30)),
            ('mean_icu_los_hours', 'ICU Length of Stay', 'Hours', (0, 500)),
        ]
        
        for metric, title, xlabel, xlim in metrics_to_plot:
            if metric not in self.df.columns:
                continue
                
            fig, axes = plt.subplots(1, 2, figsize=(14, 5))
            
            # Overall distribution
            data_to_plot = self.df[metric].dropna()
            if len(data_to_plot) > 0:
                axes[0].hist(data_to_plot, bins=50, alpha=0.7, edgecolor='black')
                axes[0].set_xlabel(xlabel)
                axes[0].set_ylabel('Number of Patients')
                axes[0].set_title(f'{title} - All Patients')
                axes[0].set_xlim(xlim)
                axes[0].grid(alpha=0.3)
            
            # Sepsis vs Non-sepsis comparison
            non_sepsis_data = self.non_sepsis_df[metric].dropna()
            sepsis_data = self.sepsis_df[metric].dropna()
            
            if len(non_sepsis_data) > 0:
                axes[1].hist(non_sepsis_data, bins=50, alpha=0.6, 
                            label='Non-Sepsis', edgecolor='black', color='blue')
            if len(sepsis_data) > 0:
                axes[1].hist(sepsis_data, bins=50, alpha=0.6, 
                            label='Sepsis', edgecolor='black', color='red')
            
            axes[1].set_xlabel(xlabel)
            axes[1].set_ylabel('Number of Patients')
            axes[1].set_title(f'{title} - Sepsis Comparison')
            axes[1].legend()
            axes[1].set_xlim(xlim)
            axes[1].grid(alpha=0.3)
            
            plt.tight_layout()
            plt.savefig(output_dir / f'dist_{metric}.png', dpi=150, bbox_inches='tight')
            plt.close()
        
        logger.info(f"Saved {len(metrics_to_plot)} distribution plots")
    
    def plot_sepsis_overview(self, output_dir: Path):
        """Generate sepsis-specific overview plots."""
        
        fig, axes = plt.subplots(2, 2, figsize=(14, 10))
        
        # Sepsis prevalence
        sepsis_counts = self.df['has_sepsis'].value_counts()
        labels = ['Non-Sepsis', 'Sepsis']
        axes[0, 0].pie(sepsis_counts, labels=labels, 
                       autopct='%1.1f%%', startangle=90, colors=['lightblue', 'coral'])
        axes[0, 0].set_title('Sepsis Prevalence in Cohort')
        
        # Sepsis type breakdown
        if len(self.sepsis_df) > 0:
            sepsis_types = self.sepsis_df['sepsis_type'].value_counts()
            if len(sepsis_types) > 0:
                axes[0, 1].bar(range(len(sepsis_types)), sepsis_types.values, color='coral')
                axes[0, 1].set_xticks(range(len(sepsis_types)))
                axes[0, 1].set_xticklabels(sepsis_types.index, rotation=45)
                axes[0, 1].set_ylabel('Number of Patients')
                axes[0, 1].set_title('Sepsis Identification Method')
                axes[0, 1].grid(alpha=0.3)
        
        # Mortality comparison
        mortality_data = pd.DataFrame({
            'Non-Sepsis': [self.non_sepsis_df['died'].sum(), len(self.non_sepsis_df) - self.non_sepsis_df['died'].sum()],
            'Sepsis': [self.sepsis_df['died'].sum(), len(self.sepsis_df) - self.sepsis_df['died'].sum()]
        }, index=['Died', 'Survived'])
        mortality_data.plot(kind='bar', ax=axes[1, 0], color=['lightblue', 'coral'])
        axes[1, 0].set_ylabel('Number of Patients')
        axes[1, 0].set_title('Mortality: Sepsis vs Non-Sepsis')
        axes[1, 0].set_xticklabels(['Died', 'Survived'], rotation=0)
        axes[1, 0].legend()
        axes[1, 0].grid(alpha=0.3)
        
        # Organ dysfunctions in sepsis patients
        all_dysfunctions = []
        for dysfunctions in self.sepsis_df['organ_dysfunctions'].dropna():
            if isinstance(dysfunctions, list):
                all_dysfunctions.extend(dysfunctions)
        
        if all_dysfunctions:
            dysfunction_counts = Counter(all_dysfunctions)
            axes[1, 1].barh(list(dysfunction_counts.keys()), list(dysfunction_counts.values()), color='coral')
            axes[1, 1].set_xlabel('Number of Patients')
            axes[1, 1].set_title('Organ Dysfunctions in Sepsis Patients')
            axes[1, 1].grid(alpha=0.3)
        else:
            axes[1, 1].text(0.5, 0.5, 'No organ dysfunction data', 
                           ha='center', va='center', transform=axes[1, 1].transAxes)
        
        plt.tight_layout()
        plt.savefig(output_dir / 'sepsis_overview.png', dpi=150, bbox_inches='tight')
        plt.close()
        
        logger.info("Saved sepsis overview plot")

# ============================================================================
# MAIN EXECUTION
# ============================================================================

def main_sepsis():
    """Main EDA pipeline."""
    
    print("=" * 80)
    print("MIMIC-IV SEPSIS-FOCUSED EXPLORATORY DATA ANALYSIS")
    print("=" * 80)
    print(f"Start time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Input directory: {CLEANED_DIR}")
    print(f"Output directory: {OUTPUT_DIR}")
    print(f"Workers: {NUM_WORKERS}")
    print("=" * 80)
    print()
    
    overall_start = time.time()
    
    # ========================================================================
    # STEP 1: COLLECT PATIENT-LEVEL STATISTICS
    # ========================================================================
    
    cache_file = CACHE_DIR / 'patient_statistics.parquet'
    
    if cache_file.exists() and False:
        print("Loading cached patient statistics...")
        patient_stats_df = pd.read_parquet(cache_file)
        print(f"✓ Loaded {len(patient_stats_df):,} patient records from cache")
    else:
        print("STEP 1: Analyzing individual patient files...")
        
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
                pool.imap(analyze_patient_wrapper, patient_files, chunksize=10),
                total=len(patient_files),
                desc="Analyzing patients"
            ))
        
        # Convert to DataFrame
        patient_stats_df = pd.DataFrame([r for r in results if r.get('subject_id') is not None])
        
        # Save cache
        patient_stats_df.to_parquet(cache_file, index=False)
        print(f"✓ Cached patient statistics to {cache_file}")
    
    print(f"\n✓ Analyzed {len(patient_stats_df):,} patients")
    print()
    
    # ========================================================================
    # STEP 2: COHORT-LEVEL ANALYSIS
    # ========================================================================
    
    print("STEP 2: Performing cohort-level analysis...")
    
    analyzer = CohortAnalyzer(patient_stats_df)
    summary = analyzer.generate_summary_statistics()
    
    # Save summary statistics
    summary_file = OUTPUT_DIR / 'cohort_summary.json'
    with open(summary_file, 'w') as f:
        json.dump(summary, f, indent=2, default=str)
    print(f"✓ Saved cohort summary to {summary_file}")
    
    # ========================================================================
    # STEP 3: GENERATE VISUALIZATIONS
    # ========================================================================
    
    print("\nSTEP 3: Generating visualizations...")
    
    analyzer.plot_distributions(FIGURES_DIR)
    analyzer.plot_sepsis_overview(FIGURES_DIR)
    
    print(f"✓ Saved figures to {FIGURES_DIR}")
    
    # ========================================================================
    # STEP 4: SEPSIS COHORT IDENTIFICATION
    # ========================================================================
    
    print("\nSTEP 4: Generating sepsis cohort files...")
    
    # Save sepsis patient list
    sepsis_patients = patient_stats_df[patient_stats_df['has_sepsis'] == True]
    if len(sepsis_patients) > 0:
        sepsis_file = OUTPUT_DIR / 'sepsis_patients.csv'
        sepsis_patients[['subject_id', 'sepsis_type', 'num_sepsis_admissions', 
                         'died']].to_csv(sepsis_file, index=False)
        print(f"✓ Saved sepsis cohort ({len(sepsis_patients):,} patients) to {sepsis_file}")
    
    # Save non-sepsis patient list
    non_sepsis_patients = patient_stats_df[patient_stats_df['has_sepsis'] == False]
    if len(non_sepsis_patients) > 0:
        non_sepsis_file = OUTPUT_DIR / 'non_sepsis_patients.csv'
        non_sepsis_patients[['subject_id', 'died']].to_csv(non_sepsis_file, index=False)
        print(f"✓ Saved non-sepsis cohort ({len(non_sepsis_patients):,} patients) to {non_sepsis_file}")
    
    # ========================================================================
    # STEP 5: PRINT SUMMARY REPORT
    # ========================================================================
    
    print("\n" + "=" * 80)
    print("ANALYSIS SUMMARY")
    print("=" * 80)
    
    print("\n📊 COHORT OVERVIEW")
    print(f"  Total patients:          {summary['cohort_overview']['total_patients']:,}")
    print(f"  Sepsis patients:         {summary['cohort_overview']['sepsis_patients']:,} "
          f"({summary['cohort_overview']['sepsis_prevalence']:.1f}%)")
    print(f"  Non-sepsis patients:     {summary['cohort_overview']['non_sepsis_patients']:,}")
    print(f"  Mortality rate:          {summary['demographics']['mortality_rate']:.1f}%")
    
    print("\n👥 DEMOGRAPHICS")
    age_mean = summary['demographics']['age_mean']
    age_std = summary['demographics']['age_std']
    if age_mean is not None and age_std is not None:
        print(f"  Mean age:                {age_mean:.1f} ± {age_std:.1f} years")
    print(f"  Gender: {summary['demographics']['gender_distribution']}")
    
    print("\n🏥 HOSPITALIZATIONS")
    print(f"  Admissions per patient:  {summary['hospitalizations']['admissions_per_patient_mean']:.1f} "
          f"(median: {summary['hospitalizations']['admissions_per_patient_median']:.0f})")
    print(f"  ICU stays per patient:   {summary['hospitalizations']['icu_stays_per_patient_mean']:.1f} "
          f"(median: {summary['hospitalizations']['icu_stays_per_patient_median']:.0f})")
    print(f"  Mean hospital LOS:       {summary['hospitalizations']['mean_los_days']:.1f} days")
    print(f"  Mean ICU LOS:            {summary['hospitalizations']['mean_icu_los_hours']:.1f} hours")
    
    print("\n🔬 CLINICAL COMPLEXITY")
    print(f"  Diagnoses per patient:   {summary['clinical_complexity']['diagnoses_per_patient_mean']:.1f} "
          f"(median: {summary['clinical_complexity']['diagnoses_per_patient_median']:.0f})")
    print(f"  Procedures per patient:  {summary['clinical_complexity']['procedures_per_patient_mean']:.1f} "
          f"(median: {summary['clinical_complexity']['procedures_per_patient_median']:.0f})")
    print(f"  Medications per patient: {summary['clinical_complexity']['medications_per_patient_mean']:.1f} "
          f"(median: {summary['clinical_complexity']['medications_per_patient_median']:.0f})")
    
    print("\n📈 DATA CHARACTERISTICS")
    print(f"  Events per patient:      {summary['data_characteristics']['events_per_patient_mean']:,.0f} "
          f"(median: {summary['data_characteristics']['events_per_patient_median']:,.0f})")
    print(f"  Lab events per patient:  {summary['data_characteristics']['lab_events_per_patient']:,.0f}")
    print(f"  Chart events per patient:{summary['data_characteristics']['chart_events_per_patient']:,.0f}")
    print(f"  Timeline per patient:    {summary['data_characteristics']['timeline_days_mean']:.1f} days")
    print(f"  Data density:            {summary['data_characteristics']['events_per_day_mean']:.1f} events/day")
    
    print("\n🦠 SEPSIS VS NON-SEPSIS COMPARISON")
    comp = summary['sepsis_comparison']
    if 'anchor_age' in comp and comp['anchor_age']['sepsis_mean'] is not None:
        print(f"  Age:                     Sepsis {comp['anchor_age']['sepsis_mean']:.1f} vs "
              f"Non-sepsis {comp['anchor_age']['non_sepsis_mean']:.1f}")
    if 'died' in comp and comp['died']['sepsis_mean'] is not None:
        print(f"  Mortality:               Sepsis {comp['died']['sepsis_mean']*100:.1f}% vs "
              f"Non-sepsis {comp['died']['non_sepsis_mean']*100:.1f}%")
    if 'num_icu_stays' in comp and comp['num_icu_stays']['ratio'] is not None:
        print(f"  ICU stays:               {comp['num_icu_stays']['ratio']:.2f}x higher in sepsis")
    if 'mean_los_days' in comp and comp['mean_los_days']['ratio'] is not None:
        print(f"  Hospital LOS:            {comp['mean_los_days']['ratio']:.2f}x longer in sepsis")
    if 'num_procedures' in comp and comp['num_procedures']['ratio'] is not None:
        print(f"  Procedures:              {comp['num_procedures']['ratio']:.2f}x more in sepsis")
    if 'num_unique_medications' in comp and comp['num_unique_medications']['ratio'] is not None:
        print(f"  Medications:             {comp['num_unique_medications']['ratio']:.2f}x more in sepsis")
    
    print("\n" + "=" * 80)
    elapsed = time.time() - overall_start
    print(f"Analysis completed in {elapsed/60:.1f} minutes")
    print(f"Results saved to: {OUTPUT_DIR}")
    print("=" * 80)

if __name__ == '__main__':
    main_sepsis()