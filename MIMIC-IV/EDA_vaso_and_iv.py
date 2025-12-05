"""
MIMIC-IV Vasopressor and IV Fluid Analysis
===========================================
Comprehensive analysis of vasopressor and IV fluid administration
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

from clinical_constants import VASOPRESSORS, IV_FLUIDS
from utils import identify_sepsis_bool

# ============================================================================
# CONFIGURATION
# ============================================================================
from config import *

OUTPUT_DIR = EDA_DIR / f"medication_analysis_v{CLEANED_VERSION}"
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
        logging.FileHandler(OUTPUT_DIR / f'medication_analysis_{datetime.now().strftime("%Y%m%d_%H%M%S")}.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# Plotting style
sns.set_style("whitegrid")
plt.rcParams['figure.figsize'] = (14, 8)
plt.rcParams['font.size'] = 10

def classify_vasopressor(drug_name: str) -> tuple:
    """
    Classify a drug as a vasopressor.
    Returns: (is_vasopressor: bool, vasopressor_type: str or None)
    """
    if not drug_name or not isinstance(drug_name, str):
        return False, None
    
    drug_lower = drug_name.lower().strip()
    
    for vaso_type, keywords in VASOPRESSORS.items():
        for keyword in keywords:
            if keyword in drug_lower:
                return True, vaso_type
    
    return False, None

def classify_iv_fluid(drug_name: str) -> tuple:
    """
    Classify a drug as an IV fluid.
    Returns: (is_iv_fluid: bool, fluid_category: str or None, fluid_type: str or None)
    """
    if not drug_name or not isinstance(drug_name, str):
        return False, None, None
    
    drug_lower = drug_name.lower().strip()
    # print(f'{drug_lower=}')
    for category, fluid_types in IV_FLUIDS.items():
        for fluid_type, keywords in fluid_types.items():
            for keyword in keywords:
                if keyword in drug_lower:
                    return True, category, fluid_type
    
    return False, None, None

def extract_patient_medications(patient_file: Path) -> dict:
    """
    Extract medication data from a patient file.
    Focus on vasopressors and IV fluids during ICU stays.
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
        
        # Extract medication events
        vasopressors = []
        iv_fluids = []
        
        for record in records:
            if not isinstance(record, dict):
                continue
            
            # Only look at prescription and input events
            event_type = record.get('event_type')
            if event_type not in ['prescription', 'input']:
                continue
            
            # Only ICU medications
#            stay_id = record.get('stay_id')
#            if stay_id not in icu_stays:
#                continue
            
            drug_name = record.get('drug') or record.get('ordercategorydescription', '')
            
            # Classify medication
            # print(f"{drug_name=}")
            # print(f"{record.get('drug')=}")
            # print(f"{record.get('ordercategorydescription')=}")
            is_vaso, vaso_type = classify_vasopressor(drug_name)
            is_fluid, fluid_cat, fluid_type = classify_iv_fluid(drug_name)
            
            if is_vaso:
                vasopressors.append({
                    'subject_id': subject_id,
                    'stay_id': stay_id,
                    'drug': drug_name,
                    'vasopressor_type': vaso_type,
                    'starttime': record.get('starttime'),
                    'stoptime': record.get('stoptime'),
                    'rate': record.get('rate'),
                    'rate_uom': record.get('rate_uom'),
                    'amount': record.get('amount'),
                    'amount_uom': record.get('amountuom'),
                    'event_type': event_type,
                })
            
            if is_fluid:
                iv_fluids.append({
                    'subject_id': subject_id,
                    'stay_id': stay_id,
                    'drug': drug_name,
                    'fluid_category': fluid_cat,
                    'fluid_type': fluid_type,
                    'starttime': record.get('starttime'),
                    'stoptime': record.get('stoptime'),
                    'rate': record.get('rate'),
                    'rate_uom': record.get('rate_uom'),
                    'amount': record.get('amount'),
                    'amount_uom': record.get('amountuom'),
                    'event_type': event_type,
                })
        
        return {
            'subject_id': subject_id,
            'has_sepsis': has_sepsis,
            'num_icu_stays': len(icu_stays),
            'vasopressors': vasopressors,
            'iv_fluids': iv_fluids,
        }
        
    except Exception as e:
        logger.error(f"Error extracting medications from {patient_file.name}: {e}")
        return {
            'subject_id': None,
            'error': str(e)
        }

def extract_medications_wrapper(patient_file: Path) -> dict:
    """Wrapper for parallel processing."""
    return extract_patient_medications(patient_file)

def extract_dose_from_rate(rate_str: str, rate_uom: str = None) -> dict:
    """
    Extract dose information from rate string.
    Returns dict with dose value and unit.
    """
    result = {'value': None, 'unit': rate_uom}
    
    if not rate_str or not isinstance(rate_str, str):
        return result
    
    try:
        # Try to extract numeric value
        import re
        numbers = re.findall(r'[-+]?\d*\.?\d+', rate_str)
        if numbers:
            result['value'] = float(numbers[0])
    except:
        pass
    
    return result

# ============================================================================
# ANALYSIS FUNCTIONS
# ============================================================================

class MedicationAnalyzer:
    """Comprehensive medication analysis."""
    
    def __init__(self, medication_data: list):
        """Initialize with list of patient medication dictionaries."""
        self.data = medication_data
        
        # Create DataFrames
        all_vasopressors = []
        all_iv_fluids = []
        
        for patient in medication_data:
            if patient.get('subject_id') is None:
                continue
            
            has_sepsis = patient.get('has_sepsis', False)
            subject_id = patient['subject_id']
            
            for vaso in patient.get('vasopressors', []):
                vaso['has_sepsis'] = has_sepsis
                all_vasopressors.append(vaso)
            
            for fluid in patient.get('iv_fluids', []):
                fluid['has_sepsis'] = has_sepsis
                all_iv_fluids.append(fluid)
        
        self.vasopressor_df = pd.DataFrame(all_vasopressors)
        self.iv_fluid_df = pd.DataFrame(all_iv_fluids)
        
        # Patient-level summary
        self.patient_summary = self._create_patient_summary()
        
        logger.info(f"Loaded {len(self.vasopressor_df):,} vasopressor records")
        logger.info(f"Loaded {len(self.iv_fluid_df):,} IV fluid records")
        logger.info(f"Patient summary: {len(self.patient_summary):,} patients")
    
    def _create_patient_summary(self) -> pd.DataFrame:
        """Create patient-level medication summary."""
        patient_data = []
        
        for patient in self.data:
            if patient.get('subject_id') is None:
                continue
            
            subject_id = patient['subject_id']
            has_sepsis = patient.get('has_sepsis', False)
            
            # Count vasopressors
            vasopressors = patient.get('vasopressors', [])
            vaso_types = Counter(v['vasopressor_type'] for v in vasopressors)
            
            # Count IV fluids
            iv_fluids = patient.get('iv_fluids', [])
            
            # Calculate total fluid amount
            total_fluid_amount = 0
            fluid_types = Counter()
            for fluid in iv_fluids:
                if fluid.get('amount'):
                    try:
                        total_fluid_amount += float(fluid['amount'])
                    except:
                        pass
                if fluid.get('fluid_type'):
                    fluid_types[fluid['fluid_type']] += 1
            
            patient_data.append({
                'subject_id': subject_id,
                'has_sepsis': has_sepsis,
                'num_icu_stays': patient.get('num_icu_stays', 0),
                'received_vasopressors': len(vasopressors) > 0,
                'num_vasopressor_administrations': len(vasopressors),
                'num_unique_vasopressors': len(vaso_types),
                'vasopressor_types': list(vaso_types.keys()),
                'received_iv_fluids': len(iv_fluids) > 0,
                'num_fluid_administrations': len(iv_fluids),
                'total_fluid_amount': total_fluid_amount if total_fluid_amount > 0 else None,
                'fluid_types': list(fluid_types.keys()),
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
            
            'vasopressor_usage': {
                'overall': {
                    'patients_received': self.patient_summary['received_vasopressors'].sum(),
                    'percentage': self.patient_summary['received_vasopressors'].mean() * 100,
                    'total_administrations': len(self.vasopressor_df),
                    'avg_per_patient': self.patient_summary['num_vasopressor_administrations'].mean(),
                },
                'sepsis': {
                    'patients_received': sepsis_patients['received_vasopressors'].sum(),
                    'percentage': sepsis_patients['received_vasopressors'].mean() * 100 if len(sepsis_patients) > 0 else 0,
                    'avg_administrations': sepsis_patients['num_vasopressor_administrations'].mean() if len(sepsis_patients) > 0 else 0,
                },
                'non_sepsis': {
                    'patients_received': non_sepsis_patients['received_vasopressors'].sum(),
                    'percentage': non_sepsis_patients['received_vasopressors'].mean() * 100 if len(non_sepsis_patients) > 0 else 0,
                    'avg_administrations': non_sepsis_patients['num_vasopressor_administrations'].mean() if len(non_sepsis_patients) > 0 else 0,
                },
            },
            
            'vasopressor_types': self._analyze_vasopressor_types(),
            
            'iv_fluid_usage': {
                'overall': {
                    'patients_received': self.patient_summary['received_iv_fluids'].sum(),
                    'percentage': self.patient_summary['received_iv_fluids'].mean() * 100,
                    'total_administrations': len(self.iv_fluid_df),
                    'avg_per_patient': self.patient_summary['num_fluid_administrations'].mean(),
                },
                'sepsis': {
                    'patients_received': sepsis_patients['received_iv_fluids'].sum(),
                    'percentage': sepsis_patients['received_iv_fluids'].mean() * 100 if len(sepsis_patients) > 0 else 0,
                    'avg_administrations': sepsis_patients['num_fluid_administrations'].mean() if len(sepsis_patients) > 0 else 0,
                    'avg_total_amount': sepsis_patients['total_fluid_amount'].mean() if len(sepsis_patients) > 0 else 0,
                },
                'non_sepsis': {
                    'patients_received': non_sepsis_patients['received_iv_fluids'].sum(),
                    'percentage': non_sepsis_patients['received_iv_fluids'].mean() * 100 if len(non_sepsis_patients) > 0 else 0,
                    'avg_administrations': non_sepsis_patients['num_fluid_administrations'].mean() if len(non_sepsis_patients) > 0 else 0,
                    'avg_total_amount': non_sepsis_patients['total_fluid_amount'].mean() if len(non_sepsis_patients) > 0 else 0,
                },
            },
            
            'iv_fluid_types': self._analyze_fluid_types(),
            
            'combined_therapy': self._analyze_combined_therapy(),
        }
        
        return summary
    
    def _analyze_vasopressor_types(self) -> dict:
        """Analyze distribution of vasopressor types."""
        result = {}
        
        if len(self.vasopressor_df) == 0:
            return result
        
        # Overall distribution
        overall_counts = self.vasopressor_df['vasopressor_type'].value_counts().to_dict()
        result['overall_distribution'] = overall_counts
        
        # By sepsis status
        sepsis_vaso = self.vasopressor_df[self.vasopressor_df['has_sepsis'] == True]
        non_sepsis_vaso = self.vasopressor_df[self.vasopressor_df['has_sepsis'] == False]
        
        result['sepsis_distribution'] = sepsis_vaso['vasopressor_type'].value_counts().to_dict()
        result['non_sepsis_distribution'] = non_sepsis_vaso['vasopressor_type'].value_counts().to_dict()
        
        return result
    
    def _analyze_fluid_types(self) -> dict:
        """Analyze distribution of IV fluid types."""
        result = {}
        
        if len(self.iv_fluid_df) == 0:
            return result
        
        # Overall distribution
        overall_counts = self.iv_fluid_df['fluid_type'].value_counts().to_dict()
        result['overall_distribution'] = overall_counts
        
        # By sepsis status
        sepsis_fluid = self.iv_fluid_df[self.iv_fluid_df['has_sepsis'] == True]
        non_sepsis_fluid = self.iv_fluid_df[self.iv_fluid_df['has_sepsis'] == False]
        
        result['sepsis_distribution'] = sepsis_fluid['fluid_type'].value_counts().to_dict()
        result['non_sepsis_distribution'] = non_sepsis_fluid['fluid_type'].value_counts().to_dict()
        
        # By category
        result['by_category'] = {
            'sepsis': sepsis_fluid['fluid_category'].value_counts().to_dict(),
            'non_sepsis': non_sepsis_fluid['fluid_category'].value_counts().to_dict(),
        }
        
        return result
    
    def _analyze_combined_therapy(self) -> dict:
        """Analyze patients receiving both vasopressors and IV fluids."""
        both_meds = self.patient_summary[
            (self.patient_summary['received_vasopressors']) & 
            (self.patient_summary['received_iv_fluids'])
        ]
        
        sepsis_both = both_meds[both_meds['has_sepsis'] == True]
        non_sepsis_both = both_meds[both_meds['has_sepsis'] == False]
        
        return {
            'overall': {
                'count': len(both_meds),
                'percentage': len(both_meds) / len(self.patient_summary) * 100 if len(self.patient_summary) > 0 else 0,
            },
            'sepsis': {
                'count': len(sepsis_both),
                'percentage': len(sepsis_both) / len(self.patient_summary[self.patient_summary['has_sepsis'] == True]) * 100 
                             if len(self.patient_summary[self.patient_summary['has_sepsis'] == True]) > 0 else 0,
            },
            'non_sepsis': {
                'count': len(non_sepsis_both),
                'percentage': len(non_sepsis_both) / len(self.patient_summary[self.patient_summary['has_sepsis'] == False]) * 100 
                             if len(self.patient_summary[self.patient_summary['has_sepsis'] == False]) > 0 else 0,
            },
        }
    
    def plot_vasopressor_distribution(self, output_dir: Path):
        """Plot vasopressor usage distribution."""
        if len(self.vasopressor_df) == 0:
            logger.warning("No vasopressor data to plot")
            return
        
        fig, axes = plt.subplots(2, 2, figsize=(16, 12))
        
        # 1. Overall vasopressor type distribution
        vaso_counts = self.vasopressor_df['vasopressor_type'].value_counts()
        axes[0, 0].barh(vaso_counts.index, vaso_counts.values, color='steelblue')
        axes[0, 0].set_xlabel('Number of Administrations')
        axes[0, 0].set_title('Vasopressor Type Distribution (All Patients)')
        axes[0, 0].grid(alpha=0.3, axis='x')
        
        # 2. Vasopressor usage: Sepsis vs Non-Sepsis
        sepsis_vaso = self.vasopressor_df[self.vasopressor_df['has_sepsis'] == True]
        non_sepsis_vaso = self.vasopressor_df[self.vasopressor_df['has_sepsis'] == False]
        
        sepsis_counts = sepsis_vaso['vasopressor_type'].value_counts()
        non_sepsis_counts = non_sepsis_vaso['vasopressor_type'].value_counts()
        
        all_types = sorted(set(list(sepsis_counts.index) + list(non_sepsis_counts.index)))
        x = np.arange(len(all_types))
        width = 0.35
        
        sepsis_vals = [sepsis_counts.get(t, 0) for t in all_types]
        non_sepsis_vals = [non_sepsis_counts.get(t, 0) for t in all_types]
        
        axes[0, 1].bar(x - width/2, sepsis_vals, width, label='Sepsis', color='coral')
        axes[0, 1].bar(x + width/2, non_sepsis_vals, width, label='Non-Sepsis', color='lightblue')
        axes[0, 1].set_xlabel('Vasopressor Type')
        axes[0, 1].set_ylabel('Number of Administrations')
        axes[0, 1].set_title('Vasopressor Usage by Sepsis Status')
        axes[0, 1].set_xticks(x)
        axes[0, 1].set_xticklabels(all_types, rotation=45, ha='right')
        axes[0, 1].legend()
        axes[0, 1].grid(alpha=0.3, axis='y')
        
        # 3. Patient-level: Percentage receiving vasopressors
        sepsis_pts = self.patient_summary[self.patient_summary['has_sepsis'] == True]
        non_sepsis_pts = self.patient_summary[self.patient_summary['has_sepsis'] == False]
        
        pct_data = pd.DataFrame({
            'Received': [
                sepsis_pts['received_vasopressors'].sum(),
                non_sepsis_pts['received_vasopressors'].sum()
            ],
            'Not Received': [
                len(sepsis_pts) - sepsis_pts['received_vasopressors'].sum(),
                len(non_sepsis_pts) - non_sepsis_pts['received_vasopressors'].sum()
            ]
        }, index=['Sepsis', 'Non-Sepsis'])
        
        pct_data.plot(kind='bar', stacked=True, ax=axes[1, 0], color=['coral', 'lightgray'])
        axes[1, 0].set_ylabel('Number of Patients')
        axes[1, 0].set_title('Vasopressor Receipt by Cohort')
        axes[1, 0].set_xticklabels(['Sepsis', 'Non-Sepsis'], rotation=0)
        axes[1, 0].legend(title='Vasopressor')
        axes[1, 0].grid(alpha=0.3, axis='y')
        
        # 4. Number of vasopressor administrations per patient
        sepsis_num = sepsis_pts[sepsis_pts['received_vasopressors']]['num_vasopressor_administrations']
        non_sepsis_num = non_sepsis_pts[non_sepsis_pts['received_vasopressors']]['num_vasopressor_administrations']
        
        axes[1, 1].hist([non_sepsis_num, sepsis_num], bins=30, label=['Non-Sepsis', 'Sepsis'], 
                       color=['lightblue', 'coral'], alpha=0.7, edgecolor='black')
        axes[1, 1].set_xlabel('Number of Vasopressor Administrations')
        axes[1, 1].set_ylabel('Number of Patients')
        axes[1, 1].set_title('Vasopressor Administrations per Patient')
        axes[1, 1].legend()
        axes[1, 1].grid(alpha=0.3, axis='y')
        axes[1, 1].set_xlim(0, 50)
        
        plt.tight_layout()
        plt.savefig(output_dir / 'vasopressor_distribution.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        logger.info("Saved vasopressor distribution plot")
    
    def plot_iv_fluid_distribution(self, output_dir: Path):
        """Plot IV fluid usage distribution."""
        if len(self.iv_fluid_df) == 0:
            logger.warning("No IV fluid data to plot")
            return

        fig, axes = plt.subplots(2, 2, figsize=(16, 12))

        # 1. Overall fluid type distribution
        fluid_counts = self.iv_fluid_df['fluid_type'].value_counts().head(10)
        axes[0, 0].barh(fluid_counts.index, fluid_counts.values, color='steelblue')
        axes[0, 0].set_xlabel('Number of Administrations')
        axes[0, 0].set_title('Top 10 IV Fluid Types (All Patients)')
        axes[0, 0].grid(alpha=0.3, axis='x')

        # 2. Fluid category: Sepsis vs Non-Sepsis
        sepsis_fluid = self.iv_fluid_df[self.iv_fluid_df['has_sepsis'] == True]
        non_sepsis_fluid = self.iv_fluid_df[self.iv_fluid_df['has_sepsis'] == False]

        sepsis_cat = sepsis_fluid['fluid_category'].value_counts()
        non_sepsis_cat = non_sepsis_fluid['fluid_category'].value_counts()

        cat_data = pd.DataFrame({
            'Sepsis': sepsis_cat,
            'Non-Sepsis': non_sepsis_cat
        }).fillna(0)

        cat_data.plot(kind='bar', ax=axes[0, 1], color=['coral', 'lightblue'])
        axes[0, 1].set_xlabel('Fluid Category')
        axes[0, 1].set_ylabel('Number of Administrations')
        axes[0, 1].set_title('IV Fluid Category by Sepsis Status')
        axes[0, 1].set_xticklabels(cat_data.index, rotation=0)
        axes[0, 1].legend()
        axes[0, 1].grid(alpha=0.3, axis='y')

        # 3. Patient-level: Percentage receiving IV fluids
        sepsis_pts = self.patient_summary[self.patient_summary['has_sepsis'] == True]
        non_sepsis_pts = self.patient_summary[self.patient_summary['has_sepsis'] == False]

        pct_data = pd.DataFrame({
            'Received': [
                int(sepsis_pts['received_iv_fluids'].sum()) if len(sepsis_pts) > 0 else 0,
                int(non_sepsis_pts['received_iv_fluids'].sum()) if len(non_sepsis_pts) > 0 else 0
            ],
            'Not Received': [
                int(len(sepsis_pts) - sepsis_pts['received_iv_fluids'].sum()) if len(sepsis_pts) > 0 else 0,
                int(len(non_sepsis_pts) - non_sepsis_pts['received_iv_fluids'].sum()) if len(non_sepsis_pts) > 0 else 0
            ]
        }, index=['Sepsis', 'Non-Sepsis'])

        pct_data.plot(kind='bar', stacked=True, ax=axes[1, 0], color=['coral', 'lightgray'])
        axes[1, 0].set_ylabel('Number of Patients')
        axes[1, 0].set_title('IV Fluid Receipt by Cohort')
        axes[1, 0].set_xticklabels(['Sepsis', 'Non-Sepsis'], rotation=0)
        axes[1, 0].legend(title='IV Fluid')
        axes[1, 0].grid(alpha=0.3, axis='y')

        # 4. Total fluid amount per patient (boxplot)
        # Use total_fluid_amount from patient_summary (may have NaNs)
        combined = self.patient_summary[['has_sepsis', 'total_fluid_amount']].dropna()
        if not combined.empty:
            sns.boxplot(x='has_sepsis', y='total_fluid_amount', data=combined.replace({True: 'Sepsis', False: 'Non-Sepsis'}),
                        ax=axes[1, 1])
            axes[1, 1].set_xlabel('Cohort')
            axes[1, 1].set_ylabel('Total Fluid Amount')
            axes[1, 1].set_title('Total Fluid Amount per Patient by Cohort')
            axes[1, 1].set_xticklabels(['Non-Sepsis', 'Sepsis'])
            axes[1, 1].grid(alpha=0.3, axis='y')
        else:
            axes[1, 1].text(0.5, 0.5, 'No total fluid amount data', ha='center', va='center')
            axes[1, 1].set_axis_off()

        plt.tight_layout()
        plt.savefig(output_dir / 'iv_fluid_distribution.png', dpi=300, bbox_inches='tight')
        plt.close()

        logger.info("Saved IV fluid distribution plot")

    def plot_combined_therapy(self, output_dir: Path):
        """Plot combined therapy (vasopressors + fluids) summary."""
        both_meds = self.patient_summary[
            (self.patient_summary['received_vasopressors']) &
            (self.patient_summary['received_iv_fluids'])
        ]

        if len(self.patient_summary) == 0:
            logger.warning("No patient summary data to plot combined therapy")
            return

        fig, ax = plt.subplots(1, 1, figsize=(8, 6))
        counts = [
            len(both_meds),
            len(self.patient_summary[self.patient_summary['received_vasopressors']]) - len(both_meds),
            len(self.patient_summary[self.patient_summary['received_iv_fluids']]) - len(both_meds),
            len(self.patient_summary) - (len(both_meds) + (len(self.patient_summary[self.patient_summary['received_vasopressors']]) - len(both_meds)) + (len(self.patient_summary[self.patient_summary['received_iv_fluids']]) - len(both_meds)))
        ]
        labels = ['Both', 'Vaso only', 'Fluids only', 'Neither']
        ax.pie(counts, labels=labels, autopct='%1.1f%%', startangle=90)
        ax.set_title('Therapy Overlap: Vasopressors and IV Fluids')
        plt.savefig(output_dir / 'combined_therapy_pie.png', dpi=300, bbox_inches='tight')
        plt.close()
        logger.info("Saved combined therapy plot")


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
    if isinstance(obj, (np.integer, np.int64)):
        return int(obj)
    if isinstance(obj, (np.floating, np.float64)):
        return float(obj)
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    return obj

def run_parallel_extraction(patient_files: list, num_workers: int = 4) -> list:
    """
    Run medication extraction in parallel using multiprocessing.
    """
    logger.info(f"Starting parallel extraction with {num_workers} workers...")
    
    start = time.time()
    with Pool(num_workers) as pool:
        results = list(tqdm(pool.imap(extract_medications_wrapper, patient_files),
                            total=len(patient_files)))
    end = time.time()
    
    logger.info(f"Completed extraction in {end - start:.2f} seconds")
    return results


def save_raw_extraction_results(results: list, cache_dir: Path):
    """
    Save raw extraction results to disk for reuse.
    """
    cache_file = cache_dir / "medication_extraction_results.json"
    
    with open(cache_file, "w") as f:
        json.dump(results, f)
    
    logger.info(f"Saved medication extraction results to {cache_file}")


def load_cached_results(cache_dir: Path) -> list | None:
    """
    Load previously saved extraction results if available.
    """
    cache_file = cache_dir / "medication_extraction_results.json"
    
    if cache_file.exists():
        logger.info(f"Loading cached extraction results from {cache_file}")
        with open(cache_file, "r") as f:
            return json.load(f)
    else:
        return None


def main_vaso_and_iv():
    logger.info("====================================================")
    logger.info("   MIMIC-IV Medication Analysis: Vasopressor & IV  ")
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
        medication_data = cached
    else:
        medication_data = run_parallel_extraction(patient_files, NUM_WORKERS)
        save_raw_extraction_results(medication_data, CACHE_DIR)
    
    # ------------------------------------------------------------
    # ANALYZE
    # ------------------------------------------------------------
    logger.info("Initializing MedicationAnalyzer...")
    analyzer = MedicationAnalyzer(medication_data)
    
    logger.info("Generating summary statistics...")
    summary_stats = analyzer.generate_summary_statistics()
    
    # Save summary statistics
    summary_file = OUTPUT_DIR / "summary_statistics.json"
    with open(summary_file, "w") as f:
        json.dump(summary_stats, f, indent=4, default=to_serializable)
#        json.dump(summary_stats, f, indent=4)
    
    logger.info(f"Saved summary statistics to {summary_file}")
    
    # ------------------------------------------------------------
    # GENERATE PLOTS
    # ------------------------------------------------------------
    logger.info("Generating plots...")
    
    try:
        analyzer.plot_vasopressor_distribution(FIGURES_DIR)
    except Exception as e:
        logger.error(f"Error plotting vasopressor distribution: {e}")
    
    try:
        analyzer.plot_iv_fluid_distribution(FIGURES_DIR)
    except Exception as e:
        logger.error(f"Error plotting IV fluid distribution: {e}")
    
    logger.info("All tasks completed successfully.")

if __name__ == "__main__":
    main_vaso_and_iv()
