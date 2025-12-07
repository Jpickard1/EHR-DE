# ============================================================================
# MIMIC-IV FILE STRUCTURE EXPLANATION
# ============================================================================
"""
MIMIC-IV File Organization:

The 'd_' prefix stands for "dictionary" or "definition" tables. These are lookup/reference tables.

Directory Structure:
- hosp/: Hospital-wide data (admissions, labs, medications, diagnoses)
  - d_labitems.csv.gz: Dictionary of lab test definitions (itemid -> lab name)
  - d_icd_diagnoses.csv.gz: ICD diagnosis code definitions
  - d_icd_procedures.csv.gz: ICD procedure code definitions
  - d_hcpcs.csv.gz: HCPCS code definitions
  - labevents.csv.gz: Actual lab results (uses itemids from d_labitems)
  - prescriptions.csv.gz: Medication prescriptions
  - etc.

- icu/: ICU-specific data (vitals, ventilator settings, input/output)
  - d_items.csv.gz: Dictionary of ICU measurements (itemid -> measurement name)
  - chartevents.csv.gz: Vital signs and assessments (uses itemids from d_items)
  - inputevents.csv.gz: Fluid/medication inputs
  - outputevents.csv.gz: Output measurements
  - procedureevents.csv.gz: Procedures performed
  - etc.

Key Concept: 
- Data tables (labevents, chartevents) reference itemids
- Dictionary tables (d_labitems, d_items) define what each itemid means
- This saves space vs storing text labels for every measurement
"""

# ============================================================================
# ITEMID MAPPINGS FOR MIMIC-IV
# ============================================================================

# Vital signs and clinical measurements
# Updated for MIMIC-IV v3.1 - removed outdated itemIDs from MIMIC-III
ITEMID_MAP = {
    # Respiratory rate (ICU)
    'respiratory_rate': [
        220045,     # Heart Rate [ICU]
        220210,     # Respiratory Rate [ICU]
        224688,     # Respiratory Rate (Set) [ICU]
        224689,     # Respiratory Rate (spontaneous) [ICU]
        224690      # Respiratory Rate (Total) [ICU]
    ],
    
    # Mean arterial pressure (ICU)
    'map': [
        220052,     # Arterial Blood Pressure mean [ICU]
        220181,     # Non Invasive Blood Pressure mean [ICU]
        225312      # ART BP Mean [ICU]
    ],
    
    # Systolic blood pressure (ICU)
    'sbp': [
        220050,     # Arterial Blood Pressure systolic [ICU]
        220179      # Non Invasive Blood Pressure systolic [ICU]
    ],
    
    # Glasgow Coma Scale components (ICU)
    'gcs_eye': [220739],      # GCS - Eye Opening [ICU]
    'gcs_verbal': [223900],   # GCS - Verbal Response [ICU]
    'gcs_motor': [223901],    # GCS - Motor Response [ICU]
    
    # Oxygen measurements (ICU)
    'fio2': [223835],         # Inspired O2 Fraction [ICU]
    'pao2': [220224],         # Arterial O2 pressure [ICU]
    'spo2': [
        220277,               # O2 saturation pulseoxymetry [ICU]
        220227                # Arterial O2 Saturation [ICU]
    ],
    
    # Laboratory values (from d_labitems - these are 5-digit codes)
    'bilirubin': [50885],     # Bilirubin, Total [LAB]
    'creatinine': [50912],    # Creatinine [LAB]
    'platelets': [51265],     # Platelet Count [LAB]
    'weight': [
        226512,     # Admission Weight (Kg) [BOTH]
        226531,     # Admission Weight (lbs.) [BOTH]
        224639,     # Daily Weight [BOTH]
    ],
    
    # Urine output (ICU output events)
    'urine_output': [
        226559,     # Foley [ICU]
        226560,     # Void [ICU]
        226561,     # Condom Cath [ICU]
        226563,     # Suprapubic [ICU]
        226564,     # R Nephrostomy [ICU]
        226584      # Ileoconduit
    ],
}

# ============================================================================
# HELPER METHODS
# ============================================================================
import pandas as pd
from typing import Dict, List, Tuple, Optional
import logging

def get_weight_kg(measure_df):
    weight_kg = []
    for idx, row in measure_df.iterrows():
        itemid = row["itemid"]
        value = row["valuenum"]
        uom = str(row.get("valueuom", "")).lower().strip()

        # Admission Weight (Kg)
        if itemid == 226512:
            print("A")
            weight_kg.append(value)

        # Admission Weight (lbs)
        elif itemid == 226531:
            print("B")
            weight_kg.append(value * 0.45359237)

        # Daily Weight (varies)
        elif itemid == 224639:
            if "kg" in uom:
                print("C")
                weight_kg.append(value)
            elif "lb" in uom:
                print("D")
                weight_kg.append(value * 0.45359237)
            else:
                print("E")
                weight_kg.append(None)

        else:
            print("F")
            weight_kg.append(None)

    measure_df['weight_kg'] = weight_kg
    return weight_kg, measure_df

# ============================================================================
# ITEMID VERIFICATION FUNCTIONS
# ============================================================================

def load_mimic_dictionaries(mimic_data_dir: str) -> Dict[str, pd.DataFrame]:
    """
    Load all MIMIC-IV dictionary files.
    
    Args:
        mimic_data_dir: Path to MIMIC-IV data directory (e.g., '/path/to/mimiciv/3.1')
        
    Returns:
        Dictionary of DataFrames with all dictionary tables
    """
    from pathlib import Path
    
    mimic_path = Path(mimic_data_dir)
    dictionaries = {}
    
    # ICU dictionary (vitals, assessments, ventilator settings)
    icu_d_items_path = mimic_path / 'icu' / 'd_items.csv.gz'
    if icu_d_items_path.exists():
        dictionaries['icu_items'] = pd.read_csv(icu_d_items_path, compression='gzip')
        print(f"✓ Loaded ICU d_items: {len(dictionaries['icu_items'])} items")
    else:
        print(f"✗ ICU d_items not found at {icu_d_items_path}")
    
    # Hospital lab dictionary (lab tests)
    hosp_d_labitems_path = mimic_path / 'hosp' / 'd_labitems.csv.gz'
    if hosp_d_labitems_path.exists():
        dictionaries['lab_items'] = pd.read_csv(hosp_d_labitems_path, compression='gzip')
        print(f"✓ Loaded d_labitems: {len(dictionaries['lab_items'])} items")
    else:
        print(f"✗ d_labitems not found at {hosp_d_labitems_path}")
    
    # ICD diagnosis codes
    d_icd_diag_path = mimic_path / 'hosp' / 'd_icd_diagnoses.csv.gz'
    if d_icd_diag_path.exists():
        dictionaries['icd_diagnoses'] = pd.read_csv(d_icd_diag_path, compression='gzip')
        print(f"✓ Loaded d_icd_diagnoses: {len(dictionaries['icd_diagnoses'])} codes")
    
    # ICD procedure codes
    d_icd_proc_path = mimic_path / 'hosp' / 'd_icd_procedures.csv.gz'
    if d_icd_proc_path.exists():
        dictionaries['icd_procedures'] = pd.read_csv(d_icd_proc_path, compression='gzip')
        print(f"✓ Loaded d_icd_procedures: {len(dictionaries['icd_procedures'])} codes")
    
    if not dictionaries:
        print(f"\n⚠️  ERROR: No dictionary files found!")
        print(f"Searched in: {mimic_path}")
        print(f"Please verify the path points to the MIMIC-IV root directory")
    
    return dictionaries


def verify_itemids_from_dictionaries(mimic_data_dir: str) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Load MIMIC-IV dictionaries to verify itemIDs.
    
    Args:
        mimic_data_dir: Path to MIMIC-IV data directory
        
    Returns:
        Tuple of (icu_items_df, lab_items_df)
    """
    dictionaries = load_mimic_dictionaries(mimic_data_dir)
    
    icu_items = dictionaries.get('icu_items', pd.DataFrame())
    lab_items = dictionaries.get('lab_items', pd.DataFrame())
    
    if len(icu_items) > 0:
        print(f"\nICU d_items columns: {icu_items.columns.tolist()}")
    
    if len(lab_items) > 0:
        print(f"Lab d_labitems columns: {lab_items.columns.tolist()}")
    
    return icu_items, lab_items


def search_itemids_by_keyword(icu_items: pd.DataFrame, 
                              lab_items: pd.DataFrame, 
                              keyword: str) -> Dict[str, pd.DataFrame]:
    """
    Search for itemIDs by keyword in both ICU and lab dictionaries.
    
    Args:
        icu_items: DataFrame from d_items.csv.gz (ICU measurements)
        lab_items: DataFrame from d_labitems.csv.gz (lab tests)
        keyword: Search term (e.g., 'respiratory rate', 'creatinine')
        
    Returns:
        Dictionary with 'icu' and 'lab' DataFrames of matches
    """
    keyword_lower = keyword.lower()
    results = {}
    
    # Search ICU items
    if len(icu_items) > 0:
        icu_matches = icu_items[
            icu_items['label'].str.lower().str.contains(keyword_lower, na=False)
        ]
        if 'abbreviation' in icu_items.columns:
            icu_matches = pd.concat([
                icu_matches,
                icu_items[icu_items['abbreviation'].str.lower().str.contains(keyword_lower, na=False)]
            ]).drop_duplicates('itemid')
        
        if len(icu_matches) > 0:
            display_cols = ['itemid', 'label']
            if 'abbreviation' in icu_matches.columns:
                display_cols.append('abbreviation')
            if 'linksto' in icu_matches.columns:
                display_cols.append('linksto')
            if 'category' in icu_matches.columns:
                display_cols.append('category')
            if 'unitname' in icu_matches.columns:
                display_cols.append('unitname')
            
            results['icu'] = icu_matches[display_cols]
    
    # Search lab items
    if len(lab_items) > 0:
        lab_matches = lab_items[
            lab_items['label'].str.lower().str.contains(keyword_lower, na=False)
        ]
        if 'fluid' in lab_items.columns:
            lab_matches = pd.concat([
                lab_matches,
                lab_items[lab_items['fluid'].str.lower().str.contains(keyword_lower, na=False)]
            ]).drop_duplicates('itemid')
        
        if len(lab_matches) > 0:
            display_cols = ['itemid', 'label']
            if 'fluid' in lab_matches.columns:
                display_cols.append('fluid')
            if 'category' in lab_matches.columns:
                display_cols.append('category')
            
            results['lab'] = lab_matches[display_cols]
    
    return results


def verify_current_itemids(icu_items: pd.DataFrame, lab_items: pd.DataFrame) -> dict:
    """
    Verify all itemIDs in ITEMID_MAP exist in dictionaries and show their labels.
    
    Args:
        icu_items: DataFrame from d_items.csv.gz
        lab_items: DataFrame from d_labitems.csv.gz
        
    Returns:
        Dictionary with verification results
    """
    results = {}
    
    print("\n" + "="*80)
    print("VERIFYING ITEMID MAPPINGS")
    print("="*80)
    print("\nNote: ICU items (vitals, assessments) from d_items.csv.gz")
    print("      Lab items (lab tests) from d_labitems.csv.gz")
    
    # Determine which items should be ICU vs Lab
    icu_measures = [
        'respiratory_rate', 'map', 'sbp', 'gcs_eye', 'gcs_verbal', 
        'gcs_motor', 'gcs_total', 'fio2', 'pao2', 'spo2', 'urine_output'
    ]
    lab_measures = ['bilirubin', 'creatinine', 'platelets']
    
    for measure_name, itemids in ITEMID_MAP.items():
        print(f"\n{measure_name.upper()}:")
        print("-" * 60)
        
        found_items = []
        missing_items = []
        
        # Determine which dictionary to search
        if measure_name in icu_measures:
            search_df = icu_items
            source = "ICU"
        elif measure_name in lab_measures:
            search_df = lab_items
            source = "LAB"
        else:
            # Try both
            search_df = pd.concat([icu_items, lab_items]) if len(icu_items) > 0 and len(lab_items) > 0 else (icu_items if len(icu_items) > 0 else lab_items)
            source = "BOTH"
        
        for itemid in itemids:
            match = search_df[search_df['itemid'] == itemid]
            
            if len(match) > 0:
                item_info = match.iloc[0]
                found_items.append({
                    'itemid': itemid,
                    'label': item_info['label'],
                    'source': source,
                    'category': item_info.get('category', 'N/A')
                })
                print(f"  ✓ {itemid}: {item_info['label']} [{source}]")
            else:
                missing_items.append(itemid)
                print(f"  ✗ {itemid}: NOT FOUND in {source} dictionary")
        
        results[measure_name] = {
            'found': found_items,
            'missing': missing_items,
            'found_count': len(found_items),
            'missing_count': len(missing_items)
        }
    
    # Summary
    print("\n" + "="*80)
    print("SUMMARY")
    print("="*80)
    
    total_found = sum(r['found_count'] for r in results.values())
    total_missing = sum(r['missing_count'] for r in results.values())
    
    print(f"Total itemIDs found: {total_found}")
    print(f"Total itemIDs missing: {total_missing}")
    
    if total_missing > 0:
        print("\n⚠️  WARNING: Some itemIDs were not found. You may need to:")
        print("   1. Check if you're using the correct MIMIC-IV version (you have 3.1)")
        print("   2. Search for alternative itemIDs using search_itemids_by_keyword()")
        print("   3. Remove invalid itemIDs from ITEMID_MAP")
        print("   4. Some itemIDs may have changed between MIMIC-IV versions")
    else:
        print("\n✓ All itemIDs verified successfully!")
    
    return results


def analyze_patient_itemids(patient_files: List[str], sample_size: int = 100) -> pd.DataFrame:
    """
    Analyze actual itemIDs present in patient data files.
    
    Args:
        patient_files: List of patient JSON file paths
        sample_size: Number of patients to sample
        
    Returns:
        DataFrame with itemID frequency analysis
    """
    import json
    from collections import Counter
    
    itemid_counter = Counter()
    value_samples = {}
    
    # Sample patient files
    if len(patient_files) > sample_size:
        import random
        patient_files = random.sample(patient_files, sample_size)
    
    print(f"Analyzing itemIDs from {len(patient_files)} patient files...")
    
    for pf in patient_files:
        try:
            with open(pf, 'r') as f:
                data = json.load(f)
            
            records = data.get('records', [])
            for record in records:
                if isinstance(record, dict):
                    itemid = record.get('itemid')
                    if itemid:
                        itemid_counter[itemid] += 1
                        
                        # Store sample values
                        if itemid not in value_samples:
                            value_samples[itemid] = {
                                'value': record.get('value'),
                                'valuenum': record.get('valuenum'),
                                'valueuom': record.get('valueuom')
                            }
        except Exception as e:
            print(f"Error reading {pf}: {e}")
            continue
    
    # Create summary DataFrame
    summary = pd.DataFrame([
        {
            'itemid': itemid,
            'count': count,
            'sample_value': value_samples[itemid]['value'],
            'sample_valuenum': value_samples[itemid]['valuenum'],
            'sample_uom': value_samples[itemid]['valueuom']
        }
        for itemid, count in itemid_counter.most_common(100)
    ])
    
    return summary


def find_alternative_itemids(icu_items: pd.DataFrame, 
                            lab_items: pd.DataFrame, 
                            current_measure: str) -> Dict[str, pd.DataFrame]:
    """
    Suggest alternative itemIDs for a given measurement.
    
    Args:
        icu_items: DataFrame from d_items.csv.gz
        lab_items: DataFrame from d_labitems.csv.gz
        current_measure: Measurement name from ITEMID_MAP keys
        
    Returns:
        Dictionary with suggested alternative itemIDs from ICU and/or Lab
    """
    search_terms = {
        'respiratory_rate': ['respiratory rate', 'resp rate', 'rr'],
        'map': ['mean arterial pressure', 'arterial bp mean', 'map', 'art bp mean'],
        'sbp': ['systolic', 'arterial bp systolic', 'sbp', 'blood pressure systolic'],
        'gcs_eye': ['gcs eye', 'glasgow eye'],
        'gcs_verbal': ['gcs verbal', 'glasgow verbal'],
        'gcs_motor': ['gcs motor', 'glasgow motor'],
        'gcs_total': ['gcs total', 'glasgow total', 'glasgow coma scale'],
        'fio2': ['fio2', 'inspired o2 fraction', 'fraction inspired oxygen'],
        'pao2': ['pao2', 'po2', 'arterial o2', 'oxygen, arterial'],
        'spo2': ['spo2', 'o2 saturation', 'oxygen saturation', 'pulse ox'],
        'bilirubin': ['bilirubin', 'total bili'],
        'creatinine': ['creatinine'],
        'platelets': ['platelet', 'plt'],
        'urine_output': ['urine', 'foley', 'void']
    }
    
    if current_measure not in search_terms:
        print(f"Unknown measure: {current_measure}")
        return {}
    
    print(f"\nSearching for alternatives for: {current_measure}")
    print("="*60)
    
    all_results = {}
    
    for term in search_terms[current_measure]:
        results = search_itemids_by_keyword(icu_items, lab_items, term)
        
        if results.get('icu') is not None and len(results['icu']) > 0:
            print(f"\nICU matches for '{term}':")
            print(results['icu'].to_string(index=False))
            if 'icu' not in all_results:
                all_results['icu'] = results['icu']
            else:
                all_results['icu'] = pd.concat([all_results['icu'], results['icu']]).drop_duplicates('itemid')
        
        if results.get('lab') is not None and len(results['lab']) > 0:
            print(f"\nLab matches for '{term}':")
            print(results['lab'].to_string(index=False))
            if 'lab' not in all_results:
                all_results['lab'] = results['lab']
            else:
                all_results['lab'] = pd.concat([all_results['lab'], results['lab']]).drop_duplicates('itemid')
    
    return all_results

def example_verification():
    """Example of how to verify itemIDs."""
    
    # Path to your MIMIC-IV data
    from config import DATA_DIR
    mimic_data_dir = str(DATA_DIR)
    
    print("="*80)
    print("MIMIC-IV ITEMID VERIFICATION TOOL")
    print("="*80)
    print("\nStep 1: Loading MIMIC-IV dictionaries...")
    print("-"*80)
    
    # 1. Load the dictionaries
    icu_items, lab_items = verify_itemids_from_dictionaries(mimic_data_dir)
    
    if len(icu_items) == 0 and len(lab_items) == 0:
        print("\n⚠️  ERROR: Could not load any dictionaries!")
        return
    
    print("\nStep 2: Verifying itemIDs in ITEMID_MAP...")
    print("-"*80)
    
    # 2. Verify all current itemIDs
    verification_results = verify_current_itemids(icu_items, lab_items)
    
    print("\nStep 3: Example searches...")
    print("-"*80)
    
    # 3. Search for specific terms
    print("\nExample A: Searching for 'respiratory rate'")
    print("="*60)
    resp_results = search_itemids_by_keyword(icu_items, lab_items, 'respiratory rate')
    if resp_results.get('icu') is not None:
        print("\nICU items:")
        print(resp_results['icu'].to_string(index=False))
    
    print("\nExample B: Searching for 'creatinine'")
    print("="*60)
    creat_results = search_itemids_by_keyword(icu_items, lab_items, 'creatinine')
    if creat_results.get('lab') is not None:
        print("\nLab items:")
        print(creat_results['lab'].to_string(index=False))
    
    # 4. Find alternatives for a specific measure
    print("\nStep 4: Finding alternatives for specific measures...")
    print("-"*80)
    print("\nExample: Finding alternatives for MAP")
    print("="*60)
    alternatives = find_alternative_itemids(icu_items, lab_items, 'map')
    
    return verification_results


if __name__ == "__main__":
    # Set up logging
    logging.basicConfig(level=logging.INFO)
    
    # Run verification example
    print("To verify itemIDs, run: example_verification()")
    print("To calculate scores, run: example_usage()")
    
    # Uncomment to run:
    example_verification()