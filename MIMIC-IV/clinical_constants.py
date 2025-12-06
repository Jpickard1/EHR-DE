# ============================================================================
# SEPSIS IDENTIFICATION
# ============================================================================

# Sepsis ICD-10 codes (comprehensive list)
SEPSIS_ICD10_CODES = {
    # Sepsis
    'A40', 'A400', 'A401', 'A403', 'A408', 'A409',  # Streptococcal sepsis
    'A41', 'A410', 'A411', 'A412', 'A413', 'A414', 'A415', 'A4150', 'A4151', 'A4152', 'A4153', 'A4159',
    'A418', 'A4181', 'A4189', 'A419',  # Other sepsis
    
    # Severe sepsis
    'R6520', 'R6521',  # Severe sepsis without/with septic shock
    
    # Septic shock
    'R6521', 'T8111', 'T8112',  # Septic shock, postprocedural septic shock
    
    # Bacteremia
    'A499', 'R780', 'R7881',  # Bacterial infection/bacteremia
}

# Sepsis ICD-9 codes (for legacy data)
SEPSIS_ICD9_CODES = {
    '038', '0380', '0381', '03810', '03811', '03812', '03819',
    '0382', '0383', '0384', '03840', '03841', '03842', '03843', '03844', '03849',
    '0388', '0389',  # Septicemia
    '99591', '99592',  # Sepsis, severe sepsis
    '78552',  # Septic shock
}

# Organ dysfunction indicators (for Sepsis-3 criteria)
ORGAN_DYSFUNCTION_ICD10 = {
    'respiratory': {'J80', 'J96', 'J9600', 'J9601', 'J9602', 'J9690', 'J9691', 'J9692'},
    'cardiovascular': {'I95', 'I959', 'R6521'},  # Shock
    'renal': {'N17', 'N170', 'N171', 'N172', 'N178', 'N179'},  # AKI
    'hepatic': {'K72', 'K7200', 'K7201', 'K7210', 'K7211', 'K7290', 'K7291'},
    'coagulation': {'D65', 'D688', 'D689'},  # DIC
    'neurological': {'G93', 'G931', 'R402', 'R404'},
}


# ============================================================================
# MEDICATION CLASSIFICATION
# ============================================================================

# Vasopressor medications (commonly used in ICU)
VASOPRESSORS = {
    'norepinephrine': ['norepinephrine', 'levophed', 'noradrenaline'],
    'epinephrine': ['epinephrine', 'adrenaline'],
    'dopamine': ['dopamine'],
    'vasopressin': ['vasopressin', 'pitressin'],
    'phenylephrine': ['phenylephrine', 'neosynephrine'],
    'dobutamine': ['dobutamine'],
    'milrinone': ['milrinone'],
    'isoproterenol': ['isoproterenol', 'isuprel'],
}

# IV Fluid types
IV_FLUIDS = {
    'crystalloid': {
        'normal_saline': ['normal saline', '0.9% sodium chloride', '0.9% nacl', 'ns', 'sodium chloride 0.9%', '0.9% Sodium Chloride'],
        'lactated_ringers': ['lactated ringer', "ringer's lactate", 'lr', 'rl'],
        'dextrose': ['dextrose', 'd5w', 'd5', 'dextrose 5%'],
        'half_normal_saline': ['0.45% sodium chloride', '0.45% nacl', 'half normal saline'],
    },
    'colloid': {
        'albumin': ['albumin'],
        'hetastarch': ['hetastarch', 'hespan'],
        'dextran': ['dextran'],
    }
}

# Steroids
STEROIDS = {
    'hydrocortisone': ['hydrocortisone', 'cortef', 'solu-cortef'],
    'methylprednisolone': ['methylprednisolone', 'solu-medrol', 'medrol'],
    'dexamethasone': ['dexamethasone', 'decadron'],
    'prednisone': ['prednisone', 'deltasone'],
    'prednisolone': ['prednisolone'],
    'fludrocortisone': ['fludrocortisone', 'florinef'],
}


"""
Minimal ICD Code Verification for MIMIC-IV
Verifies that ICD codes match their expected descriptions
"""

import pandas as pd
from pathlib import Path
from typing import Dict, Set, List



#!/usr/bin/env python3
"""
Minimal ICD Code Verification for MIMIC-IV
Verifies that ICD codes match their expected descriptions
"""

import pandas as pd
from pathlib import Path
from typing import Dict, Set, List

# ============================================================================
# VERIFICATION FUNCTIONS
# ============================================================================

def load_icd_dictionaries(mimic_data_dir: str) -> Dict[str, pd.DataFrame]:
    """Load ICD diagnosis dictionaries from MIMIC-IV."""
    mimic_path = Path(mimic_data_dir)
    dicts = {}
    
    # Load ICD-10 diagnoses
    icd10_path = mimic_path / 'hosp' / 'd_icd_diagnoses.csv.gz'
    if icd10_path.exists():
        df = pd.read_csv(icd10_path, compression='gzip')
        dicts['icd10'] = df[df['icd_version'] == 10]
        dicts['icd9'] = df[df['icd_version'] == 9]
        print(f"✓ Loaded {len(dicts['icd10'])} ICD-10 codes")
        print(f"✓ Loaded {len(dicts['icd9'])} ICD-9 codes")
    else:
        print(f"✗ Dictionary not found at {icd10_path}")
        return {}
    
    return dicts


def verify_icd_codes(codes: Set[str], 
                     icd_df: pd.DataFrame, 
                     version: str,
                     category: str = "") -> pd.DataFrame:
    """
    Verify ICD codes exist and show their descriptions.
    
    Args:
        codes: Set of ICD codes to verify
        icd_df: DataFrame from d_icd_diagnoses.csv.gz
        version: 'ICD-9' or 'ICD-10'
        category: Optional category label (e.g., 'Sepsis', 'Respiratory')
    """
    results = []
    
    print(f"\n{'='*80}")
    print(f"{category.upper()} - {version}" if category else version)
    print(f"{'='*80}")
    
    for code in sorted(codes):
        # Try exact match first
        match = icd_df[icd_df['icd_code'] == code]
        
        # If no exact match, try without dots (ICD-10 can have dots)
        if len(match) == 0:
            code_nodot = code.replace('.', '')
            match = icd_df[icd_df['icd_code'] == code_nodot]
        
        if len(match) > 0:
            desc = match.iloc[0]['long_title']
            results.append({
                'code': code,
                'found': True,
                'description': desc
            })
            print(f"✓ {code:8s} {desc}")
        else:
            results.append({
                'code': code,
                'found': False,
                'description': 'NOT FOUND'
            })
            print(f"✗ {code:8s} NOT FOUND IN MIMIC DICTIONARY")
    
    df_results = pd.DataFrame(results)
    found_count = df_results['found'].sum()
    total_count = len(df_results)
    
    print(f"\nSummary: {found_count}/{total_count} codes verified")
    
    return df_results


def verify_all_codes(mimic_data_dir: str) -> Dict[str, pd.DataFrame]:
    """Main verification function."""
    
    print("="*80)
    print("ICD CODE VERIFICATION FOR SEPSIS IDENTIFICATION")
    print("="*80)
    
    # Load dictionaries
    dicts = load_icd_dictionaries(mimic_data_dir)
    if not dicts:
        print("\n⚠️  ERROR: Could not load ICD dictionaries!")
        return {}
    
    results = {}
    
    # Verify ICD-10 sepsis codes
    results['sepsis_icd10'] = verify_icd_codes(
        SEPSIS_ICD10_CODES,
        dicts['icd10'],
        'ICD-10',
        'SEPSIS CODES'
    )
    
    # Verify ICD-9 sepsis codes
    results['sepsis_icd9'] = verify_icd_codes(
        SEPSIS_ICD9_CODES,
        dicts['icd9'],
        'ICD-9',
        'SEPSIS CODES'
    )
    
    # Verify organ dysfunction codes by system
    for system, codes in ORGAN_DYSFUNCTION_ICD10.items():
        results[f'organ_{system}'] = verify_icd_codes(
            codes,
            dicts['icd10'],
            'ICD-10',
            f'ORGAN DYSFUNCTION - {system.upper()}'
        )
    
    # Overall summary
    print("\n" + "="*80)
    print("OVERALL SUMMARY")
    print("="*80)
    
    total_codes = sum(len(df) for df in results.values())
    total_found = sum(df['found'].sum() for df in results.values())
    
    print(f"Total codes checked: {total_codes}")
    print(f"Total codes found: {total_found}")
    print(f"Total codes missing: {total_codes - total_found}")
    
    if total_found == total_codes:
        print("\n✓ ALL CODES VERIFIED SUCCESSFULLY!")
    else:
        print(f"\n⚠️  WARNING: {total_codes - total_found} codes not found")
        print("Missing codes may need to be:")
        print("  1. Corrected (typo in code)")
        print("  2. Updated (different in MIMIC-IV)")
        print("  3. Removed (not in MIMIC-IV)")
    
    return results


# ============================================================================
# USAGE EXAMPLE
# ============================================================================

if __name__ == "__main__":
    # Update this path to your MIMIC-IV data directory
    # Path to your MIMIC-IV data
    from config import DATA_DIR
    MIMIC_DATA_DIR = str(DATA_DIR)
        
    results = verify_all_codes(MIMIC_DATA_DIR)
    
    # Optional: Save results to CSV
    # for name, df in results.items():
    #     df.to_csv(f'icd_verification_{name}.csv', index=False)