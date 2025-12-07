# ============================================================================
# SEPSIS IDENTIFICATION
# ============================================================================

# Sepsis ICD-10 codes (comprehensive list)
SEPSIS_ICD10_CODES = {
    'A40',      # Streptococcal sepsis
    'A400',     # Sepsis due to streptococcus, group A
    'A401',     # Sepsis due to streptococcus, group B
    'A403',     # Sepsis due to Streptococcus pneumoniae
    'A408',     # Other streptococcal sepsis
    'A409',     # Streptococcal sepsis, unspecified
    'A41',      # Other sepsis
    'A410',     # Sepsis due to Staphylococcus aureus
    'A411',     # Sepsis due to other specified staphylococcus
    'A412',     # Sepsis due to unspecified staphylococcus
    'A413',     # Sepsis due to Hemophilus influenzae
    'A414',     # Sepsis due to anaerobes
    'A415',     # Sepsis due to other Gram-negative organisms
    'A4150',    # Gram-negative sepsis, unspecified
    'A4151',    # Sepsis due to Escherichia coli [E. coli]
    'A4152',    # Sepsis due to Pseudomonas
    'A4153',    # Sepsis due to Serratia
    'A4159',    # Other Gram-negative sepsis
    'A418',     # Other specified sepsis
    'A4181',    # Sepsis due to Enterococcus
    'A4189',    # Other specified sepsis
    'A419',     # Sepsis, unspecified organism
    'A499',     # Bacterial infection, unspecified
    'R6520',    # Severe sepsis without septic shock
    'R6521',    # Severe sepsis with septic shock
    'R7881',    # Bacteremia
    'T8111',    # Postprocedural cardiogenic shock
    'T8112'     # Postprocedural septic shock
}

# Sepsis ICD-9 codes (for legacy data)
SEPSIS_ICD9_CODES = {
    '0380',     # Streptococcal septicemia
    '03810',    # Staphylococcal septicemia, unspecified
    '03811',    # Methicillin susceptible Staphylococcus aureus septicemia
    '03812',    # Methicillin resistant Staphylococcus aureus septicemia
    '03819',    # Other staphylococcal septicemia
    '0382',     # Pneumococcal septicemia [Streptococcus pneumoniae septicemia]
    '0383',     # Septicemia due to anaerobes
    '03840',    # Septicemia due to gram-negative organism, unspecified
    '03841',    # Septicemia due to hemophilus influenzae [H. influenzae]
    '03842',    # Septicemia due to escherichia coli [E. coli]
    '03843',    # Septicemia due to pseudomonas
    '03844',    # Septicemia due to serratia
    '03849',    # Other septicemia due to gram-negative organisms
    '0388',     # Other specified septicemias
    '0389',     # Unspecified septicemia
    '78552',    # Septic shock
    '99591',    # Sepsis
    '99592',    # Severe sepsis
}

# Organ dysfunction indicators (for Sepsis-3 criteria)
ORGAN_DYSFUNCTION_ICD10 = {
    'respiratory': {
        'J80',      # Acute respiratory distress syndrome
        'J96',      # Respiratory failure, not elsewhere classified
        'J9600',    # Acute respiratory failure, unspecified whether with hypoxia or hypercapnia
        'J9601',    # Acute respiratory failure with hypoxia
        'J9602',    # Acute respiratory failure with hypercapnia
        'J9690',    # Respiratory failure, unspecified, unspecified whether with hypoxia or hypercapnia
        'J9691',    # Respiratory failure, unspecified with hypoxia
        'J9692',    # Respiratory failure, unspecified with hypercapnia'
    },
    'cardiovascular': {
        # TODO: check if these count as organi dysfunction codes
        'I95',      # Hypotension
        'I959',     # Hypotension, unspecified
    },
    'renal': {
        'N17',      # Acute kidney failure
        'N170',     # Acute kidney failure with tubular necrosis
        'N171',     # Acute kidney failure with acute cortical necrosis
        'N172',     # Acute kidney failure with medullary necrosis
        'N178',     # Other acute kidney failure
        'N179'      # Acute kidney failure, unspecified
    },
    'hepatic': {
        'K72',      # Hepatic failure, not elsewhere classified
        'K7200',    # Acute and subacute hepatic failure without coma
        'K7201',    # Acute and subacute hepatic failure with coma
        'K7210',    # Chronic hepatic failure without coma
        'K7211',    # Chronic hepatic failure with coma
        'K7290',    # Hepatic failure, unspecified without coma
        'K7291'     # Hepatic failure, unspecified with coma
    },
    'coagulation': {
        'D65',      # Disseminated intravascular coagulation [defibrination syndrome]
        'D688',     # Other specified coagulation defects
        'D689'      # Coagulation defect, unspecified
    },
    'neurological': {
        'G93',      # Other disorders of brain
        'G931',     # Anoxic brain damage, not elsewhere classified
        'R402',     # Coma
        'R404'      # Transient alteration of awareness
    },
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