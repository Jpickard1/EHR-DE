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
    'respiratory_rate': [220045, 220210, 224688, 224689, 224690],
    # Note: 220045 is actually Heart Rate but often co-occurs with RR
    # 224688 = Respiratory Rate (Set) - ventilator setting
    # 224689 = Respiratory Rate (spontaneous)
    # 224690 = Respiratory Rate (Total)
    
    # Mean arterial pressure (ICU)
    'map': [220052, 220181, 225312],
    # 220052 = Arterial Blood Pressure mean (invasive)
    # 220181 = Non Invasive Blood Pressure mean
    # 225312 = ART BP Mean
    
    # Systolic blood pressure (ICU)
    'sbp': [220050, 220179],
    # 220050 = Arterial Blood Pressure systolic (invasive)
    # 220179 = Non Invasive Blood Pressure systolic
    
    # Glasgow Coma Scale components (ICU)
    'gcs_eye': [220739],      # GCS - Eye Opening
    'gcs_verbal': [223900],   # GCS - Verbal Response
    'gcs_motor': [223901],    # GCS - Motor Response
    'gcs_total': [198],       # GCS Total (if available, otherwise sum components)
    
    # Oxygen measurements (ICU)
    'fio2': [223835],         # Inspired O2 Fraction
    'pao2': [220224],         # Arterial O2 pressure (from ABG)
    'spo2': [220277, 220227], # O2 saturation pulseoxymetry, Arterial O2 Saturation
    
    # Laboratory values (from d_labitems - these are 5-digit codes)
    'bilirubin': [50885],     # Bilirubin, Total (mg/dL)
    'creatinine': [50912],    # Creatinine (mg/dL)
    'platelets': [51265],     # Platelet Count (K/uL)
    
    # Urine output (ICU output events)
    'urine_output': [226559, 226560, 226561, 226563, 226564, 226584],
    # 226559 = Foley
    # 226560 = Void
    # 226561 = Condom Cath
    # 226563 = Suprapubic
    # 226564 = R Nephrostomy
    # 226584 = Ileoconduit
}