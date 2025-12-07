"""
SOFA and qSOFA Score Calculator for MIMIC-IV
=============================================
Computes Sepsis-related Failure Assessment (SOFA) and quick SOFA (qSOFA) 
scores at each time point from patient trajectory data.

Sepsis-related Organ Failure Assesment
--------------------------------------
The SOFA (Sequential Organ Failure Assessment) score is used to track the
extent of a patient's organ dysfunction over time in an intensive care unit
(ICU) setting, helping predict mortality risk based on six organ systems. The
qSOFA (quick SOFA) score is a simpler bedside prompt that uses three clinical
criteria (respiratory rate, blood pressure, and mental status) to rapidly
identify patients with suspected infection outside the ICU who are at a higher
risk of poor outcomes.

See Vincent et al, 1996: https://github.com/Jpickard1/EHR-DE/blob/main/docs/references/41786868.pdf

IMPORTANT NOTES FOR MIMIC-IV v3.1:
----------------------------------
1. ItemIDs have changed from MIMIC-III and earlier MIMIC-IV versions
2. ICU measurements (vitals, assessments) use d_items.csv.gz
3. Lab measurements use d_labitems.csv.gz
4. GCS Total (itemid 198) may not be present; calculate from components if needed
5. Some measurements may be in chartevents OR labevents depending on context
   (e.g., PaO2 from ABG is typically in labevents, not chartevents)

Verified ItemIDs for MIMIC-IV v3.1 - see ITEMID_MAP below.
"""

import pandas as pd
import numpy as np
from datetime import datetime, timedelta
from typing import Dict, List, Tuple, Optional
import logging
from mimic_constants import ITEMID_MAP, get_weight_kg
from clinical_constants import VASOPRESSORS

logger = logging.getLogger(__name__)

# ============================================================================
# SOFA SCORING FUNCTIONS
# ============================================================================

def calculate_pao2_fio2_ratio(pao2: float, fio2: float) -> Optional[float]:
    """
    Calculate PaO2/FiO2 ratio.
    
    Args:
        pao2: Arterial oxygen partial pressure (mmHg)
        fio2: Fraction of inspired oxygen (0.21-1.0 or 21-100)
    
    Returns:
        PaO2/FiO2 ratio or None if invalid
    """
    if pd.isna(pao2) or pd.isna(fio2):
        return None
    
    # Convert FiO2 to fraction if given as percentage
    if fio2 > 1:
        fio2 = fio2 / 100
    
    if fio2 < 0.21 or fio2 > 1.0:
        return None
    
    return pao2 / fio2


def calculate_spo2_fio2_ratio(spo2: float, fio2: float) -> Optional[float]:
    """
    Calculate SpO2/FiO2 ratio (surrogate for PaO2/FiO2).
    
    Args:
        spo2: Oxygen saturation (%)
        fio2: Fraction of inspired oxygen (0.21-1.0 or 21-100)
    
    Returns:
        SpO2/FiO2 ratio or None if invalid
    """
    if pd.isna(spo2) or pd.isna(fio2):
        return None
    
    # Convert FiO2 to fraction if given as percentage
    if fio2 > 1:
        fio2 = fio2 / 100
    
    if fio2 < 0.21 or fio2 > 1.0 or spo2 < 0 or spo2 > 100:
        return None
    
    return spo2 / fio2


def sofa_respiratory(pao2_fio2: Optional[float], 
                     spo2_fio2: Optional[float],
                     mechanical_ventilation: bool = False) -> int:
    """
    SOFA Respiratory score (0-4).
    
    Criteria:
    0: PaO2/FiO2 >= 400
    1: PaO2/FiO2 < 400
    2: PaO2/FiO2 < 300
    3: PaO2/FiO2 < 200 with ventilation
    4: PaO2/FiO2 < 100 with ventilation
    
    # TODO: mechanical ventilation is not currently checked for

    Uses SpO2/FiO2 as surrogate if PaO2/FiO2 unavailable.
    """
    ratio = pao2_fio2
    
    # Use SpO2/FiO2 if PaO2/FiO2 unavailable (Rice et al. 2007 conversion)
    #   Comparison of the Spo2/Fio2 Ratio and the Pao2/Fio2 Ratio in Patients With Acute Lung Injury or ARDS
    #   TODO: verify these conversions between Spo2/Fio2 to Pao2/Fio2
    if ratio is None and spo2_fio2 is not None:
        # Approximate conversion: SpO2/FiO2 <= 235 ~ PaO2/FiO2 <= 200
        if spo2_fio2 >= 512:
            ratio = 400  # Approximate normal
        elif spo2_fio2 >= 357:
            ratio = 350
        elif spo2_fio2 >= 235:
            ratio = 250
        elif spo2_fio2 >= 150:
            ratio = 150
        else:
            ratio = 75
    
    if ratio is None:
        return 0  # Cannot calculate
    
    if ratio >= 400:
        return 0
    elif ratio >= 300:
        return 1
    elif ratio >= 200:
        return 2
    elif ratio >= 100:
        return 3 if mechanical_ventilation else 2
    else:
        return 4 if mechanical_ventilation else 3


def sofa_coagulation(platelets: Optional[float]) -> int:
    """
    SOFA Coagulation score (0-4).
    
    Criteria (platelets x 10³/μL):
    0: >= 150
    1: < 150
    2: < 100
    3: < 50
    4: < 20
    """
    if pd.isna(platelets):
        return 0
    
    if platelets >= 150:
        return 0
    elif platelets >= 100:
        return 1
    elif platelets >= 50:
        return 2
    elif platelets >= 20:
        return 3
    else:
        return 4


def sofa_liver(bilirubin: Optional[float]) -> int:
    """
    SOFA Liver score (0-4).
    
    Criteria (bilirubin mg/dL):
    0: < 1.2
    1: 1.2-1.9
    2: 2.0-5.9
    3: 6.0-11.9
    4: >= 12.0
    """
    if pd.isna(bilirubin):
        return 0
    
    if bilirubin < 1.2:
        return 0
    elif bilirubin < 2.0:
        return 1
    elif bilirubin < 6.0:
        return 2
    elif bilirubin < 12.0:
        return 3
    else:
        return 4


def sofa_cardiovascular(map_value: Optional[float],
                       vasopressor_dose: Dict[str, float] = None) -> int:
    """
    SOFA Cardiovascular score (0-4).
    
    Criteria:
    0: MAP >= 70 mmHg
    1: MAP < 70 mmHg
    2: Dopamine <= 5 or dobutamine (any dose)
    3: Dopamine 5.1-15 or epi/norepi <= 0.1
    4: Dopamine > 15 or epi/norepi > 0.1
    
    Vasopressor doses in μg/kg/min
    """
    if vasopressor_dose is None:
        vasopressor_dose = {}
    
    # Check for vasopressor use
    dopamine = vasopressor_dose.get('dopamine', 0)
    dobutamine = vasopressor_dose.get('dobutamine', 0)
    epinephrine = vasopressor_dose.get('epinephrine', 0)
    norepinephrine = vasopressor_dose.get('norepinephrine', 0)
    
    # Score 4: High-dose vasopressors
    if dopamine > 15 or epinephrine > 0.1 or norepinephrine > 0.1:
        return 4
    
    # Score 3: Moderate-dose vasopressors
    if (5 < dopamine <= 15) or (0 < epinephrine <= 0.1) or (0 < norepinephrine <= 0.1):
        return 3
    
    # Score 2: Low-dose dopamine or any dobutamine
    if dopamine > 0 or dobutamine > 0:
        return 2
    
    # Score based on MAP
    if pd.isna(map_value):
        return 0
    
    if map_value < 70:
        return 1
    
    return 0


def sofa_neurological(gcs: Optional[float]) -> int:
    """
    SOFA Neurological score (0-4).
    
    Criteria (Glasgow Coma Scale):
    0: GCS 15
    1: GCS 13-14
    2: GCS 10-12
    3: GCS 6-9
    4: GCS < 6
    """
    if pd.isna(gcs):
        return 0
    
    if gcs == 15:
        return 0
    elif gcs >= 13:
        return 1
    elif gcs >= 10:
        return 2
    elif gcs >= 6:
        return 3
    else:
        return 4


def sofa_renal(creatinine: Optional[float], 
               urine_output: Optional[float]) -> int:
    """
    SOFA Renal score (0-4).
    
    Criteria:
    0: Creatinine < 1.2 mg/dL
    1: Creatinine 1.2-1.9 mg/dL
    2: Creatinine 2.0-3.4 mg/dL
    3: Creatinine 3.5-4.9 mg/dL or UO < 500 mL/day
    4: Creatinine >= 5.0 mg/dL or UO < 200 mL/day
    
    Args:
        creatinine: Serum creatinine (mg/dL)
        urine_output: 24-hour urine output (mL)
    """
    score = 0
    
    # Score based on creatinine
    if not pd.isna(creatinine):
        if creatinine < 1.2:
            score = 0
        elif creatinine < 2.0:
            score = 1
        elif creatinine < 3.5:
            score = 2
        elif creatinine < 5.0:
            score = 3
        else:
            score = 4
    
    # Adjust based on urine output if available
    if not pd.isna(urine_output):
        if urine_output < 200:
            score = max(score, 4)
        elif urine_output < 500:
            score = max(score, 3)
    
    return score


# ============================================================================
# qSOFA SCORING
# ============================================================================

def qsofa_score(respiratory_rate: Optional[float],
                sbp: Optional[float],
                gcs: Optional[float]) -> int:
    """
    Quick SOFA (qSOFA) score (0-3).
    
    Criteria (1 point each):
    - Respiratory rate >= 22/min
    - Altered mentation (GCS < 15)
    - Systolic BP <= 100 mmHg
    
    Score >= 2 suggests high risk for poor outcomes.

    See: https://www.mdcalc.com/calc/2654/qsofa-quick-sofa-score-sepsis
    """
    score = 0
    
    if not pd.isna(respiratory_rate) and respiratory_rate >= 22:
        score += 1
    
    if not pd.isna(sbp) and sbp <= 100:
        score += 1
    
    if not pd.isna(gcs) and gcs < 15:
        score += 1
    
    return score


# ============================================================================
# DATA EXTRACTION AND AGGREGATION
# ============================================================================

def extract_vital_signs(df: pd.DataFrame, time_window: str = '6H') -> pd.DataFrame:
    """
    Extract and aggregate vital signs and lab values over time windows.
    
    Args:
        df: Patient trajectory dataframe
        time_window: Time window for aggregation (e.g., '6H', '24H')
    
    Returns:
        DataFrame with aggregated values per time window
    """
    # Filter to chart events and lab events
    chart_df = df[df['event_type'].isin(['chart', 'lab'])].copy()
    
    if len(chart_df) == 0:
        return pd.DataFrame()
    
    # Use charttime as primary time, fall back to time
    chart_df['event_time'] = pd.to_datetime(
        chart_df['charttime'].fillna(chart_df['time']), 
        errors='coerce'
    )
    chart_df = chart_df.dropna(subset=['event_time'])
    
    # Extract measurements by itemid
    measurements = {}
    
    for measure_name, itemids in ITEMID_MAP.items():
        measure_df = chart_df[chart_df['itemid'].isin(itemids)].copy()
        
        if len(measure_df) > 0:
            # Use valuenum for numeric values
            measure_df = measure_df[measure_df['valuenum'].notna()]
            
            # Group by time window and take mean/min/max as appropriate
            if measure_name in ['map', 'sbp', 'respiratory_rate', 'spo2']:
                # For vitals, take mean within window
                grouped = measure_df.groupby(pd.Grouper(key='event_time', freq=time_window))['valuenum'].mean()
            elif measure_name in ['gcs_total', 'gcs_eye', 'gcs_verbal', 'gcs_motor']:
                # For GCS, take minimum (worst)
                grouped = measure_df.groupby(pd.Grouper(key='event_time', freq=time_window))['valuenum'].min()
            elif measure_name in ['weight']:
                _, measure_df = get_weight_kg(measure_df)
                grouped = measure_df.groupby(pd.Grouper(key='event_time', freq=time_window))['valuenum'].last()
                grouped = grouped.ffill().bfill()
            else:
                # For labs, take last value in window
                grouped = measure_df.groupby(pd.Grouper(key='event_time', freq=time_window))['valuenum'].last()
            
            measurements[measure_name] = grouped
    
    # Combine into single dataframe
    if measurements:
        result = pd.DataFrame(measurements)
        result.index.name = 'time'
        return result.reset_index()
    
    return pd.DataFrame()


def calculate_gcs_total(gcs_eye: Optional[float], 
                       gcs_verbal: Optional[float],
                       gcs_motor: Optional[float],
                       gcs_total: Optional[float]) -> Optional[float]:
    """Calculate total GCS from components or use direct total."""
    if not pd.isna(gcs_total):
        return gcs_total
    
    if not pd.isna(gcs_eye) and not pd.isna(gcs_verbal) and not pd.isna(gcs_motor):
        return gcs_eye + gcs_verbal + gcs_motor
    
    return None


def calculate_24h_urine_output(df: pd.DataFrame, current_time: pd.Timestamp) -> Optional[float]:
    """
    Calculate 24-hour urine output looking back from current time.
    
    Args:
        df: Patient trajectory dataframe
        current_time: Current timepoint
    
    Returns:
        Total urine output in mL over past 24 hours
    """
    # Filter urine output events
    urine_df = df[
        (df['event_type'] == 'outputevents') & 
        (df['itemid'].isin(ITEMID_MAP['urine_output']))
    ].copy()
    
    if len(urine_df) == 0:
        return None
    
    # Parse times
    urine_df['event_time'] = pd.to_datetime(urine_df['charttime'].fillna(urine_df['time']), errors='coerce')
    urine_df = urine_df.dropna(subset=['event_time'])
    
    # Get 24-hour window
    start_time = current_time - pd.Timedelta(hours=24)
    window_df = urine_df[
        (urine_df['event_time'] >= start_time) & 
        (urine_df['event_time'] <= current_time)
    ]
    
    if len(window_df) == 0:
        return None
    
    # Sum amounts
    total = window_df['value'].astype(float).sum()
    return total if total > 0 else None


def extract_vasopressor_doses(df: pd.DataFrame, current_time: pd.Timestamp) -> Dict[str, float]:
    """
    Extract vasopressor doses at current timepoint.
    
    TODO: need to implement this function

    Returns doses in μg/kg/min for SOFA scoring.
    """
    # Determine the patients weight
    

    # This is a simplified version - would need patient weight and rate conversions
    # For now, return empty dict (indicating no vasopressors detected)
    return {}

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


# ============================================================================
# MAIN SCORING FUNCTION
# ============================================================================

def calculate_sofa_qsofa_scores(patient_df: pd.DataFrame, 
                                time_window: str = '6H') -> pd.DataFrame:
    """
    Calculate SOFA and qSOFA scores at each timepoint for a patient.
    
    Args:
        patient_df: Patient trajectory dataframe (single patient)
        time_window: Time window for aggregating measurements
    
    Returns:
        DataFrame with columns:
            - time: Timepoint
            - qsofa_score: qSOFA score (0-3)
            - sofa_respiratory: Respiratory component (0-4)
            - sofa_coagulation: Coagulation component (0-4)
            - sofa_liver: Liver component (0-4)
            - sofa_cardiovascular: Cardiovascular component (0-4)
            - sofa_neurological: Neurological component (0-4)
            - sofa_renal: Renal component (0-4)
            - sofa_total: Total SOFA score (0-24)
    """
    # Extract vital signs and labs
    vitals_df = extract_vital_signs(patient_df, time_window=time_window)
    
    if len(vitals_df) == 0:
        logger.warning("No vital signs found for SOFA/qSOFA calculation")
        return pd.DataFrame()
    
    # Calculate scores for each timepoint
    scores = []
    
    for idx, row in vitals_df.iterrows():
        current_time = row['time']
        
        # Calculate GCS total
        gcs = calculate_gcs_total(
            row.get('gcs_eye'),
            row.get('gcs_verbal'),
            row.get('gcs_motor'),
            row.get('gcs_total')
        )
        
        # Calculate PaO2/FiO2 and SpO2/FiO2 ratios
        pao2_fio2 = calculate_pao2_fio2_ratio(
            row.get('pao2'),
            row.get('fio2')
        )
        spo2_fio2 = calculate_spo2_fio2_ratio(
            row.get('spo2'),
            row.get('fio2')
        )
        
        # Get 24-hour urine output
        urine_24h = calculate_24h_urine_output(patient_df, current_time)
        
        # TODO: this does not work
        # Get vasopressor doses
        vasopressor_doses = extract_vasopressor_doses(patient_df, current_time)
        
        # Calculate qSOFA
        qsofa = qsofa_score(
            row.get('respiratory_rate'),
            row.get('sbp'),
            gcs
        )
        
        # Calculate SOFA components
        sofa_resp = sofa_respiratory(pao2_fio2, spo2_fio2)
        sofa_coag = sofa_coagulation(row.get('platelets'))
        sofa_liv = sofa_liver(row.get('bilirubin'))
        sofa_cardio = sofa_cardiovascular(row.get('map'), vasopressor_doses)
        sofa_neuro = sofa_neurological(gcs)
        sofa_ren = sofa_renal(row.get('creatinine'), urine_24h)
        
        sofa_total = (sofa_resp + sofa_coag + sofa_liv + 
                     sofa_cardio + sofa_neuro + sofa_ren)
        
        scores.append({
            'time': current_time,
            'qsofa_score': qsofa,
            'sofa_respiratory': sofa_resp,
            'sofa_coagulation': sofa_coag,
            'sofa_liver': sofa_liv,
            'sofa_cardiovascular': sofa_cardio,
            'sofa_neurological': sofa_neuro,
            'sofa_renal': sofa_ren,
            'sofa_total': sofa_total,
            # Include raw values for reference
            'gcs': gcs,
            'map': row.get('map'),
            'respiratory_rate': row.get('respiratory_rate'),
            'sbp': row.get('sbp'),
            'pao2_fio2': pao2_fio2,
            'spo2_fio2': spo2_fio2,
            'platelets': row.get('platelets'),
            'bilirubin': row.get('bilirubin'),
            'creatinine': row.get('creatinine'),
            'urine_24h': urine_24h,
        })
    
    return pd.DataFrame(scores)
