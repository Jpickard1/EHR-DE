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