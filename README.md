# Data Engineering for Electronic Health Records

Data Engineering for Electronic Health Records

## MIMIC: Medical Information Mart for Intensive Care

**1. Data Access:**
```
cd <path/to/data>
wget -r -N -c -np --user <username> --ask-password https://physionet.org/files/mimiciv/3.1/
```

**2. Processing Patient Trajectories**

The goal of this step is to transform the raw MIMIC-IV CSV tables into **patient-level trajectory files**. A *patient trajectory* is a chronological sequence of clinical events (labs, vitals, diagnoses, medications, ICU stays, etc.) grouped by patient and prepared for downstream machine learning or analytic tasks.

Before running the processing script, you must update the configuration values inside:


*Run As*
```
tmux
conda activate mimic4
cd MIMIC-IV
python build_patient_files.py
```

Required Configuration

- **`DATA_DIR`**  
  The absolute path to the folder where your downloaded MIMIC-IV dataset is stored.  
  This directory should contain subfolders such as `core/`, `hosp/`, `icu/`, `note/`, etc.

- **`TRAJECTORY_VERSION`**  
  A string identifier for this processing run (e.g., `"v1"` or `"exp_2025_01"`).  
  This is useful for tracking different preprocessing pipelines or versions of the output trajectory format.

- **`MAX_CPUS`** (default: `20`)  
  Controls the number of parallel worker processes.  
  Increase this value if you have many CPU cores; decrease it if memory is limited.

*What the Script Does*

1. **Loads raw MIMIC-IV tables**  
   Reads relevant CSV files from `DATA_DIR`, including patient demographics, admissions, ICU stays, labs, charted vitals, procedures, diagnoses, and medication administrations.

2. **Merges data across tables**  
   Joins tables by `subject_id`, `hadm_id`, and `stay_id` to combine multiple sources of clinical information.

3. **Cleans and standardizes event fields**  
   Ensures timestamps, units, and codes are consistent and usable for downstream analysis.

4. **Writes output trajectory files**  
   Saves the processed trajectories—typically one file per patient or per admission—tagged using your chosen `TRAJECTORY_VERSION`.


