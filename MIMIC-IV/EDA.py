"""
Exploratory Data Analysis on the MIMIC-IV ICU Cohort
"""

from EDA_sepsis import main_sepsis
from EDA_vaso_and_iv import main_vaso_and_iv
from EDA_times import main_time

def main():
    main_sepsis()
    main_vaso_and_iv()
    main_time()

if __name__ == "__main__":
    main()
