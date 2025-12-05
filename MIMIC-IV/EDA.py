"""
Exploratory Data Analysis on the MIMIC-IV ICU Cohort
"""

from EDA_sepsis import main_sepsis
from EDA_vaso_and_iv import main_vaso_and_iv

def main():
    main_sepsis()
    main_vaso_and_iv()

if __name__ == "__main__":
    main()
