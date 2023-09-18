# Sepsis-3 criteria in AmsterdamUMCdb

This repository contains code accompanying [our manuscript](http://dx.doi.org/10.46471/gigabyte.45) published in the journal Gigabyte.

To run this code, download the AmsterdamUMCdb repository and move the folder `concepts/sepsis3` to the `concepts` folder in your local version of AmsterdamUMCdb. Then run all python scripts in this order, changing `ADB_PATH` in the below to the relevant file path for AmsterdamUMCdb, as follows:

```
mv sepsis3_gigabyte/concepts/sepsis3 ADB_PATH/concepts/
mkdir ADB_PATH/data/additional_files
cd ADB_PATH/concepts/
python reason_for_admission.py
python sepsis3_amsterdamumcdb.py
python sepsis_comparison.py
```

This will create additional .csv files in a new folder `data/additional_files`. These data files will be created in this repository:
+ `sofa.csv` contains daily SOFA scores (component and total scores) for each admission. This is based on [AmsterdamUMCdb SOFA score SQL code](https://github.com/AmsterdamUMC/AmsterdamUMCdb/blob/master/concepts/severityscores/sofa.ipynb), but is written using pandas from the base .csv data files.
+ `sepsis3.csv` contains total SOFA score, antibiotic escalation, prophylactic antibiotic usage, suspected infection and Sepsis-3 diagnosis (including septic shock) for each admission, daily.
+ `combined_diagnoses.csv` is a pandas rewrite of the [AmsterdamUMCdb reason for admission SQL code](https://github.com/AmsterdamUMC/AmsterdamUMCdb/blob/master/concepts/diagnosis/reason_for_admission.ipynb).
+ `sepsis3_latex_tables.txt` is a text file containing the LaTeX code for tables included in the manuscript accompanying this work.

This will reproduce the process described in the Gigabyte manuscript, as well as the Tables contained in this manuscript. Our implementation of the Sepsis-3 criteria was slightly updated in our work on sepsis epidemiology in AmsterdamUMCdb. The original code is included in this folder for posterity, but please consider using the updated python scripts in the folder `../sepsis3_epidemiology/concepts/sepsis3`. This folder includes a summary of the changes made to the implementation of the Sepsis-3 criteria.
