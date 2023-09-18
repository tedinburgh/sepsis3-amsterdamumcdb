# Sepsis-3 epidemiology in AmsterdamUMCdb

This repository contains files for reproducing figures and tables in a manuscript describing sepsis epidemiology in the freely-accessible Amsterdam University Medical Centers Database. This manuscript is currently under peer review. The scripts in this repository update those in `../sepsis3_gigabyte'.

To run this code, download the repository and simply move the folder `concepts/sepsis3` to the `concepts` folder in your local version of AmsterdamUMCdb and run all python scripts in this order, changing `ADB_PATH` in the below to the relevant file path for AmsterdamUMCdb, as follows:

```
mv concepts/sepsis3 ADB_PATH/concepts/
mkdir ADB_PATH/data/additional_files
mkdir ADB_PATH/data/additional_files/exclude_noradrenaline_6hr
cd ADB_PATH/concepts/
python reason_for_admission.py
python sepsis3_amsterdamumcdb.py
python sepsis3_amsterdamumcdb.py --exclude_noradrenaline --output_file_path '../../data/additional_files/exclude_noradrenaline_6hr'
python sepsis3_tables.py
python sepsis3_figures.py
```

This will create a number of additional .csv files in a new folder `data/additional_files`. These data files will be created in this repository:
+ `sofa.csv` contains daily SOFA scores (component and total scores) for each admission. This is based on [AmsterdamUMCdb SOFA score SQL code](https://github.com/AmsterdamUMC/AmsterdamUMCdb/blob/master/concepts/severityscores/sofa.ipynb), but is written using pandas from the base .csv data files.
+ `sepsis3.csv` contains total SOFA score, antibiotic escalation, prophylactic antibiotic usage, suspected infection and Sepsis-3 diagnosis (including septic shock) for each admission, daily.
+ `combined_diagnoses.csv` is a pandas rewrite of the [AmsterdamUMCdb reason for admission SQL code](https://github.com/AmsterdamUMC/AmsterdamUMCdb/blob/master/concepts/diagnosis/reason_for_admission.ipynb).
+ `sepsis_all.csv' contains more detailed information for each admission, with all variables involved in computing sepsis incidence.
+ `sofa_platelets.csv', `creatinine.csv', `sofa_cardiovascular.csv', `sofa_cns.csv', `sofa_bilirubin.csv', `sofa_cardiovascular_meds.csv', `sofa_renal.csv' and `sofa_respiration.csv' contain daily SOFA component scores, following [AmsterdamUMCdb SOFA score SQL code](https://github.com/AmsterdamUMC/AmsterdamUMCdb/blob/master/concepts/severityscores/sofa.ipynb).
+ `drugitems_abx.csv' contains all instances of non-prophylactic antibiotic administration, which informs the antibiotic escalation definition. 
+ `drugitems_abx_all.csv' contains all instances of antibiotic administration, including prophylactic antibiotics.
+ `tables/' contains comma-separated Tables included in the accompanying manuscript, plus additional unused tables.
+ `figures/' contains figures included in the accompanying manuscript, plus additional unused figures. 
