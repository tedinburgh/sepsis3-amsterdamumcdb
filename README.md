# Sepsis-3 criteria in AmsterdamUMCdb

This repository contains files for implementing the Sepsis-3 definition in the freely-accessible Amsterdam University Medical Centers Database.

You must first have access to AmsterdamUMCdb, including the data. Further instructions for this can be found on the [Amsterdam Medical Data Science](https://amsterdammedicaldatascience.nl/) website or [AmsterdamUMCdb GitHub](https://github.com/AmsterdamUMC/AmsterdamUMCdb). To run this code, download the repository and simply move the folder `concepts/sepsis3` to the `concepts` folder in your local version of AmsterdamUMCdb and run all python scripts in this order, changing `ADB_PATH` in the below to the relevant file path for AmsterdamUMCdb, as follows:

```
mv concepts/sepsis3 ADB_PATH/concepts/
mkdir ADB_PATH/data/additional_files
cd ADB_PATH/concepts/
python reason_for_admission.py
python sepsis3_amsterdamumcdb.py
python sepsis_comparison.py
```

This will create additional .csv files in a new folder `data/additional_files`. These data files cannot be included in this repository 
+ `sofa.csv` contains daily SOFA scores (component and total scores) for each admission. This is based on [AmsterdamUMCdb SOFA score SQL code](https://github.com/AmsterdamUMC/AmsterdamUMCdb/blob/master/concepts/severityscores/sofa.ipynb), but is written using pandas from the base .csv data files.
+ `sepsis3.csv` contains total SOFA score, antibiotic escalation, prophylactic antibiotic usage, suspected infection and Sepsis-3 diagnosis (including septic shock) for each admission, daily.
+ `combined_diagnoses.csv` is a pandas rewrite of the [AmsterdamUMCdb reason for admission SQL code](https://github.com/AmsterdamUMC/AmsterdamUMCdb/blob/master/concepts/diagnosis/reason_for_admission.ipynb).
+ `sepsis3_latex_tables.txt` is a text file containing the LaTeX code for tables included in the manuscript accompanying this work.

The database is primarily in Dutch. It is generally easier to work with itemids rather than the (Dutch) item names, in any table. The database contains a lookup dictionary of all terms, instructions for installing the `amsterdamumcdb` package that contains this is given [here](https://github.com/AmsterdamUMC/AmsterdamUMCdb/tree/master/setup-amsterdamumcdb). It may be worth familiarising yourself with some common terms (e.g. Neurochirurgie = Neurosurgery, Vaatchirurgie = Vascular Surgery) ahead of time.

We encourage other researchers (or anyone interested in this code from a data or clinical perspective) to submit an [issue](https://github.com/tedinburgh/sepsis3-amsterdamumcdb/issues) on this repository if they would like to report or fix any bugs or to contribute improvements.
