## Author: Tom Edinburgh
## v1: date 03/10/2021.

## This script is to read the .csv files containing the sepsis definition
## as per the AmsterdamUMCdb GitHub and the Sepsis-3 definition, comparing the
## identification of sepsis in each case via a confusion table, which is
## generated as a string for using in LaTeX.
## https://github.com//AmsterdamUMC//AmsterdamUMCdb/blob/master/concepts/diagnosis/reason_for_admission.ipynb
## https://github.com//AmsterdamUMC//AmsterdamUMCdb/blob/master/concepts/sepsis3/reason_for_admission_sepsis.ipynb
## https://github.com/tedinburgh/AmsterdamUMCdb/blob/master/concepts/sepsis3/sepsis3_amsterdamumcdb.ipynb

################################################################################

import pandas as pd
import numpy as np
import re

## The directory 'data' should contain the base AmsterdamUMCdb .csv files.
## This is not directly available from the AmsterdamUMCdb GitHub page, and
## access must be specifically requested from the database owners. The
## directory 'data/additional_files' will contain output .csv files from this
## script, as well as the ATC codes and antibiotics .csv files mentioned later.
file_path = '../../data/'
additional_file_path = '../../data/additional_files/'

admissions_df = pd.read_csv(file_path + 'admissions.csv')

################################################################################
## Compare to AmsterdamUMCdb 'reason for admission' sepsis definition

combined_diagnoses = pd.read_csv(
    additional_file_path + 'combined_diagnoses.csv')

sofa = pd.read_csv(additional_file_path + 'sofa.csv')
sepsis = pd.read_csv(additional_file_path + 'sepsis3.csv')

combined_diagnoses_add = combined_diagnoses[['admissionid', 'sepsis']]
combined_diagnoses_add.rename(columns={'sepsis': 'old_sepsis'}, inplace=True)

sepsis_add = sepsis.loc[(sepsis['time'] <= 0) & (sepsis['time'] >= -1)]
sepsis_add = sepsis_add[['admissionid', 'sepsis_episode']]

sepsis_patients = sepsis_add.groupby(['admissionid']).agg(
        sepsis=pd.NamedAgg(column='sepsis_episode', aggfunc='any'),
    ).reset_index()

admissions_df = pd.merge(admissions_df, combined_diagnoses_add,
    on='admissionid', how='left')
admissions_df = pd.merge(admissions_df, sepsis_patients,
    on='admissionid', how='left')

temp_ind = (
    admissions_df['sepsis'].notna() &
    admissions_df['old_sepsis'].notna())
all_admissions_tp = (
    temp_ind &
    (admissions_df['sepsis'] == 1) &
    (admissions_df['old_sepsis'] == 1)).sum()
all_admissions_fn = (
    temp_ind &
    (admissions_df['sepsis'] == 1) &
    (admissions_df['old_sepsis'] == 0)).sum()
all_admissions_fp = (
    temp_ind &
    (admissions_df['sepsis'] == 0) &
    (admissions_df['old_sepsis'] == 1)).sum()
all_admissions_tn = (
    temp_ind &
    (admissions_df['sepsis'] == 0) &
    (admissions_df['old_sepsis'] == 0)).sum()

temp_ind_first = temp_ind & ~admissions_df['patientid'].duplicated()
first_admissions_tp = (
    temp_ind_first &
    (admissions_df['sepsis'] == 1) &
    (admissions_df['old_sepsis'] == 1)).sum()
first_admissions_fn = (
    temp_ind_first &
    (admissions_df['sepsis'] == 1) &
    (admissions_df['old_sepsis'] == 0)).sum()
first_admissions_fp = (
    temp_ind_first &
    (admissions_df['sepsis'] == 0) &
    (admissions_df['old_sepsis'] == 1)).sum()
first_admissions_tn = (
    temp_ind_first &
    (admissions_df['sepsis'] == 0) &
    (admissions_df['old_sepsis'] == 0)).sum()

def get_latex(df):
    df = df.astype(str)
    latex_str = r'\begin{table}' + '\n\t' + r'\begin{tabular}{l|'
    latex_str += 'r' * df.columns.size + r'}'
    latex_str += '\n\t\t' + r'\toprule'
    latex_str += '\n\t\t & ' + ' & '.join(df.columns) + r' \\'
    latex_str += '\n\t\t' + r'\midrule'
    for row in list(df.index):
        row_list = [x.split('.')[0] for x in df.loc[row]]
        latex_str += '\n\t\t' + ' & '.join(row_list) + r' \\'
    latex_str += '\n\t\t' + r'\bottomrule' + '\n\t' + r'\end{tabular}'
    latex_str += '\n' + r'\end{table}'
    latex_str = latex_str.replace('%', '\%')
    return latex_str

tableI = r'\begin{table}' + '\n\t'
tableI += r'\begin{tabular}{l|l|rr|r}' + '\n\t\t'
tableI += r'\multicolumn{2}{c|}{Unique first} & \multicolumn{2}{c}{Current}'
tableI += r' \\' +  '\n\t\t'
tableI += r'\multicolumn{2}{c|}{admissions} & True & False & \\' + '\n\t\t'
tableI += r'\midrule' + '\n\t\t'
tableI += r'Sepsis-3 & True & ' + str(first_admissions_tp)
tableI += ' & ' + str(first_admissions_fn) + r' \\' + '\n\t\t'
tableI += ' & False & ' + str(first_admissions_fp)
tableI += ' & ' + str(first_admissions_tn) + r' \\' + '\n\t\t'
tableI += r'\bottomrule' + '\n\t'
tableI += r'\end{tabular}~~~~~~~~' + '\n\t'
tableI += r'\begin{tabular}{l|l|rr|r}' + '\n\t\t'
tableI += r'\multicolumn{2}{c|}{All} & ' + r'\multicolumn{2}{c}{Current}'
tableI += r' \\' +  '\n\t\t'
tableI += r'\multicolumn{2}{c|}{admissions} & True & False & \\' + '\n\t\t'
tableI += r'\midrule' + '\n\t\t'
tableI += r'Sepsis-3 & True & ' + str(all_admissions_tp)
tableI += ' & ' + str(all_admissions_fn) + r' \\' + '\n\t\t'
tableI += ' & False & ' + str(all_admissions_fp) + r' \\' + '\n\t\t'
tableI += r'\bottomrule' + '\n\t'
tableI += r'\end{tabular}' + '\n'
tableI += r'\end{table}'

tableII = get_latex(sofa.loc[sofa['admissionid'].isin([0,1,2,3])])
tableIII = get_latex(sepsis.loc[sepsis['admissionid'].isin([0,1,2,3])])

output_txt = 'Table 1 in LaTeX:\n\n' + tableI + '\n\n\n'
output_txt += 'Table 2 in LaTeX:\n\n' + tableII + '\n\n\n'
output_txt += 'Table 3 in LaTeX:\n\n' + tableIII + '\n'

with open(additional_file_path + 'sepsis3_latex_tables.txt', 'w') as f:
    f.write(output.txt)
