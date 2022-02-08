#!/usr/bin/env python
"""
Author: Tom Edinburgh
v2: date 04/02/2022.

This script is to read the .csv files containing the sepsis definition as per
the AmsterdamUMCdb GitHub and the Sepsis-3 definition, comparing the
identification of sepsis in each case via a confusion table, which is generated
as a string for using in LaTeX. See files:
https://github.com/AmsterdamUMC/AmsterdamUMCdb/blob/master/concepts/diagnosis/reason_for_admission.ipynb
https://github.com/AmsterdamUMC/AmsterdamUMCdb/blob/master/concepts/sepsis3/reason_for_admission_sepsis.ipynb
https://github.com/tedinburgh/AmsterdamUMCdb/blob/master/concepts/sepsis3/sepsis3_amsterdamumcdb.ipynb
"""
###############################################################################

import argparse
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(
        description='''Compare the number of cases of sepsis at admission for
            for the Sepsis-3 and the previous sepsis definition.''')
    parser.add_argument(
        '--data_file_path',
        default='../../data/',
        help='''File path to the directory that contains the base
            AmsterdamUMCdb .csv files. These files are not directly available
            from the AmsterdamUMCdb GitHub page, and access must be
            specifically requested from Amsterdam UMC.
            (default: %(default)s)''',
        type=str)
    parser.add_argument(
        '--additional_file_path',
        default='../../data/additional_files/',
        help='''File path to the directory that contains the output of the
            script sepsis3_amsterdamumcdb.py (default: %(default)s)''',
        type=str)
    args = parser.parse_args()
    return args


inputs = parse_args()

###############################################################################
# Compare to AmsterdamUMCdb 'reason for admission' sepsis definition

admissions_df = pd.read_csv(inputs.data_file_path + 'admissions.csv')

combined_diagnoses = pd.read_csv(
    inputs.additional_file_path + 'combined_diagnoses.csv')

sofa = pd.read_csv(inputs.additional_file_path + 'sofa.csv')
sepsis = pd.read_csv(inputs.additional_file_path + 'sepsis3.csv')

combined_diagnoses_add = combined_diagnoses[['admissionid', 'sepsis']]
combined_diagnoses_add.rename(columns={'sepsis': 'old_sepsis'}, inplace=True)

# Sepsis-3 at admission is defined as Sepsis-3 on day 0  or day -1.
sepsis_add = sepsis.loc[(sepsis['time'] <= 0) & (sepsis['time'] >= -1)]
sepsis_add = sepsis_add[['admissionid', 'sepsis_episode']]

sepsis_patients = sepsis_add.groupby(['admissionid']).agg(
        sepsis=pd.NamedAgg(column='sepsis_episode', aggfunc='any'),
    ).reset_index()

admissions_df = pd.merge(
    admissions_df, combined_diagnoses_add, on='admissionid', how='left')
admissions_df = pd.merge(
    admissions_df, sepsis_patients, on='admissionid', how='left')

###############################################################################


def get_confusion_matrix(df, col0, col1):
    df = df.loc[df[col0].notna() & df[col1].notna()]
    index = ['T', 'F']
    cm = pd.DataFrame(None, index=index, columns=index)
    cm.loc['T', 'T'] = ((df[col0] == 1) & (df[col1] == 1)).sum()
    cm.loc['T', 'F'] = ((df[col0] == 1) & (df[col1] == 0)).sum()
    cm.loc['F', 'T'] = ((df[col0] == 0) & (df[col1] == 1)).sum()
    cm.loc['F', 'F'] = ((df[col0] == 0) & (df[col1] == 0)).sum()
    return cm


def get_latex(df):
    df = df.astype(str)
    latex_str = r'\begin{tabular}{l|' + 'r' * df.columns.size + r'}'
    latex_str += '\n\t' + r'\toprule'
    latex_str += '\n\t & ' + ' & '.join(df.columns) + r' \\'
    latex_str += '\n\t' + r'\midrule'
    for row in list(df.index):
        row_list = [x.split('.')[0] for x in df.loc[row]]
        latex_str += '\n\t' + ' & '.join(row_list) + r' \\'
    latex_str += '\n\t' + r'\bottomrule' + '\n' + r'\end{tabular}'
    latex_str = latex_str.replace('%', r'\%')
    return latex_str


cm_all = get_confusion_matrix(admissions_df, 'sepsis', 'old_sepsis')
cm_first = get_confusion_matrix(
    admissions_df[~admissions_df['patientid'].duplicated()],
    'sepsis', 'old_sepsis')

tableIa = r'\begin{tabular}{ll|rr}' + '\n\t'
tableIa += r'\multicolumn{2}{c|}{Unique first} & \multicolumn{2}{c}{Current}'
tableIa += r' \\' + '\n\t'
tableIa += r'\multicolumn{2}{c|}{admissions} & True & False \\' + '\n\t'
tableIa += r'\midrule' + '\n\t'
tableIa += r'Sepsis-3'
tableIa += ' & True & ' + ' & '.join([str(x) for x in cm_first.loc['T', :]])
tableIa += r' \\' + '\n\t'
tableIa += ' & False & ' + ' & '.join([str(x) for x in cm_first.loc['F', :]])
tableIa += r' \\' + '\n\t'
tableIa += r'\bottomrule' + '\n'
tableIa += r'\end{tabular}'

tableIb = r'\begin{tabular}{ll|rr}' + '\n\t'
tableIb += r'\multicolumn{2}{c|}{All} & ' + r'\multicolumn{2}{c}{Current}'
tableIb += r' \\' + '\n\t'
tableIb += r'\multicolumn{2}{c|}{admissions} & True & False \\' + '\n\t'
tableIb += r'\midrule' + '\n\t'
tableIb += r'Sepsis-3'
tableIb += ' & True & ' + ' & '.join([str(x) for x in cm_all.loc['T', :]])
tableIb += r' \\' + '\n\t'
tableIb += ' & False & ' + ' & '.join([str(x) for x in cm_all.loc['F', :]])
tableIb += r' \\' + '\n\t'
tableIb += r'\bottomrule' + '\n'
tableIb += r'\end{tabular}'

tableII = get_latex(sofa.loc[sofa['admissionid'].isin([0, 1, 2, 3])])
tableIII = get_latex(sepsis.loc[sepsis['admissionid'].isin([0, 1, 2, 3])])

tableII = tableII.replace('nan', 'NaN')
tableIII = tableIII.replace('nan', 'NaN')

output_txt = 'Table 1 in LaTeX:\n\n' + tableIa + '~' * 7 + '\n'
output_txt += tableIb + '\n' * 3
output_txt += 'Table 2 in LaTeX:\n\n' + tableII + '\n' * 3
output_txt += 'Table 3 in LaTeX:\n\n' + tableIII

with open(inputs.additional_file_path + 'sepsis3_latex_tables.txt', 'w') as f:
    f.write(output_txt)
