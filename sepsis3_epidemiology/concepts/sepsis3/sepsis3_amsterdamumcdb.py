#!/usr/bin/env python
'''
Author: Tom Edinburgh
v4: date 31/08/2022.

This script is to calculate SOFA scores from the AmsterdamUMCdb database and
write them to a csv file (that also contains other admission information).

It's adapted from this AmsterdamUMCdb script, which calculates SOFA scores:
https://github.com/AmsterdamUMC/AmsterdamUMCdb/blob/master/concepts/severityscores/sofa.ipynb)
But this code using pandas, working directly from the AmsterdamUMCdb .csv
data files, instead of via SQL queries (as in the GitHub above). These SQL
queries are at the bottom of the script for reference (but aren't used in the
script).
'''

###############################################################################

import pandas as pd
import numpy as np
import re
import argparse
import amsterdamumcdb as adb

###############################################################################


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
        help='''File path to the directory that will contain the output of
            intermediary scripts, which should be the same output file path as
            in the script reason_for_admission.py (default: %(default)s)''',
        type=str)
    parser.add_argument(
        '--output_file_path',
        default='../../data/additional_files/',
        help='''File path to the directory that will contain the output of this
            script (default: %(default)s)''',
        type=str)
    parser.add_argument(
        '--MIN_TIME',
        default=0,
        help='''In order to ease computation time, we limit the range of days
            in which we identify Sepsis-3 episodes to >=MIN_TIME days. This can
            be set to 'None' to include all data but doing so may result in
            memory costs that kill your session (default: %(default)s)''',
        type=int)
    parser.add_argument(
        '--MAX_TIME',
        default=14,
        help='''In order to ease computation time, we limit the range of days
            in which we identify Sepsis-3 episodes to >=MAX_TIME days. This can
            be set to 'None' to include all data but doing so may result in
            memory costs that kill your session (default: %(default)s)''',
        type=int)
    parser.add_argument(
        '--exclude_noradrenaline',
        default=None,
        help='''As a sensitivity analysis, we may wish to exclude from
            medication data any administration of noradrenaline that lasts less
            than e.g. 360 mins (6 hrs). If set to a number, this excludes
            noradrenaline administration of less than this number of minutes
            (default: %(default)s)''',
        type=float)
    parser.add_argument(
        '--suspected_infection_by_cultures',
        default=True,
        help='''As a sensitivity analysis, we may wish to identify suspected
            infections by the clinician administering antibiotics and taking
            cultures, instead of using antibiotic escalation for defining a
            suspected infection. This is likely to lead to false negatives,
            since in practice the clinician may decide against taking a
            culture, if they are sure it is sepsis (cultures are specific but
            not particularly sensitive, so the clinician may just keep a
            broader spectrum antibiotic). (default: %(default)s)''',
        type=bool)
    parser.add_argument(
        '--save_intermediate_files',
        default=True,
        help='''Save intermediate .csv files for easier further analysis
        (default: %(default)s)''',
        type=bool)
    args = parser.parse_args()
    return args


inputs = parse_args()
MAX_TIME = inputs.MAX_TIME
MIN_TIME = inputs.MIN_TIME

###############################################################################
# Load the data


def numerics_read(itemid_list, start=None, end=None, admissions_df=None):
    '''
    Load numerics by chunk and retains rows with itemid in the given list. If
    either start/end are given as inputs to this function, then only rows with
    start < t < end are retained (where t refers to time since admission if
    admissions_df is given, else is just the given measurement time)
    '''
    numerics_list = []
    ii = 0
    if admissions_df is not None:
        admissions_df = admissions_df[['admissionid', 'admittedat']]
    file_name = inputs.data_file_path + 'numericitems.csv'
    with pd.read_csv(file_name, **numerics_csv) as reader:
        for chunk in reader:
            if ((ii % 100) == 0):
                print(ii)
            chunk = chunk.loc[chunk['itemid'].isin(itemid_list)]
            if admissions_df is not None:
                chunk = pd.merge(
                    chunk, admissions_df,
                    on='admissionid', how='left')
                chunk['diff_measuredat'] = (
                    chunk['measuredat'] - chunk['admittedat'])
            else:
                chunk['diff_measuredat'] = chunk['measuredat']
            if start is not None:
                chunk = chunk.loc[chunk['diff_measuredat'] >= start]
            if end is not None:
                chunk = chunk.loc[chunk['diff_measuredat'] <= end]
            numerics_list.append(chunk)
            ii += 1
    numerics = pd.concat(numerics_list)
    # This is the 'day' since admission, i.e. 24hr period zeroed at admission.
    numerics['time'] = numerics['diff_measuredat'] // (1000*60*60*24)
    numerics.drop(columns=['diff_measuredat'], inplace=True)
    numerics.reset_index(drop=True, inplace=True)
    return numerics


dictionary = adb.get_dictionary()

list_columns = ['admissionid', 'itemid', 'valueid', 'value']
list_columns += ['updatedat', 'measuredat', 'registeredby']
listitems = pd.read_csv(
    inputs.data_file_path + 'listitems.csv',
    usecols=list_columns, encoding='latin-1')

drug_columns = ['admissionid', 'itemid', 'item', 'duration', 'rate']
drug_columns += ['rateunit', 'start', 'stop', 'dose', 'doserateperkg']
drug_columns += ['doseunitid', 'doserateunitid', 'ordercategoryid']
drugitems = pd.read_csv(
    inputs.data_file_path + 'drugitems.csv',
    usecols=drug_columns, encoding='latin-1')

freetextitems = pd.read_csv(
    inputs.data_file_path + 'freetextitems.csv', encoding='latin-1')

procedureorderitems = pd.read_csv(
    inputs.data_file_path + 'procedureorderitems.csv', encoding='latin-1')

admissions_df = pd.read_csv(
    inputs.data_file_path + 'admissions.csv', encoding='latin-1')

# From the script 'reason_for_admission.py', which must be run first
combined_diagnoses = pd.read_csv(
    inputs.additional_file_path + 'combined_diagnoses.csv')

# Out of the AmsterdamUMCdb data files, numerics is by far the biggest. So
# we have to be a bit more careful loading it in (by chunks).
numerics_columns = ['admissionid', 'itemid', 'item', 'value', 'unitid']
numerics_columns += ['measuredat', 'registeredby', 'islabresult', 'fluidout']
numerics_dtypes = ['int64', 'int64', 'str', 'float64', 'int64']
numerics_dtypes += ['int64', 'str', 'bool', 'float64']
numerics_dtypes = dict(zip(numerics_columns, numerics_dtypes))
numerics_csv = dict(
    encoding='latin-1', usecols=numerics_columns, dtype=numerics_dtypes,
    chunksize=10**6)

# This list is taken from the AmsterdamUMCdb SOFA scores script.
numerics_sofa_itemid = [8845]  # O2 l/min
numerics_sofa_itemid += [10387]  # Zuurstof toediening (bloed)
numerics_sofa_itemid += [18587]  # Zuurstof toediening
numerics_sofa_itemid += [6699]  # FiO2 %: setting on Evita ventilator
numerics_sofa_itemid += [12279]  # 12279 O2 concentratie Servo-i/Servo-U vent.
numerics_sofa_itemid += [12369]  # SET %O2: used with BiPap Vision ventilator
numerics_sofa_itemid += [16246]  # Zephyros FiO2: Non-invasive ventilation
numerics_sofa_itemid += [8794]  # UrineCAD
numerics_sofa_itemid += [8796]  # UrineSupraPubis
numerics_sofa_itemid += [8798]  # UrineSpontaan
numerics_sofa_itemid += [8800]  # UrineIncontinentie
numerics_sofa_itemid += [8803]  # UrineUP
numerics_sofa_itemid += [10743]  # Nefrodrain li Uit
numerics_sofa_itemid += [10745]  # Nefrodrain re Uit
numerics_sofa_itemid += [19921]  # UrineSplint Li
numerics_sofa_itemid += [19922]  # UrineSplint Re]
numerics_sofa_itemid += [6846]  # PCO2
numerics_sofa_itemid += [9990]  # pCO2 (bloed)
numerics_sofa_itemid += [21213]  # PCO2 (bloed) - kPa
numerics_sofa_itemid += [7433]  # PO2
numerics_sofa_itemid += [9996]  # PO2 (bloed)
numerics_sofa_itemid += [21214]  # PO2 (bloed) - kPa
numerics_sofa_itemid += [9964]  # Thrombo's (bloed)
numerics_sofa_itemid += [6797]  # Thrombocyten
numerics_sofa_itemid += [10409]  # Thrombo's citr. bloed (bloed)
numerics_sofa_itemid += [14252]  # Thrombo CD61 (bloed)
numerics_sofa_itemid += [6813]  # Bili Totaal
numerics_sofa_itemid += [9945]  # Bilirubine (bloed)

numerics_abp_itemid = [6642]  # ABP gemiddeld
numerics_abp_itemid += [6679]  # Niet invasieve bloeddruk gemiddeld
numerics_abp_itemid += [8843]  # ABP gemiddeld II

numerics_creatinine_itemid = [6836]  # 6836: Kreatinine µmol/l ...
# (erroneously documented as µmol)
numerics_creatinine_itemid += [9941]  # Kreatinine (bloed) µmol/l
numerics_creatinine_itemid += [14216]  # KREAT enzym. (bloed) µmol/l

# Convert specified times from days to milliseconds (as required by data)
end_time = 1000*60*60*24*MAX_TIME if MAX_TIME is not None else None
start_time = 1000*60*60*24*MIN_TIME if MIN_TIME is not None else None

# Create dataframes for SOFA score calculation, baseline Cr and lactate values
numerics_sofa = numerics_read(
    numerics_sofa_itemid + numerics_abp_itemid, admissions_df=admissions_df,
    end=end_time, start=start_time)
# We need a baseline creatinine, so look back further.
numerics_creatinine = numerics_read(
    numerics_creatinine_itemid, admissions_df=admissions_df,
    end=end_time, start=-1000*60*60*24*365)
# Lactate (for septic shock: max lactate of 2mmol/L + cardiovascular SOFA score
# of >=3 which corresponds to use of vasopressors)
lactate_itemids = [10053, 6837, 9580]
numerics_lactate = numerics_read(
    lactate_itemids, admissions_df=admissions_df,
    end=end_time, start=start_time)

###############################################################################
# Add admission time etc. to dataframes, and crop drug/list etc. tables to
# min/max times

# Convert to time since admission (in discrete 24hr periods)
admissions_add = admissions_df[['admissionid', 'admittedat', 'dischargedat']]
listitems = pd.merge(listitems, admissions_add, on='admissionid', how='left')
listitems['time'] = (listitems['measuredat'] - listitems['admittedat'])
listitems['time'] //= (1000*60*60*24)  # Convert to 'day'

drugitems = pd.merge(drugitems, admissions_add, on='admissionid', how='left')
drugitems['start_time'] = (drugitems['start'] - drugitems['admittedat'])
drugitems['stop_time'] = (drugitems['stop'] - drugitems['admittedat'])
drugitems[['start_time', 'stop_time']] //= (1000*60*60*24)  # Convert to 'day'

freetextitems = pd.merge(
    freetextitems, admissions_add, on='admissionid', how='left')
freetextitems['time'] = (
    freetextitems['measuredat'] - freetextitems['admittedat'])
freetextitems['time'] //= (1000*60*60*24)

procedureorderitems = pd.merge(
    procedureorderitems, admissions_add, on='admissionid', how='left')
procedureorderitems['time'] = (
    procedureorderitems['registeredat'] - procedureorderitems['admittedat'])
procedureorderitems['time'] //= (1000*60*60*24)

# Cut only to the time window specified earlier
if MIN_TIME is not None:
    listitems = listitems.loc[(listitems['time'] >= MIN_TIME)]
    drugitems = drugitems.loc[(drugitems['start_time'] >= MIN_TIME)]
    freetextitems = freetextitems.loc[(freetextitems['time'] >= MIN_TIME)]
    procedureorderitems = procedureorderitems.loc[
        (procedureorderitems['time'] >= MIN_TIME)]
if MAX_TIME is not None:
    listitems = listitems.loc[(listitems['time'] <= MAX_TIME)]
    # Keep this as start time rather than end time (i.e. started drug treatment
    # on that particular 'day')
    drugitems = drugitems.loc[(drugitems['start_time'] <= MAX_TIME)]
    freetextitems = freetextitems.loc[(freetextitems['time'] <= MAX_TIME)]
    procedureorderitems = procedureorderitems.loc[
        (procedureorderitems['time'] <= MAX_TIME)]


def systeem_flag(df):
    return df['registeredby'].str.contains('systeem', flags=re.IGNORECASE)


###############################################################################
# Identify noradrenaline administration to exclude (if required)
# We exclude any continuous administration of noradrenaline that lasts less
# than a specified length of time (e.g. 6hrs)

if inputs.exclude_noradrenaline is not None:
    noradrenaline = drugitems.loc[(drugitems['itemid'] == 7229)].reset_index()
    noradrenaline.sort_values(by=['admissionid', 'start'], inplace=True)
    # Changes in dose are new entries (rows), so we need to sum up the total
    # duration in time across dose changes.
    noradrenaline['continued'] = (
        noradrenaline['start'].diff(periods=-1) ==
        noradrenaline[['start', 'stop']].diff(periods=-1, axis=1)['start'])
    noradrenaline.loc[(
            noradrenaline['continued'].astype(int).diff(periods=1) == -1),
        'continued'] = True
    noradrenaline['start_n'] = (
        noradrenaline['continued'].astype(int).diff() == -1).cumsum()
    noradrenaline = pd.merge(
        noradrenaline,
        noradrenaline.groupby(['admissionid', 'start_n']).agg(
            total_duration=pd.NamedAgg(column='duration', aggfunc='sum')
            ).reset_index(),
        on=['admissionid', 'start_n'])
    # Default is not to exclude
    noradrenaline['exclude'] = False
    # If less than the specified duration, this switches to True
    noradrenaline.loc[(
            noradrenaline['total_duration'] < inputs.exclude_noradrenaline),
        'exclude'] = True
    noradrenaline.set_index('index', inplace=True)
    # Match the noradrenaline exclusions to their entries in the main
    # drugitems table and discard these entries
    drugitems['exclude'] = False
    drugitems.loc[
        noradrenaline.loc[noradrenaline['exclude']].index, 'exclude'] = True
    drugitems = drugitems.loc[(drugitems['exclude'] == 0)]
    drugitems.drop(columns='exclude', inplace=True)

###############################################################################
# Begin the SOFA scores (again following AmsterdamUMCdb!)

########################
# Respiration score

# Get PaO2/FiO2 ratio
oxy_flow_listitems = listitems.loc[
    (listitems['itemid'] == 8189),
    ['admissionid', 'valueid', 'value', 'measuredat', 'time']]
oxy_dev = numerics_sofa.loc[numerics_sofa['itemid'].isin([8845, 10387, 18587])]
oxy_flow = pd.merge(
    oxy_flow_listitems, oxy_dev[['admissionid', 'measuredat', 'value']],
    on=['admissionid', 'measuredat'], how='left')
oxy_flow.rename(columns={'value_x': 'O2_device'}, inplace=True)
oxy_flow.rename(columns={'value_y': 'O2_flow'}, inplace=True)
oxy_flow.head()

# Get PaO2 and FiO2 values
# Simultaneously retrieve PaCO2 and the 'nearest' FiO2 from the ventilator or
# estimated FiO2 based on applied oxygen device. Ideally documentation of
# measurements should be at the same time, but since this is not guaranteed
# allow a window.
# In more recent data PaCO2 and PaO2 were documented in kPa instead of mmHg.
fio2_itemid = [8845, 10387, 18587, 6699, 12279, 12369, 16246]
fio2_table = numerics_sofa.loc[
    (numerics_sofa['value'] > 0) &
    (numerics_sofa['itemid'].isin(fio2_itemid))]
fio2_table = pd.merge(
    fio2_table, oxy_flow_listitems.drop(columns=['time']),
    on=['admissionid', 'measuredat'], how='left')
fio2_table.rename(columns={'value_y': 'O2_device'}, inplace=True)
fio2_table.rename(columns={'value_x': 'value'}, inplace=True)
fio2_table['ventilatory_support'] = False
fio2_table.loc[(
        fio2_table['itemid'].isin([6699, 12279, 12369, 16246])),
    'ventilatory_support'] = True

fio2_table['fio2'] = 0.21
fio2_ind = (
    fio2_table['ventilatory_support'] & (fio2_table['value'].isnull() == 0))
fio2_table.loc[fio2_ind, 'fio2'] = fio2_table.loc[fio2_ind, 'value']

valueid1 = [1]  # Diep Nasaal
valueid1 += [2]  # Nasaal
valueid1 += [3]  # Kapje
valueid1 += [4]  # Kunstneus
valueid1 += [7]  # O2-bril
valueid1 += [8]  # Kinnebak
valueid1 += [9]  # Nebulizer
valueid1 += [18]  # Spreekcanule
valueid1 += [19]  # Spreekklepje
fio2_ind1 = (
    (fio2_table['ventilatory_support'] == 0) &
    (fio2_table['valueid'].isin(valueid1)))
fio2_table.loc[(
        fio2_ind1 & (fio2_table['value'] >= 1) & (fio2_table['value'] < 2)),
    'fio2'] = 0.22
fio2_table.loc[(
        fio2_ind1 & (fio2_table['value'] >= 2) & (fio2_table['value'] < 3)),
    'fio2'] = 0.25
fio2_table.loc[(
        fio2_ind1 & (fio2_table['value'] >= 3) & (fio2_table['value'] < 4)),
    'fio2'] = 0.27
fio2_table.loc[(
        fio2_ind1 & (fio2_table['value'] >= 4) & (fio2_table['value'] < 5)),
    'fio2'] = 0.30
fio2_table.loc[fio2_ind1 & (fio2_table['value'] >= 5), 'fio2'] = 0.35

valueid2 = [1]  # Diep Nasaal
valueid2 += [3]  # Kapje
valueid2 += [4]  # Kunstneus
valueid2 += [8]  # Kinnebak
valueid2 += [9]  # Nebulizer
valueid2 += [18]  # Spreekcanule
valueid2 += [19]  # Spreekklepje
fio2_ind2 = (
    (fio2_table['ventilatory_support'] == 0) &
    (fio2_table['valueid'].isin(valueid2)))
fio2_table.loc[(
        fio2_ind2 & (fio2_table['value'] >= 6) & (fio2_table['value'] < 7)),
    'fio2'] = 0.40
fio2_table.loc[(
        fio2_ind2 & (fio2_table['value'] >= 7) & (fio2_table['value'] < 8)),
    'fio2'] = 0.45
fio2_table.loc[fio2_ind2 & (fio2_table['value'] >= 8), 'fio2'] = 0.50

valueid3 = [10]  # Waterset
valueid3 += [11]  # Trach.stoma
valueid3 += [13]  # Ambu
valueid3 += [14]  # Guedel
valueid3 += [15]  # DL-tube
valueid3 += [16]  # CPAP
valueid3 += [17]  # Non-Rebreathing masker
fio2_ind3 = (
    (fio2_table['ventilatory_support'] == 0) &
    (fio2_table['valueid'].isin(valueid3)))
fio2_table.loc[(
        fio2_ind3 & (fio2_table['value'] >= 6) & (fio2_table['value'] < 7)),
    'fio2'] = 0.60
fio2_table.loc[(
        fio2_ind3 & (fio2_table['value'] >= 7) & (fio2_table['value'] < 8)),
    'fio2'] = 0.70
fio2_table.loc[(
        fio2_ind3 & (fio2_table['value'] >= 8) & (fio2_table['value'] < 9)),
    'fio2'] = 0.80
fio2_table.loc[(
        fio2_ind3 & (fio2_table['value'] >= 9) & (fio2_table['value'] < 10)),
    'fio2'] = 0.85
fio2_table.loc[fio2_ind3 & (fio2_table['value'] >= 10), 'fio2'] = 0.90
fio2_table.rename(columns={'measuredat': 'fio2_measuredat'}, inplace=True)
fio2_columns = ['admissionid', 'fio2_measuredat', 'fio2']
fio2_columns += ['ventilatory_support', 'time']

# This is initially pao2 and then merged with fio2 from above.
po2_itemid = [7433]  # PO2
po2_itemid += [9996]  # PO2 (bloed)
po2_itemid += [21214]  # PO2 (bloed) - kPa
oxygenation_po2 = numerics_sofa.loc[numerics_sofa['itemid'].isin(po2_itemid)]
# Conversion from kPa to mmHg
oxygenation_po2.loc[oxygenation_po2['unitid'] == 152, 'value'] *= 7.50061683
oxygenation_po2.rename(columns={'value': 'pao2'}, inplace=True)
oxygenation_po2['manual_entry'] = True
oxygenation_po2.loc[systeem_flag(oxygenation_po2), 'manual_entry'] = False

f = freetextitems.loc[freetextitems['itemid'] == 11646]
f = f[['admissionid', 'measuredat', 'value']]
f.rename(columns={'value': 'specimen_source'}, inplace=True)
oxygenation_po2_f = pd.merge(
    oxygenation_po2, f, on=['admissionid', 'measuredat'], how='left')

oxygenation_po2_f = oxygenation_po2_f.loc[
    oxygenation_po2_f['specimen_source'].isnull() |
    oxygenation_po2_f['specimen_source'].str.contains(
        'art', flags=re.IGNORECASE)]

oxygenation = pd.merge(
    oxygenation_po2_f, fio2_table[fio2_columns],
    on=['admissionid', 'time'], how='left')
oxygenation['FiO2_time_difference'] = (
    oxygenation['fio2_measuredat'] - oxygenation['measuredat'])
# Keep fio2 only if no earlier than 60 minutes before pao2 measurement
oxygenation = oxygenation.loc[(
    oxygenation['FiO2_time_difference'] > -1000*60*60)]
# and no later than 15 minutes after pao2 measuremen
oxygenation = oxygenation.loc[oxygenation['FiO2_time_difference'] < 1000*60*15]
# Convert to days (not discrete)
oxygenation['FiO2_time_difference'] /= (1000*60*60*24)
oxygenation['abs_measuredat'] = oxygenation['FiO2_time_difference'].abs()

# Sort by the smallest fio2 time difference for each patient and timestamp
oxygenation = oxygenation.sort_values(
    by=['admissionid', 'measuredat', 'abs_measuredat'])
# Discard duplicates of the same patient and timestamp (keeping only the
# smallest fio2 time difference)
oxygenation = oxygenation.loc[(
    (oxygenation[['admissionid', 'measuredat']].duplicated() == 0))]
oxygenation['priority'] = 1

sofa_resp_columns = ['admissionid', 'pao2', 'specimen_source', 'manual_entry']
sofa_resp_columns += ['time', 'fio2', 'ventilatory_support']
sofa_resp_columns += ['FiO2_time_difference', 'priority']
sofa_respiration = oxygenation[sofa_resp_columns]
sofa_respiration.head()

# Remove extreme outliers (in the AmsterdamUMCdb script, histograms are plotted
# to identify these outliers by eye, here we just copy those values)
sofa_respiration.loc[(sofa_respiration['fio2'] > 100), 'fio2'] = np.nan
sofa_respiration.loc[(
        (sofa_respiration['fio2'] < 20) &
        (sofa_respiration['fio2'] >= 1)),
    'fio2'] = np.nan
# Convert FiO2 in % to fraction
sofa_respiration.loc[(
        (sofa_respiration['fio2'] <= 100) &
        (sofa_respiration['fio2'] >= 20)),
    'fio2'] /= 100
# Remove lower outliers, most likely incorrectly labeled as 'arterial' instead
# of '(mixed/central) venous'
sofa_respiration.loc[sofa_respiration['pao2'] < 50, 'pao2'] = np.nan
sofa_respiration = sofa_respiration.dropna(subset=['pao2'])

# Get the PF ratio
sofa_respiration.loc[:, 'pf_ratio'] = (
    sofa_respiration['pao2'] / sofa_respiration['fio2'])
# Some entries may be 'None', need to make these False instead
sofa_respiration['ventilatory_support'].fillna(False, inplace=True)

# Calculate SOFA respiration score:
sofa_respiration['sofa_respiration_score'] = 0
sofa_respiration.loc[(
        (sofa_respiration['pf_ratio'] < 400) &
        (sofa_respiration['pf_ratio'] >= 300)),
    'sofa_respiration_score'] = 1
sofa_respiration.loc[(
        sofa_respiration['pf_ratio'] < 300),
    'sofa_respiration_score'] = 2
sofa_respiration.loc[(
        (sofa_respiration['pf_ratio'] < 200) &
        (sofa_respiration['pf_ratio'] >= 100) &
        sofa_respiration['ventilatory_support']),
    'sofa_respiration_score'] = 3
sofa_respiration.loc[(
        (sofa_respiration['pf_ratio'] < 100) &
        sofa_respiration['ventilatory_support']),
    'sofa_respiration_score'] = 4

sofa_respiration.head()

########################
# Coagulation score (platelets (thrombocytes))

platelet_itemid = [9964]  # Thrombo's (bloed)
platelet_itemid += [6797]  # Thrombocyten
platelet_itemid += [10409]  # Thrombo's citr. bloed (bloed)
platelet_itemid += [14252]  # Thrombo CD61 (bloed)
sofa_platelets = numerics_sofa.loc[
    numerics_sofa['itemid'].isin(platelet_itemid)]
sofa_platelets['manual_entry'] = True
sofa_platelets.loc[systeem_flag(sofa_platelets), 'manual_entry'] = False
sofa_platelets_columns = ['admissionid', 'itemid', 'item', 'value']
sofa_platelets_columns += ['registeredby', 'manual_entry', 'time']
sofa_platelets = sofa_platelets[sofa_platelets_columns]

# Calculate SOFA coagulation score:
sofa_platelets['sofa_coagulation_score'] = 0
sofa_platelets.loc[(
        (sofa_platelets['value'] < 150) & (sofa_platelets['value'] >= 100)),
    'sofa_coagulation_score'] = 1
sofa_platelets.loc[(
        (sofa_platelets['value'] < 100) & (sofa_platelets['value'] >= 50)),
    'sofa_coagulation_score'] = 2
sofa_platelets.loc[(
        (sofa_platelets['value'] < 50) & (sofa_platelets['value'] >= 20)),
    'sofa_coagulation_score'] = 3
sofa_platelets.loc[(
        sofa_platelets['value'] < 20),
    'sofa_coagulation_score'] = 4

sofa_platelets.head()

########################
# Liver score (bilirubin)

bilirubin_itemid = [6813]
bilirubin_itemid += [9945]
sofa_bilirubin = numerics_sofa.loc[(
    numerics_sofa['itemid'].isin(bilirubin_itemid))]
sofa_bilirubin['manual_entry'] = True
sofa_bilirubin.loc[systeem_flag(sofa_bilirubin), 'manual_entry'] = False
sofa_bilirubin_columns = ['admissionid', 'itemid', 'item', 'value']
sofa_bilirubin_columns += ['registeredby', 'manual_entry', 'time']
sofa_bilirubin = sofa_bilirubin[sofa_bilirubin_columns]

# Calculate SOFA liver score:
sofa_bilirubin['sofa_liver_score'] = 0
sofa_bilirubin.loc[(
        (sofa_bilirubin['value'] >= 20) & (sofa_bilirubin['value'] < 33)),
    'sofa_liver_score'] = 1
sofa_bilirubin.loc[(
        (sofa_bilirubin['value'] >= 33) & (sofa_bilirubin['value'] < 102)),
    'sofa_liver_score'] = 2
sofa_bilirubin.loc[(
        (sofa_bilirubin['value'] >= 102) & (sofa_bilirubin['value'] < 204)),
    'sofa_liver_score'] = 3
sofa_bilirubin.loc[(
        sofa_bilirubin['value'] >= 204),
    'sofa_liver_score'] = 4

sofa_bilirubin.head()

########################
# Cardiovascular score

cv_drug_itemid = [7179]  # Dopamine (Inotropin)
cv_drug_itemid += [7178]  # Dobutamine (Dobutrex)
cv_drug_itemid += [6818]  # Adrenaline (Epinefrine)
cv_drug_itemid += [7229]  # Noradrenaline (Norepinefrine)

sofa_cardiovascular = drugitems.loc[(
    drugitems['itemid'].isin(cv_drug_itemid) &
    (drugitems['rate'] >= 0.1) &  # changed from Amsterdam's script!
    (drugitems['ordercategoryid'] == 65))]

admissions_add = admissions_df[['admissionid', 'admittedat', 'weightgroup']]
sofa_cardiovascular = pd.merge(
    sofa_cardiovascular, admissions_add, on='admissionid', how='left')

# Replace weight group category with its numeric midpoint
weight_group_dict = {
    '59-': 55, '60-69': 65, '70-79': 75, '80-89': 85, '90-99': 95,
    '100-109': 105, '110+': 115, np.nan: 80}
sofa_cardiovascular['patientweight'] = (
    sofa_cardiovascular['weightgroup'].replace(weight_group_dict))
# If missing, then (crudely) impute with a midpoint of 80
sofa_cardiovascular.loc[(
    ~(sofa_cardiovascular['patientweight'] > 0)), 'patientweight'] = 80

# Want to add extra rows to the dataframe, for when drug administration
# happened over consecutive 'days' (i.e. make an entry for each 'day' that
# the drug administration window overlaps with)
n_days = sofa_cardiovascular['stop_time'] - sofa_cardiovascular['start_time']
n_days += 1
sofa_cardiovascular = sofa_cardiovascular.loc[
        sofa_cardiovascular.index.repeat(n_days)
    ].reset_index(drop=True)
sofa_cardiovascular['time'] = np.hstack([np.arange(x) for x in n_days])
sofa_cardiovascular['time'] += sofa_cardiovascular['start_time']
sofa_cardiovascular = sofa_cardiovascular.loc[
    (sofa_cardiovascular['time'] <= MAX_TIME)]
sofa_cardiovascular.drop(
    columns=['ordercategoryid', 'start', 'stop', 'weightgroup'], inplace=True)
sofa_cardiovascular.head()

# Calculate gamma, as per AmsterdamUMCdb script
# Gamma is ug/kg/min, so divide by patient weight
sofa_cardiovascular['gamma'] = (
    sofa_cardiovascular['dose'] / sofa_cardiovascular['patientweight'])
# Unless specified the dose already factors this
sofa_cardiovascular.loc[sofa_cardiovascular['doserateperkg'] == 1, 'gamma'] = (
    sofa_cardiovascular.loc[sofa_cardiovascular['doserateperkg'] == 1, 'dose'])
# Convert to ug/kg/min, depending on the doseunitid (10 is mg, 11 is ug)
sofa_cardiovascular.loc[(
        sofa_cardiovascular['doseunitid'] == 10),
    'gamma'] *= 1000
# Convert to ug/kg/min, depending on the doserateunitid (4 is min, 5 is hr)
sofa_cardiovascular.loc[(
        sofa_cardiovascular['doserateunitid'] == 5),
    'gamma'] /= 60

# Mean ABP
mean_abp = numerics_sofa.loc[numerics_sofa['itemid'].isin(numerics_abp_itemid)]
mean_abp['validated'] = True
mean_abp.loc[mean_abp['registeredby'].isnull(), 'validated'] = False
mean_abp.head()

# Remove extreme outliers, most likely data entry errors or measurement errors
mean_abp.loc[(mean_abp['value'] > 165), 'value'] = np.nan
mean_abp.loc[(mean_abp['value'] <= 30), 'value'] = np.nan
mean_abp_columns = ['admissionid', 'itemid', 'item', 'value']
mean_abp_columns += ['validated', 'time']
mean_abp = mean_abp[mean_abp_columns]
mean_abp = mean_abp.dropna(subset=['value'])
# Use mean_abp 'cleansed' dataframe
cv_groupby_columns = ['admissionid', 'itemid', 'item', 'time']
sofa_cardiovascular_map = mean_abp.groupby(cv_groupby_columns).agg(
        lowest_mean_abp=pd.NamedAgg(column='value', aggfunc='min')
    ).reset_index()

# Calculate SOFA cardiovascular score:
sofa_cardiovascular_map['sofa_cardiovascular_score'] = 0
# MAP < 70
sofa_cardiovascular_map.loc[(
        sofa_cardiovascular_map['lowest_mean_abp'] < 70),
    'sofa_cardiovascular_score'] = 1
sofa_cardiovascular_map.head()

sofa_cardiovascular_meds = sofa_cardiovascular.groupby(cv_groupby_columns).agg(
        total_duration=pd.NamedAgg(column='duration', aggfunc='sum'),
        max_gamma=pd.NamedAgg(column='gamma', aggfunc='max')
    ).reset_index()
sofa_cardiovascular_meds.head()

sofa_cardiovascular_meds['sofa_cardiovascular_score'] = 0
# Dopamine (itemid 7179) <= 5 or dobutamine (itemid 7178) any dose
sofa_cardiovascular_meds.loc[(
        ((sofa_cardiovascular_meds['itemid'] == 7179) &
            (sofa_cardiovascular_meds['max_gamma'] <= 5)) |
        (sofa_cardiovascular_meds['itemid'] == 7178)),
    'sofa_cardiovascular_score'] = 2
# Dopamine (itemid 7179) > 5, epinephrine (itemid 6818) <= 0.1,
# norepinephrine (itemid 7229) <= 0.1
sofa_cardiovascular_meds.loc[(
        ((sofa_cardiovascular_meds['itemid'] == 7179) &
            (sofa_cardiovascular_meds['max_gamma'] > 5) &
            (sofa_cardiovascular_meds['max_gamma'] < 15)) |
        ((sofa_cardiovascular_meds['itemid'] == 6818) &
            (sofa_cardiovascular_meds['max_gamma'] <= 0.1)) |
        ((sofa_cardiovascular_meds['itemid'] == 7229) &
            (sofa_cardiovascular_meds['max_gamma'] <= 0.1))),
    'sofa_cardiovascular_score'] = 3
# Dopamine (itemid 7179) > 15, epinephrine (itemid 6818) > 0.1,
# norepinephrine (itemid 7229) > 0.1
sofa_cardiovascular_meds.loc[(
        ((sofa_cardiovascular_meds['itemid'] == 7179) &
            (sofa_cardiovascular_meds['max_gamma'] > 15)) |
        ((sofa_cardiovascular_meds['itemid'] == 6818) &
            (sofa_cardiovascular_meds['max_gamma'] > 0.1)) |
        ((sofa_cardiovascular_meds['itemid'] == 7229) &
            (sofa_cardiovascular_meds['max_gamma'] > 0.1))),
    'sofa_cardiovascular_score'] = 4
sofa_cardiovascular_meds.head()

# Combine the scores from MAP and cardiovascular medication
# This concatenates the dataframes by row, so each admission will potentially
# have a score of 1 from MAP and a score of 2/3/4 from medications
# Later, we take the maximum by admissionid and time, retaining only one value
sofa_cardiovascular = pd.concat(
        [sofa_cardiovascular_map, sofa_cardiovascular_meds], sort=False,
    ).sort_values(by=['admissionid', 'time']).reset_index(drop=True)

sofa_cardiovascular.head()

########################
# Glasgow Coma Scale score

eyes_itemids = [6732]  # Actief openen van de ogen
eyes_itemids += [13077]  # A_Eye
eyes_itemids += [14470]  # RA_Eye
eyes_itemids += [16628]  # MCA_Eye
eyes_itemids += [19635]  # E_EMV_NICE_24uur
eyes_itemids += [19638]  # E_EMV_NICE_Opname
motor_itemids = [6734]  # Beste motore reactie van de armen
motor_itemids += [13072]  # A_Motoriek
motor_itemids += [14476]  # RA_Motoriek
motor_itemids += [16634]  # MCA_Motoriek
motor_itemids += [19636]  # M_EMV_NICE_24uur
motor_itemids += [19639]  # M_EMV_NICE_Opname
verbal_itemids = [6735]  # Beste verbale reactie
verbal_itemids += [13066]  # A_Verbal
verbal_itemids += [14482]  # RA_Verbal
verbal_itemids += [16640]  # MCA_Verbal
verbal_itemids += [19637]  # V_EMV_NICE_24uur
verbal_itemids += [19640]  # V_EMV_NICE_Opname

# GCS eyes component
gcs_components = listitems.loc[listitems['itemid'].isin(eyes_itemids)]
gcs_components['eyes_score'] = 1
# Actief openen van de ogen
gcs_components.loc[gcs_components['itemid'] == 6732, 'eyes_score'] = (
    5 - gcs_components.loc[gcs_components['itemid'] == 6732, 'valueid'])
# A_Eye
gcs_components.loc[gcs_components['itemid'] == 13077, 'eyes_score'] = (
    gcs_components.loc[gcs_components['itemid'] == 13077, 'valueid'])
# RA_Eye, MCA_Eye, E_EMV_NICE_24uur
gcs_components.loc[
        gcs_components['itemid'].isin([14470, 16628, 19635]), 'eyes_score'] = (
    gcs_components.loc[
        gcs_components['itemid'].isin([14470, 16628, 19635]), 'valueid'] - 4)
# E_EMV_NICE_Opname
gcs_components.loc[gcs_components['itemid'] == 19638, 'eyes_score'] = (
    gcs_components.loc[gcs_components['itemid'] == 19638, 'valueid'] - 8)

# Preference, ranked by discipline
gcs_components['preference'] = 8
gcs_preferences = {
    'ICV_Medisch Staflid': 1,
    'ICV_Medisch': 2,
    'ANES_Anesthesiologie': 3,
    'ICV_Physician assistant': 4,
    'ICH_Neurochirurgie': 5,
    'ICV_IC-Verpleegkundig': 6,
    'ICV_MC-Verpleegkundig': 7}
gcs_ind = gcs_components['registeredby'].isin(gcs_preferences.keys())
gcs_components.loc[gcs_ind, 'preference'] = (
    gcs_components.loc[gcs_ind, 'registeredby'].replace(gcs_preferences))
gcs_components.sort_values(
    by=['admissionid', 'time', 'preference', 'eyes_score'], inplace=True)
# Only keep the lowest score for the discipline of smallest rank
# (i.e. ICV_Medisch Staflid supercedes ICV_Medisch etc.)
# The limitation of this here is that we trust the 'preference' first, rather
# than the lowest score. If they happen to be more flat temporarily or the
# operator isn’t giving a proper stimulus, the patient will momentarily present
# a lower score than they actually have for that day
gcs_components = gcs_components.loc[(
    (gcs_components[['admissionid', 'time']].duplicated() == 0))]
gcs_components.drop(
    columns=['itemid', 'valueid', 'value', 'admittedat', 'registeredby'],
    inplace=True)

gcs_columns = ['admissionid', 'measuredat', 'itemid', 'valueid']
gcs_columns += ['registeredby']
# Add GCS motor score
gcs_components = pd.merge(
    gcs_components,
    listitems.loc[listitems['itemid'].isin(motor_itemids), gcs_columns],
    on=['admissionid', 'measuredat'],
    how='left')
gcs_components['motor_score'] = 1
# 6734 Beste motore reactie van de armen
gcs_components.loc[gcs_components['itemid'] == 6734, 'motor_score'] = (
    7 - gcs_components.loc[gcs_components['itemid'] == 6734, 'valueid'])
# 13072 A_Motoriek
gcs_components.loc[gcs_components['itemid'] == 13072, 'motor_score'] = (
    gcs_components.loc[gcs_components['itemid'] == 13072, 'valueid'])
m_itemid = [14476]  # RA_Motoriek
m_itemid += [16634]  # MCA_Motoriek
m_itemid += [19636]  # M_EMV_NICE_24uur
gcs_components.loc[gcs_components['itemid'].isin(m_itemid), 'motor_score'] = (
    gcs_components.loc[gcs_components['itemid'].isin(m_itemid), 'valueid'] - 6)
# 19639 M_EMV_NICE_Opname
gcs_components.loc[gcs_components['itemid'] == 19639, 'motor_score'] = (
    gcs_components.loc[gcs_components['itemid'] == 19639, 'valueid'] - 12)

# As above, add in preference by discipline
gcs_components['preference'] = 8
gcs_ind = gcs_components['registeredby'].isin(gcs_preferences.keys())
gcs_components.loc[gcs_ind, 'preference'] = (
    gcs_components.loc[gcs_ind, 'registeredby'].replace(gcs_preferences))
# Give higher preference by discipline (eye score should be the same for each
# admission id+time here)
gcs_components.sort_values(
    by=['admissionid', 'time', 'preference', 'motor_score'], inplace=True)
# Only keep the highest score for the discipline of smallest rank
gcs_components = gcs_components.loc[(
    (gcs_components[['admissionid', 'time']].duplicated() == 0))]
gcs_components.drop(
    columns=['itemid', 'valueid', 'registeredby'], inplace=True)
# Motor score is a float (due to pandas merge, so convert to int)
gcs_components['motor_score'] = gcs_components['motor_score'].astype(int)

# Add GCS verbal score
gcs_components = pd.merge(
    gcs_components,
    listitems.loc[listitems['itemid'].isin(verbal_itemids), gcs_columns],
    on=['admissionid', 'measuredat'],
    how='left')
gcs_components['verbal_score'] = 1
# 6735 Beste verbale reactie
gcs_components.loc[gcs_components['itemid'] == 6735, 'verbal_score'] = (
    6 - gcs_components.loc[gcs_components['itemid'] == 6735, 'valueid'])
# 13066 A_Verbal
gcs_components.loc[gcs_components['itemid'] == 13066, 'verbal_score'] = (
    gcs_components.loc[gcs_components['itemid'] == 13066, 'valueid'])
v_itemid = [14482]  # RA_Verbal
v_itemid += [16640]  # MCA_Verbal
gcs_components.loc[gcs_components['itemid'].isin(v_itemid), 'verbal_score'] = (
    gcs_components.loc[gcs_components['itemid'].isin(v_itemid), 'valueid'] - 5)
# 19637 V_EMV_NICE_24uur
gcs_components.loc[gcs_components['itemid'] == 19637, 'verbal_score'] = (
    gcs_components.loc[gcs_components['itemid'] == 19637, 'valueid'] - 9)
# 19640 V_EMV_NICE_Opname
gcs_components.loc[gcs_components['itemid'] == 19640, 'verbal_score'] = (
    gcs_components.loc[gcs_components['itemid'] == 19640, 'valueid'] - 15)

# As above, add in preference by discipline
gcs_components['preference'] = 8
gcs_ind = gcs_components['registeredby'].isin(gcs_preferences.keys())
gcs_components.loc[gcs_ind, 'preference'] = (
    gcs_components.loc[gcs_ind, 'registeredby'].replace(gcs_preferences))
# Give higher preference by discipline
gcs_components.sort_values(
    by=['admissionid', 'time', 'preference', 'verbal_score'], inplace=True)
# Only keep the highest score for the discipline of smallest rank
gcs_components = gcs_components.loc[(
    (gcs_components[['admissionid', 'time']].duplicated() == 0))]
gcs_components.drop(
    columns=['itemid', 'valueid', 'registeredby'], inplace=True)
# Verbal score minimum is 1
gcs_components.loc[gcs_components['verbal_score'] < 1, 'verbal_score'] = 1
gcs_components['verbal_score'] = gcs_components['verbal_score'].astype(int)

# Combine the scores for total GCS score
gcs_components['min_gcs'] = (
    gcs_components['eyes_score'] +
    gcs_components['motor_score'] +
    gcs_components['verbal_score'])
gcs_components.head()
sofa_cns = gcs_components[['admissionid', 'time', 'min_gcs']]

# Calculate SOFA cardiovascular score:
sofa_cns['sofa_cns_score'] = 0
# MAP < 70
sofa_cns.loc[(
        (sofa_cns['min_gcs'] >= 13) & (sofa_cns['min_gcs'] < 15)),
    'sofa_cns_score'] = 1
sofa_cns.loc[(
        (sofa_cns['min_gcs'] >= 10) & (sofa_cns['min_gcs'] < 13)),
    'sofa_cns_score'] = 2
sofa_cns.loc[(
        (sofa_cns['min_gcs'] >= 6) & (sofa_cns['min_gcs'] < 10)),
    'sofa_cns_score'] = 3
sofa_cns.loc[(sofa_cns['min_gcs'] < 6), 'sofa_cns_score'] = 4

sofa_cns.head()

########################
# Renal score

# Get urineoutput
sofa_urine_output_itemid = [8794]
sofa_urine_output_itemid += [8796]
sofa_urine_output_itemid += [8798]
sofa_urine_output_itemid += [8800]
sofa_urine_output_itemid += [8803]
sofa_urine_output_itemid += [10743]
sofa_urine_output_itemid += [10745]
sofa_urine_output_itemid += [19921]
sofa_urine_output_itemid += [19922]
# This df is called sofa_renal_urine_output in AmsterdamUMCdb sofa script
sofa_urine_output = numerics_sofa.loc[
    numerics_sofa['itemid'].isin(sofa_urine_output_itemid)]
sofa_urine_output.drop(
    columns=['unitid', 'measuredat', 'islabresult', 'fluidout', 'admittedat'],
    inplace=True)
sofa_urine_output.head()

# Probably decimal error when entering volumes > 2500
sofa_urine_output.loc[(sofa_urine_output['value'] > 2500), 'value'] /= 10
# Remove extreme outliers, most likely data entry error)
sofa_urine_output.loc[(sofa_urine_output['value'] > 4500), 'value'] = np.nan
sofa_urine_output = sofa_urine_output.dropna(subset=['value'])
# Get urine output per 24 hours
sofa_d_urine_output = sofa_urine_output.groupby(['admissionid', 'time']).agg(
        daily_urine_output=pd.NamedAgg(column='value', aggfunc='sum')
    ).reset_index()
sofa_d_urine_output.head()

# Calculate SOFA renal score for urine output:
sofa_d_urine_output['sofa_renal_score'] = 0
# Urine output < 500 ml/day
sofa_d_urine_output.loc[(
        (sofa_d_urine_output['daily_urine_output'] < 500) &
        (sofa_d_urine_output['daily_urine_output'] > 200)),
     'sofa_renal_score'] = 3
# Urine output < 200 ml/day
sofa_d_urine_output.loc[(
        (sofa_d_urine_output['daily_urine_output'] < 200)),
    'sofa_renal_score'] = 4
sofa_d_urine_output.head()

# Get serum creatinine (baseline from -365 days from admission)
baseline_creatinine = numerics_creatinine.groupby(['admissionid']).agg(
        baseline_creatinine=pd.NamedAgg(column='value', aggfunc='min')
    ).reset_index()
# Max creatinine on each day (but only from MIN_TIME rather than -365 days)
max_creatinine = numerics_creatinine.copy()
if MIN_TIME is not None:
    max_creatinine = max_creatinine.loc[max_creatinine['time'] >= MIN_TIME]
max_creatinine = max_creatinine.groupby(['admissionid', 'time']).agg(
        max_creatinine=pd.NamedAgg(column='value', aggfunc='max')
    ).reset_index()
# Merge baseline on admissionid only and max on both admissionid and time
creatinine = pd.merge(
    numerics_creatinine, baseline_creatinine, on='admissionid', how='left')
creatinine = pd.merge(
    creatinine, max_creatinine, on=['admissionid', 'time'], how='right')

creatinine['manual_entry'] = True
creatinine.loc[systeem_flag(creatinine), 'manual_entry'] = False

creatinine['acute_renal_failure'] = False
# AKI definition: 3 fold increase
creatinine.loc[(
        (creatinine['baseline_creatinine'] > 0) &
        (creatinine['max_creatinine'] /
            creatinine['baseline_creatinine'] > 3)),
    'acute_renal_failure'] = True
# AKI definition: increase to >= 354 umol/l AND at least 44 umol/l increase
creatinine.loc[(
        (creatinine['max_creatinine'] >= 354) &
        ((creatinine['max_creatinine'] -
            creatinine['baseline_creatinine']) >= 44)),
    'acute_renal_failure'] = True

creatinine.drop(
    columns=['unitid', 'measuredat', 'islabresult', 'fluidout', 'admittedat'],
    inplace=True)

creatinine.head()

# Looking at the data it's relevatively easy to spot most lab collection errors
# (i.e. single outliers between relatively normal values)
# Remove extreme outliers, most likely data entry errors (manual_entry = True)
creatinine.loc[(
        (creatinine['value'] < 30) & (creatinine['manual_entry'])),
    'value'] = np.nan
creatinine = creatinine.dropna(subset=['value'])

# Get highest creatinine per 24 hours
# Use creatinine 'cleansed' dataframe from APACHE score
sofa_renal_creatinine = creatinine.groupby(['admissionid', 'time']).agg(
        max_creatinine=pd.NamedAgg(column='value', aggfunc='max')
    ).reset_index()
sofa_renal_creatinine.head()
# Calculate SOFA renal score for creatinine:
sofa_renal_creatinine['sofa_renal_score'] = 0
# Creatinine 110-170 umol/l
sofa_renal_creatinine.loc[(
        (sofa_renal_creatinine['max_creatinine'] >= 110) &
        (sofa_renal_creatinine['max_creatinine'] < 171)),
     'sofa_renal_score'] = 1
# Creatinine 171-299 umol/l
sofa_renal_creatinine.loc[(
        (sofa_renal_creatinine['max_creatinine'] >= 171) &
        (sofa_renal_creatinine['max_creatinine'] < 300)),
    'sofa_renal_score'] = 2
# Creatinine 300-440 umol/l
sofa_renal_creatinine.loc[(
        (sofa_renal_creatinine['max_creatinine'] >= 300) &
        (sofa_renal_creatinine['max_creatinine'] <= 440)),
     'sofa_renal_score'] = 3
# Creatinine >440 umol/l
sofa_renal_creatinine.loc[(
        (sofa_renal_creatinine['max_creatinine'] > 440)),
    'sofa_renal_score'] = 4
sofa_renal_creatinine.head()

# Combine the scores from creatinine and urine output
sofa_renal = pd.concat(
        [sofa_renal_creatinine, sofa_d_urine_output], sort=False,
    ).sort_values(by=['admissionid', 'time'])

sofa_renal.head()

########################
# Final SOFA scores


def join_sofa(sofa, scores):
    '''
    Function to merge new sub-component scores to the full SOFA dataframe.
    We only want to keep columns where the 'day' (24hr periods post-admission)
    match, but some 'days' in the new SOFA component score may not be in the
    previous SOFA component scores (and we want to keep this)
    So for these rows, set all previous SOFA component scores to nan and then
    only keep matching 'days'
    '''
    scores = scores.sort_values(by=['admissionid', 'time']).reset_index()
    sofa = pd.concat([sofa, scores[['admissionid', 'time']]])
    sofa = sofa.loc[(sofa[['admissionid', 'time']].duplicated() == 0)]
    sofa = sofa.sort_values(by=['admissionid', 'time']).reset_index(drop=True)
    sofa = pd.merge(sofa, scores, on=['admissionid', 'time'], how='left')
    return sofa


# Merge the scores
sofa = admissions_df['admissionid']

# Max respiration score (don't need some steps in the function the first time)
scores = sofa_respiration.groupby(['admissionid', 'time']).agg(
    sofa_respiration_score=pd.NamedAgg(
        column='sofa_respiration_score', aggfunc='max'))
scores = scores.sort_values(by=['admissionid', 'time']).reset_index()
sofa = pd.merge(sofa, scores, on='admissionid', how='left')

# Max coagulation score
scores = sofa_platelets.groupby(['admissionid', 'time']).agg(
    sofa_coagulation_score=pd.NamedAgg(
        column='sofa_coagulation_score', aggfunc='max'))
sofa = join_sofa(sofa, scores)

# Max liver score
scores = sofa_bilirubin.groupby(['admissionid', 'time']).agg(
    sofa_liver_score=pd.NamedAgg(column='sofa_liver_score', aggfunc='max'))
sofa = join_sofa(sofa, scores)

# Max cardiovascular score
scores = sofa_cardiovascular.groupby(['admissionid', 'time']).agg(
    sofa_cardiovascular_score=pd.NamedAgg(
        column='sofa_cardiovascular_score', aggfunc='max'))
sofa = join_sofa(sofa, scores)

# Max central nervous system score
scores = sofa_cns.groupby(['admissionid', 'time']).agg(
    sofa_cns_score=pd.NamedAgg(column='sofa_cns_score', aggfunc='max'))
sofa = join_sofa(sofa, scores)

# Max renal score
scores = sofa_renal.groupby(['admissionid', 'time']).agg(
    sofa_renal_score=pd.NamedAgg(column='sofa_renal_score', aggfunc='max'))
sofa = join_sofa(sofa, scores)

# Calculate total score (add al values in columns)
total_scores = sofa.set_index(['admissionid', 'time']).sum(
    axis=1, skipna=True).to_frame('sofa_total_score')
sofa = pd.merge(sofa, total_scores, on=['admissionid', 'time'], how='left')
sofa.head()

sofa.dropna(subset=['time'], inplace=True)

# # save as .csv file
# sofa.to_csv(inputs.output_file_path + 'sofa.csv', index=False)

# SOFA scores (as per AmsterdamUMCdb) completed

###############################################################################
# Next for Sepsis-3 definition is antibiotics escalation

########################
# Antibiotics

antibiotics = pd.DataFrame(columns=['itemid', 'rank'])
antibiotics.loc[0] = [7185, 1]  # Doxycycline (Vibramycine)
antibiotics.loc[1] = [9142, 1]  # Tetracycline # Was 2 previously
antibiotics.loc[2] = [19764, 4]  # Tigecycline (Tygacil)
antibiotics.loc[3] = [9047, 2]  # Chlooramfenicol
antibiotics.loc[4] = [6847, 1]  # Amoxicilline (Clamoxyl/Flemoxin)
antibiotics.loc[5] = [9128, 3]  # Piperacilline (Pipcil) # Note, no tazobactam
antibiotics.loc[6] = [6871, 1]  # Benzylpenicilline (Penicilline)
antibiotics.loc[7] = [9037, 1]  # Feneticilline (Broxil)
antibiotics.loc[8] = [9029, 2]  # Amoxicilline/Clavulaanzuur (Augmentin)
antibiotics.loc[9] = [7123, 4]  # Imipenem (Tienam)
antibiotics.loc[10] = [9151, 2]  # Cefuroxim (Zinacef) # Was 1 previously
antibiotics.loc[11] = [6917, 3]  # Ceftazidim (Fortum) # Was 2
antibiotics.loc[12] = [9133, 2]  # Ceftriaxon (Rocephin)
antibiotics.loc[13] = [9030, 4]  # Aztreonam (Azactam)
antibiotics.loc[14] = [8127, 4]  # Meropenem (Meronem) # Shah have as 4
antibiotics.loc[15] = [7208, 2]  # Erythromycine (Erythrocine)
antibiotics.loc[16] = [8546, 2]  # Claritromycine (Klacid)
antibiotics.loc[17] = [13057, 2]  # Azitromycine (Zithromax)
antibiotics.loc[18] = [6958, 2]  # Clindamycine (Dalacin)
antibiotics.loc[19] = [7044, 1]  # Tobramycine (Obracin) # Changed from 2
antibiotics.loc[20] = [7235, 3]  # Gentamicine (Garamycin) # Changed from 2
antibiotics.loc[21] = [9109, 4]  # Neomycine sulfaat # Chris suggests 4
antibiotics.loc[22] = [6834, 4]  # Amikacine (Amukin) # Changed from 2
antibiotics.loc[23] = [6948, 2]  # Ciprofloxacine (Ciproxin)
antibiotics.loc[24] = [9117, 2]  # Norfloxacine (Noroxin)
antibiotics.loc[25] = [12398, 2]  # Levofloxacine (Tavanic)
antibiotics.loc[26] = [12956, 2]  # Moxifloxacin (Avelox)
antibiotics.loc[27] = [7064, 3]  # Vancomycine # rank 1 if not IV!
antibiotics.loc[28] = [8549, 4]  # Belcomycine (Colistinesulfaat) 4 x dgs
antibiotics.loc[29] = [10584, 4]  # Belcomycine (Colistinesulfaat) 6 x dgs
antibiotics.loc[30] = [20175, 4]  # Colistine
antibiotics.loc[31] = [20176, 4]  # Colistine Inhalatie
antibiotics.loc[32] = [8942, 1]  # Metronidazol-Flagyl
antibiotics.loc[33] = [7187, 1]  # Metronidazol (Flagyl)
antibiotics.loc[34] = [14236, 1]  # Nitrofurantoïne ( Furadantine) # Was 1
antibiotics.loc[35] = [19137, 4]  # Linezolid (Zyvoxid)
antibiotics.loc[36] = [19773, 4]  # Daptomycine (Cubicin)
antibiotics.loc[37] = [8394, 1]  # Co-Trimoxazol (Bactrimel)
antibiotics.loc[38] = [9052, 1]  # Co-trimoxazol forte (Bactrimel)
antibiotics.loc[39] = [6919, 2]  # Cefotaxim (Claforan)
# The following are not included (various reasons, documented):
# 7231 Fluconazol (Diflucan): an antifungal/antiviral
# 9152 Cefazoline (Kefzol): rarely used as infection treatment in Amsterdam
# 9070 Zilversulfadiazine (Flammazine): topical
# 6932 Chlooramfenicol (Globenicol): cream/ointment
# 9070 Flucloxacilline (Stafoxil/Floxapen): cream/ointment
# 13102: Dexamethason/gentamicine oogzalf (Dexamytrex): eye drops
# 13094: Tobramycine oogzalf (Tobrex): eye drops
# 12997 Ofloxacine (Trafloxal) oogdruppels:eye drops
# 9075 Fusidinezuur (Fucidin): prophylactic after cardiothoracic surgery
# 13045 Fusidinezuur oogdruppels (Fusithalmic): eye drops

# As part of the selective digestive decontamination, patients expected to stay
# at least 24-48hrs in ICU receive 16 (at least between 10 and 20) doses of
# cefotaxime across 4 days. If the clinician suspects an infection, this should
# be switched to ceftriaxone (which has a similar spectrum). If cefotaxime is
# continued after these initial doses, assume the clinician has suspected an
# infection and kept cefotaxime.
# There is perhaps no straightforward answer of how to convert this into code
# here. However, if this is the case, we discard any administration of
# cefotaxime that ends before day 4. If a course of cefotaxime is continued
# beyond day 4, we include it in our antibiotic escalation as if it were
# started on day 4 itself (not necessarily at the point the clinician suspected
# an infection). E.g. if cefotaxime is given until day 6, we consider any
# escalation that involves this to happen at day 4.
# Similarly, for cardiothoracic surgery patients, antibiotic usage will be
# prophylactic (vancomycin/fusidinezuur and cefazolin).
# Fusidinezuur is always prophylactic, vancomycin prophylactic in day 0.
# Prophylactic antibiotic usage is picked up again later in the script.

drugitems_abx = drugitems.loc[drugitems['itemid'].isin(antibiotics['itemid'])]

# Sepsis-3 (Shah et al.) say antibiotics considered during 24hr period if
# administration occurred within 24hr period or within 12hrs before 24hr period
# So an antibiotic at time 18hr on day 0 would also be considered as
# relevant to day 1. We handle this by extending the stop time by 12hr.
drugitems_abx['start_time'] = (
    drugitems_abx['start'] - drugitems_abx['admittedat'])
drugitems_abx.loc[(
        drugitems_abx['start_time'] <= 1000*60*60*24*MIN_TIME),
    'start_time'] = 1000*60*60*24*MIN_TIME
drugitems_abx['start_time'] //= (1000*60*60*24)
# The same as before
drugitems_abx['stop_time'] = (
    drugitems_abx['stop'] - drugitems_abx['admittedat'] + 1000*60*60*12)
drugitems_abx['stop_time'] //= (1000*60*60*24)
# As with the SOFA cardiovascular score, we want to add extra rows to the
# dataframe, for when drug administration happened over consecutive 'days'
# (i.e. make an entry for each 'day' tha the drug administration window
# overlaps with)
n_days = drugitems_abx['stop_time'] - drugitems_abx['start_time'] + 1
drugitems_abx = drugitems_abx.loc[
    drugitems_abx.index.repeat(n_days)].reset_index(drop=True)
drugitems_abx['time'] = np.hstack([np.arange(x) for x in n_days])
drugitems_abx['time'] += drugitems_abx['start_time']
drugitems_abx = drugitems_abx.loc[drugitems_abx['time'] <= MAX_TIME]

drugitems_abx['prophylactic'] = False

# Cefotaxim is always prophylactic (SDD)
drugitems_abx.loc[(
        (drugitems_abx['time'] < 4) &  # changed from stop_time
        (drugitems_abx['itemid'] == 6919)),
    'prophylactic'] = True

# Vancomycine is prophylactic after cardiac surgery
cardiosurg_admissions = combined_diagnoses.loc[(
        combined_diagnoses['cardiosurg'] == 1), 'admissionid']
drugitems_abx.loc[(
        (drugitems_abx['admissionid'].isin(cardiosurg_admissions)) &
        # (drugitems_abx['start_time'] < 1) &  # changed from stop_time
        (drugitems_abx['itemid'] == 7064)),
    'prophylactic'] = True

# Erythromycine in low doses (250mg) is only used for gastric retention in
# AmsterdamUMCdb, not for treating infection.
drugitems_abx.loc[(
        (drugitems_abx['itemid'] == 7208) &
        (drugitems_abx['doseunitid'] == 10) &
        (drugitems_abx['dose'] <= 250)),
    'prophylactic'] = True
drugitems_abx.loc[(
        (drugitems_abx['itemid'] == 7208) &
        (drugitems_abx['doseunitid'] == 9) &
        (drugitems_abx['dose'] <= 0.25)),
    'prophylactic'] = True

# We also assume antibiotic usage in the first 24hr after elective surgery was
# prophylactic these may have been given, and discard all antibiotics for
# such patients for the first 24hr after ICU admission. This should include
# almost all of the above (vancomycin after cardiothoracic surgery), unless
# that surgery was emergency
surgical_columns = ['admissionid', 'surgical', 'urgency']
drugitems_abx = pd.merge(
    drugitems_abx, combined_diagnoses[surgical_columns],
    on='admissionid', how='left')
drugitems_abx['elective_surgery'] = (
    (drugitems_abx['surgical'] == 1) &
    (drugitems_abx['urgency'] == 0))
drugitems_abx['emergency_surgery'] = (
    (drugitems_abx['surgical'] == 1) &
    (drugitems_abx['urgency'] == 1))

drugitems_abx.loc[(
        (drugitems_abx['elective_surgery'] == 1) &
        (drugitems_abx['start_time'] < 1)),
    'prophylactic'] = True

drugitems_abx = pd.merge(
    drugitems_abx, antibiotics[['itemid', 'rank']], on='itemid', how='left')

drugitems_abx['intravenous'] = (
    drugitems_abx['ordercategoryid'].isin([15, 65, 55]))

drugitems_abx_all = drugitems_abx.copy()
drugitems_abx = drugitems_abx.loc[(drugitems_abx['prophylactic'] == 0)]

abx_columns = ['admissionid', 'time', 'intravenous', 'rank', 'item', 'itemid']
drugitems_abx = drugitems_abx[abx_columns]
drugitems_abx = drugitems_abx.loc[(drugitems_abx.duplicated() == 0)]

max_rank = drugitems_abx.groupby(['admissionid', 'time']).agg(
        max_rank=pd.NamedAgg(column='rank', aggfunc='max')
    ).reset_index()
drugitems_abx = pd.merge(
    drugitems_abx, max_rank, on=['admissionid', 'time'], how='left')
drugitems_abx['rank_diff'] = (
    drugitems_abx['max_rank'] - drugitems_abx['rank'])

# Now compute the antibiotic escalation from day to day. Antibiotic escalation
# occurs if a new drug of higher ranking is administered (i.e. the max_rank
# increases) or if the number of drugs at the highest ranking is increased
# (i.e. n_max_rank increases), when this is accompanied by at least one
# antibiotics given intravenously
abx_escalation = drugitems_abx.groupby(['admissionid', 'time']).agg(
        intravenous=pd.NamedAgg(column='intravenous', aggfunc='any'),
        max_rank=pd.NamedAgg(column='rank', aggfunc='max'),
        n_max_rank=pd.NamedAgg(
            column='rank_diff', aggfunc=lambda x: (x == 0).sum())
    ).reset_index()

abx_escalation['max_rank_increase'] = (abx_escalation['max_rank'].diff() > 0)
abx_escalation['n_max_rank_increase'] = (
    abx_escalation['n_max_rank'].diff() > 0)
abx_escalation['time_diff'] = abx_escalation['time'].diff()
# If two consecutive rows correspond to different patients, then we don't want
# antibiotics escalation defined from the first patient to the second!
abx_escalation.loc[(
        (abx_escalation['admissionid'].duplicated() == 0)),
    ['max_rank_increase', 'n_max_rank_increase', 'time_diff']] = np.nan
abx_escalation['antibiotic_escalation'] = False
# Any new entry is assumed to be the first antibiotics given to that patient
abx_escalation.loc[(
        (abx_escalation['admissionid'].duplicated() == 0)),
    'antibiotic_escalation'] = True
abx_escalation.loc[(
        (abx_escalation['max_rank_increase'] > 0) |
        (abx_escalation['n_max_rank_increase'] > 0) |
        (abx_escalation['time_diff'] != 1)),
    'antibiotic_escalation'] = True

abx_escalation.drop(
    columns=['max_rank_increase', 'n_max_rank_increase', 'time_diff'],
    inplace=True)

abx_escalation['antibiotics_any'] = True

###############################################################################
# Now we are in a position to compute sepsis episodes according to the Sepsis-3
# definition of increase in SOFA score accompanied by antibiotics escalation

# We need max lactate values for each day of each admission, for defining
# septic shock
lactate_max = numerics_lactate.groupby(['admissionid', 'time']).agg(
    max_lactate=pd.NamedAgg(column='value', aggfunc='max')).reset_index()

# Similar to SOFA scores, we want to include antibiotic escalation where no
# SOFA scores exist yet, so make placeholder rows for these, then remove
# duplicates
sepsis = pd.concat([sofa, abx_escalation[['admissionid', 'time']]])
sepsis = sepsis.loc[(sepsis[['admissionid', 'time']].duplicated() == 0)]
sepsis = sepsis.sort_values(by=['admissionid', 'time']).reset_index(drop=True)

# Retain only days within the input-specified range, and expand entries to
# between the minimum and maximum for each admission (within the range)
sepsis_add = sepsis.loc[
    (sepsis['time'] <= MAX_TIME) & (sepsis['time'] >= MIN_TIME)]
n_time = sepsis_add.groupby('admissionid').agg(
    n_time=pd.NamedAgg(column='time', aggfunc=lambda x: x.max() - x.min() + 1),
    min_time=pd.NamedAgg(column='time', aggfunc='min'))
sepsis = n_time.reset_index()
sepsis = sepsis.loc[
    sepsis.index.repeat(sepsis['n_time'])].reset_index(drop=True)
sepsis['time'] = np.hstack([np.arange(x) for x in n_time['n_time']])
sepsis['time'] = (sepsis['time'] + sepsis['min_time'])
sepsis = sepsis[['admissionid', 'time']]
sepsis = pd.merge(sepsis, sepsis_add, on=['admissionid', 'time'], how='left')

# Shah et al. (DOI: 10.1097/CCM.0000000000005169) in a Sepsis-3 review remove
# any patients with less than 3 SOFA scores in the first 24hrs, we mark these
# patients in case we want to do the same
sofa_components = ['sofa_respiration_score', 'sofa_coagulation_score']
sofa_components += ['sofa_liver_score', 'sofa_cardiovascular_score']
sofa_components += ['sofa_cns_score', 'sofa_renal_score']
sepsis['n_sofa_scores'] = sepsis[sofa_components].notna().sum(axis=1)

# If any pre-ICU data is used, and there is no SOFA score, then we impute as 0
sepsis.loc[(
        (sepsis['sofa_total_score'].isna() & (sepsis['time'] < 0))),
    'sofa_total_score'] = 0

# For patients who experience ICU death (defined as death within 24hr of
# discharge or if explicitly recorded as on ICU), if they are missing SOFA
# scores on the day of death, then impute with 24 (the maximum value)
combined_diagnoses['overleden_IC'] = (
    (combined_diagnoses['destination'] == 'Overleden') |
    (combined_diagnoses['dateofdeath'] -
        combined_diagnoses['dischargedat'] < (1000*60*60*24)))

combined_diagnoses['overleden_time'] = np.nan
combined_diagnoses.loc[(
        (combined_diagnoses['overleden_IC'] == 1)), 'overleden_time'] = (
    combined_diagnoses['lengthofstay'] // 24)
overleden_add = combined_diagnoses[['admissionid', 'overleden_time']]
overleden_add = overleden_add.dropna(subset=['overleden_time'])
overleden_add['overleden_IC'] = 1
overleden_add.rename(columns={'overleden_time': 'time'}, inplace=True)

sepsis = pd.merge(
    sepsis, overleden_add, on=['admissionid', 'time'], how='left')

sepsis.loc[(
        (sepsis['n_sofa_scores'] < 3) &
        (sepsis['overleden_IC'].notna())),
    'sofa_total_score'] = 24

# We will discard admissions for whom less than 3 SOFA scores were computed
# at admission
sepsis['discard'] = False
sepsis.loc[(
        (sepsis['n_sofa_scores'] < 3) &
        (sepsis['time'] == 0)),
    'discard'] = True
discard_ind = sepsis['admissionid'].isin(
    sepsis.loc[sepsis['discard'], 'admissionid'])
sepsis.loc[discard_ind, 'discard'] = True

# sepsis.loc[sepsis['n_sofa_scores'] < 3, 'sofa_total_score'] = np.nan

# Merge sofa scores with antibiotic escalation ressults
sepsis = pd.merge(
    sepsis, abx_escalation, on=['admissionid', 'time'], how='left')
sepsis['intravenous'] = sepsis['intravenous'].fillna(False)
sepsis['antibiotic_escalation'] = sepsis['antibiotic_escalation'].fillna(False)
sepsis['antibiotics_any'] = sepsis['antibiotics_any'].fillna(False)

# SOFA increase corresponds to three time periods: previous and current day,
# current and subsequent day, previous and subsequent day, where 'current day'
# here is the day of antibiotic escalation (see Shah et al. for more details)
# Change in SOFA between previous and current 24hr periods
sepsis['sofa_diff0'] = sepsis['sofa_total_score'].diff()
# If this is the first entry for the patient, then there is no 'previous day'
# So the difference in these cases is simply the (total score - 0). We assume a
# SOFA score of 0 on admission)
sepsis.loc[(sepsis['admissionid'].duplicated() == 0), 'sofa_diff0'] = (
    sepsis.loc[(sepsis['admissionid'].duplicated() == 0), 'sofa_total_score'])
# Change in SOFA between current and subsequent 24hr periods
sepsis['sofa_diff1'] = sepsis['sofa_diff0'].shift(-1)
# Change in SOFA between previous and subsequent 24hr periods
sepsis['sofa_diff2'] = sepsis['sofa_total_score'].diff(periods=2).shift(-1)

# The above is agnostic to admissionid. On the last day for that patient,
# sofa_diff1 and sofa_diff2 should not exist, so reset to NaN
sepsis.loc[(
        sepsis['admissionid'].duplicated(keep='last') == 0),
    ['sofa_diff1', 'sofa_diff2']] = np.nan
# If the number of SOFA scores on the 'current' day is less than 3, then
# the SOFA increase from that day to the next can't happen
sepsis.loc[(sepsis['n_sofa_scores'] < 3), 'sofa_diff1'] = np.nan
# Similarly if the time difference is >1 day, we will disregard those (this
# shouldn't be the case however)
sepsis.loc[(
        (sepsis['time'].diff() > 1) &
        (sepsis['admissionid'].duplicated() == 1)),
    ['sofa_diff0', 'sofa_diff1', 'sofa_diff2']] = np.nan

# A sepsis episode is defined as antibiotics escalation ('infection')
# accompanied by SOFA increase of 2 or more.
sepsis['sepsis_episode'] = (
    (sepsis['antibiotic_escalation'] == 1) &
    ((sepsis['sofa_diff0'] >= 2) |
        (sepsis['sofa_diff1'] >= 2) |
        (sepsis['sofa_diff2'] >= 2)))

# If cefotaxime is continued after 4 days, we assume that there has been
# a suspected infection sometime between day 0 and day 3, and that cefotaxime
# is deliberately continued instead of a switch to ceftriaxone. We therefore
# search for possible sepsis episodes before day 4 in this set of patients,
# with a stricter set of conditions: SOFA score greater than that at admission,
# SOFA increase of >=2 on consecutive days (rather than e.g. on previous
# and next day)
cefotaxime_admissionids = drugitems_abx.loc[
    (drugitems_abx['itemid'] == 6919), 'admissionid']
cefotaxime_sepsis = sepsis.loc[(
        (sepsis['admissionid'].isin(cefotaxime_admissionids)) &
        (sepsis['time'] <= 4)
    )].groupby('admissionid').agg(
        sepsis_in_sdd=pd.NamedAgg(column='sepsis_episode', aggfunc='any')
    ).reset_index()
cefotaxime_admissionids = cefotaxime_sepsis.loc[
    (cefotaxime_sepsis['sepsis_in_sdd'] == 0), 'admissionid']
# Reset index to retain the original index within the sepsis dataframe
cefotaxime_sepsis = sepsis.loc[(
    (sepsis['admissionid'].isin(cefotaxime_admissionids)) &
    (sepsis['time'] <= 4))].reset_index()
cefotaxime_sofa_admission = cefotaxime_sepsis.loc[(
    (cefotaxime_sepsis['time'] == 0)), ['admissionid', 'sofa_total_score']]
cefotaxime_sofa_admission.rename(
    columns={'sofa_total_score': 'sofa_admission'}, inplace=True)
cefotaxime_sepsis = pd.merge(
    cefotaxime_sepsis, cefotaxime_sofa_admission, on='admissionid', how='left')
cefotaxime_sepsis['sepsis_episode'] = (
    (cefotaxime_sepsis['sofa_total_score'] >
        cefotaxime_sepsis['sofa_admission']) &
    (cefotaxime_sepsis['sofa_total_score'].diff() > 2) &
    (cefotaxime_sepsis['admissionid'].duplicated() == 1))
cefotaxime_sepsis = cefotaxime_sepsis.loc[(
    (cefotaxime_sepsis['sepsis_episode'] == 1))]
# Match to the sepsis dataframe using the retained index
sepsis.loc[cefotaxime_sepsis['index'], 'sepsis_episode'] = True

# We set a cutoff so that potential new sepsis episodes can only occur >72hr
# after a previous sepsis episode (otherwise it is considered part of the
# same sepsis episode)
new_sepsis = sepsis.loc[sepsis['sepsis_episode']]
new_sepsis['time_diff'] = new_sepsis['time'].diff()
new_sepsis.loc[(
    (new_sepsis['admissionid'].duplicated() == 0)), 'time_diff'] = np.nan
new_sepsis['new_sepsis_episode'] = (
    (new_sepsis['admissionid'].duplicated() == 0) |
    (new_sepsis['time_diff'] >= 3))
new_sepsis_add = new_sepsis[['admissionid', 'time', 'new_sepsis_episode']]
sepsis = pd.merge(
    sepsis, new_sepsis_add, on=['admissionid', 'time'], how='left')
sepsis['new_sepsis_episode'] = sepsis['new_sepsis_episode'].fillna(False)

sepsis = pd.merge(
    sepsis, lactate_max[['admissionid', 'time', 'max_lactate']],
    on=['admissionid', 'time'], how='left')

sepsis['septic_shock'] = (
    (sepsis['sepsis_episode'] == 1) &
    (sepsis['max_lactate'] > 2) &
    (sepsis['sofa_cardiovascular_score'] >= 2))  # Should be 2 not 3!!

# We set a cutoff so that potential new sepsis episodes can only occur >72hr
# after a previous sepsis episode (otherwise it is considered part of the
# same sepsis episode)
new_septic_shock = sepsis.loc[sepsis['septic_shock']]
new_septic_shock['time_diff'] = new_septic_shock['time'].diff()
new_septic_shock.loc[(
    (new_septic_shock['admissionid'].duplicated() == 0)), 'time_diff'] = np.nan
new_septic_shock['new_septic_shock'] = (
    (new_septic_shock['admissionid'].duplicated() == 0) |
    (new_septic_shock['time_diff'] >= 3))
new_septic_shock_add = new_septic_shock[
    ['admissionid', 'time', 'new_septic_shock']]
sepsis = pd.merge(
    sepsis, new_septic_shock_add, on=['admissionid', 'time'], how='left')
sepsis['new_septic_shock'] = sepsis['new_septic_shock'].fillna(False)

# We also want to identify a step down (i.e. from septic shock to sepsis
# without shock, from sepsis without shock to antibiotics without sepsis)
sepsis['sepsis_baseline_sofa'] = np.nan
baseline_ind0 = (sepsis['new_sepsis_episode'] & (sepsis['sofa_diff0'] >= 2))
sepsis.loc[baseline_ind0, 'sepsis_baseline_sofa'] = (
    sepsis['sofa_total_score'] - sepsis['sofa_diff0']).loc[baseline_ind0]
baseline_ind1 = (
    sepsis['new_sepsis_episode'] &
    sepsis['sepsis_baseline_sofa'].isna() &
    (sepsis['sofa_diff1'] > sepsis['sofa_diff2']))
sepsis.loc[baseline_ind1, 'sepsis_baseline_sofa'] = (
    sepsis.loc[baseline_ind1, 'sofa_total_score'])
baseline_ind2 = (
    sepsis['new_sepsis_episode'] &
    sepsis['sepsis_baseline_sofa'].isna() &
    (sepsis['sofa_diff2'] > sepsis['sofa_diff1']))
sepsis.loc[baseline_ind2, 'sepsis_baseline_sofa'] = (
    sepsis['sofa_total_score'].shift(-1).loc[baseline_ind2])

# Replaces NaNs with the last non-null value
sepsis['sepsis_baseline_sofa'] = sepsis['sepsis_baseline_sofa'].ffill()
# Fills the NaN values before the first non-null value to 0
sepsis['sepsis_baseline_sofa'] = sepsis['sepsis_baseline_sofa'].fillna(0)

first_sepsis = sepsis.loc[(
        (sepsis['new_sepsis_episode'] == 1))].groupby('admissionid').agg(
    time=pd.NamedAgg(column='time', aggfunc='min')).reset_index()
first_sepsis['had_sepsis'] = True

sepsis = pd.merge(sepsis, first_sepsis, on=['admissionid', 'time'], how='left')
sepsis.loc[(
        (sepsis['admissionid'].duplicated() == 0) &
        (sepsis['had_sepsis'].isna())),
    'had_sepsis'] = False
# ffill replaces NaNs with the last non-null value
sepsis['had_sepsis'] = sepsis['had_sepsis'].ffill()

sepsis_episode_number = sepsis.groupby('admissionid').agg(
    n_sepsis_episodes=pd.NamedAgg(column='new_sepsis_episode', aggfunc='sum'))
sepsis_episode_number = np.hstack(
    [np.arange(x) for x in sepsis_episode_number.values])
sepsis['sepsis_episode_number'] = np.nan
sepsis.loc[(
        (sepsis['new_sepsis_episode'] == 1)),
    'sepsis_episode_number'] = sepsis_episode_number

sepsis.loc[(sepsis['had_sepsis'] == 1), 'sepsis_episode_number'] = sepsis.loc[
    (sepsis['had_sepsis'] == 1), 'sepsis_episode_number'].ffill()

sepsis['possible_continued_sepsis'] = (
    (sepsis['sofa_total_score'] > sepsis['sepsis_baseline_sofa']) &
    (sepsis['antibiotics_any'] == 1) &
    (sepsis['had_sepsis'] == 1))

sepsis.loc[(
    (sepsis['new_sepsis_episode'] == 1)), 'possible_continued_sepsis'] = True

end_sepsis = sepsis[(
        (sepsis['possible_continued_sepsis'] == 1))].groupby(
            ['admissionid', 'sepsis_episode_number']).agg(
    time=pd.NamedAgg(column='time', aggfunc='max')).reset_index()
end_sepsis.drop(columns=['sepsis_episode_number'], inplace=True)
end_sepsis['continued_sepsis_episode'] = True
sepsis = pd.merge(sepsis, end_sepsis, on=['admissionid', 'time'], how='left')

# Turn the entry immediately following a True value in 'end_sepsis_episode'
# from NaN to False.
sepsis.loc[(
        (sepsis['continued_sepsis_episode'].isna() == 1) &
        (sepsis['continued_sepsis_episode'] == 1).shift(1).fillna(False) == 1),
    'continued_sepsis_episode'] = False

sepsis.loc[(
        (sepsis['time'] == 0) & (sepsis['new_sepsis_episode'] == 0)),
    'continued_sepsis_episode'] = False

sepsis.loc[(
    (sepsis['new_sepsis_episode'] == 1)), 'continued_sepsis_episode'] = True
sepsis['continued_sepsis_episode'] = sepsis['continued_sepsis_episode'].ffill()
sepsis['continued_sepsis_episode'] = (
    sepsis['continued_sepsis_episode'].fillna(False))

sepsis.drop(columns=['possible_continued_sepsis', 'had_sepsis'], inplace=True)

# Step down from septic shock to sepsis without shock if
# (lactate<2) AND (SOFA cardiovascular <= 3)
first_septic_shock = sepsis.loc[(
        (sepsis['septic_shock'] == 1))].groupby('admissionid').agg(
    time=pd.NamedAgg(column='time', aggfunc='min')).reset_index()
first_septic_shock['had_septic_shock'] = True

sepsis = pd.merge(
    sepsis, first_septic_shock, on=['admissionid', 'time'], how='left')
sepsis.loc[(
        (sepsis['admissionid'].duplicated() == 0) &
        (sepsis['had_septic_shock'].isna())),
    'had_septic_shock'] = False
# ffill replaces NaNs with the last non-null value
sepsis['had_septic_shock'] = sepsis['had_septic_shock'].ffill()

septic_shock_episode_number = sepsis.groupby('admissionid').agg(
    n_septic_shock_episodes=pd.NamedAgg(column='septic_shock', aggfunc='sum'))
septic_shock_episode_number = np.hstack(
    [np.arange(x) for x in septic_shock_episode_number.values])
sepsis['septic_shock_episode_number'] = np.nan
sepsis.loc[(
        (sepsis['septic_shock'] == 1)),
    'septic_shock_episode_number'] = septic_shock_episode_number

sepsis.loc[(
        (sepsis['had_septic_shock'] == 1)),
    'septic_shock_episode_number'] = (
        sepsis.loc[(
                (sepsis['had_septic_shock'] == 1)),
            'septic_shock_episode_number'].ffill())

sepsis['possible_septic_shock'] = (
    (sepsis['had_septic_shock'] == 1) &
    (sepsis['continued_sepsis_episode'] == 1) &
    ((sepsis['max_lactate'] > 2) |
        (sepsis['sofa_cardiovascular_score'] >= 3)))

end_septic_shock = sepsis[(
        (sepsis['possible_septic_shock'] == 1))].groupby(
            ['admissionid', 'septic_shock_episode_number']).agg(
    time=pd.NamedAgg(column='time', aggfunc='max')).reset_index()
end_septic_shock.drop(columns=['septic_shock_episode_number'], inplace=True)
end_septic_shock['continued_septic_shock'] = True
sepsis = pd.merge(
    sepsis, end_septic_shock, on=['admissionid', 'time'], how='left')

# Turn the entry immediately following a True value in 'end_sepsis_episode'
# from NaN to False.
sepsis.loc[(
        (sepsis['continued_septic_shock'].isna() == 1) &
        (sepsis['continued_septic_shock'] == 1).shift(1).fillna(False) == 1),
    'continued_septic_shock'] = False

sepsis.loc[(
        (sepsis['time'] == 0) & (sepsis['septic_shock'] == 0)),
    'continued_septic_shock'] = False

sepsis.loc[(
    (sepsis['septic_shock'] == 1)), 'continued_septic_shock'] = True
sepsis['continued_septic_shock'] = sepsis['continued_septic_shock'].ffill()
sepsis['continued_septic_shock'] = (
    sepsis['continued_septic_shock'].fillna(False))

sepsis.drop(columns=['possible_septic_shock'], inplace=True)

###############################################################################
# Cultures taken (for suspected infection)

if inputs.suspected_infection_by_cultures:
    # The following are used in the Amsterdam script for sepsis at admission
    cultures_itemid = [9189]  # Bloedkweken afnemen
    cultures_itemid += [9190]  # Cathetertipkweek afnemen
    cultures_itemid += [9200]  # Wondkweek afnemen
    cultures_itemid += [9202]  # Ascitesvochtkweek afnemen
    cultures_itemid += [9205]  # Legionella sneltest (urine)
    # The following are not used in the Amsterdam script
    cultures_itemid += [8097]  # Sputumkweek afnemen
    # cultures_v2_itemid += [9204]  # SDD-kweken afnemen
    # cultures_v2_itemid += [13024]  # SDD Inventarisatiekweken afnemen
    cultures_itemid += [9194]  # Liquorkweek afnemen
    cultures_itemid += [9192]  # Faeceskweek afnemen
    cultures_itemid += [8418]  # Urinekweek afnemen
    cultures_itemid += [9191]  # Drainvochtkweek afnemen
    cultures_itemid += [9203]  # Keelkweek afnemen
    cultures_itemid += [8588]  # MRSA kweken afnemen
    cultures_itemid += [9195]  # Neuskweek afnemen
    cultures_itemid += [9198]  # Rectumkweek afnemen
    cultures_itemid += [9197]  # Perineumkweek afnemen
    cultures_itemid += [9193]  # X-Kweek nader te bepalen

    cultures = procedureorderitems.loc[
        procedureorderitems['itemid'].isin(cultures_itemid)]
    cultures = cultures[['admissionid', 'time']]
    cultures = cultures.loc[(cultures.duplicated() == 0)]

    drugitems_cultures = pd.merge(
        drugitems_abx, cultures,
        on='admissionid', how='outer', suffixes=['', '_c'])

    # Keep only rows in which the culture was ordered within a day either side
    # of the antibiotics
    drugitems_cultures = drugitems_cultures.loc[(
        (drugitems_cultures['time'] -
            drugitems_cultures['time_c']).abs() <= 1)]

    drugitems_cultures['suspected'] = True

    drugitems_cultures_add = drugitems_cultures[
        ['admissionid', 'time', 'suspected']]
    drugitems_cultures_add = drugitems_cultures_add.loc[(
        (drugitems_cultures_add.duplicated() == 0))]

    # Merge on 'left' here, instead of 'outer'
    # This will ignore rows from drugitems_cultures_add that do not match those
    # in sepsis. As there is no corresponding row in sepsis, these would not be
    # marked as having sepsis under the culture and antibiotic definition,
    # since they automatically have no SOFA
    sepsis = pd.merge(
        sepsis, drugitems_cultures_add, on=['admissionid', 'time'], how='left')
    sepsis['suspected'] = sepsis['suspected'].fillna(False)

    sepsis['sepsis_episode_cultures'] = (
        (sepsis['suspected']) &
        ((sepsis['sofa_diff0'] >= 2) |
            (sepsis['sofa_diff1'] >= 2) |
            (sepsis['sofa_diff2'] >= 2)))

###############################################################################
# Save outputs

# This gives a binary variable for each patient if they had at least one
# sepsis episode.
sepsis_patients = sepsis.groupby(['admissionid']).agg(
        sepsis=pd.NamedAgg(column='sepsis_episode', aggfunc='any'),
        septic_shock=pd.NamedAgg(column='septic_shock', aggfunc='any')
    ).reset_index()

sepsis_columns = ['admissionid', 'time', 'discard', 'sofa_total_score']
sepsis_columns += ['antibiotic_escalation', 'sepsis_episode', 'septic_shock']
sepsis_output = sepsis[sepsis_columns]

# Save to csv
sofa.to_csv(inputs.output_file_path + 'sofa.csv', index=False)
sepsis_output.to_csv(inputs.output_file_path + 'sepsis3.csv', index=False)

# Including subscores
if inputs.save_intermediate_files:
    sepsis.to_csv(inputs.output_file_path + 'sepsis_all.csv', index=False)
    sofa_platelets.to_csv(
        inputs.output_file_path + 'sofa_platelets.csv', index=False)
    creatinine.to_csv(
        inputs.output_file_path + 'creatinine.csv', index=False)
    sofa_cardiovascular.to_csv(
        inputs.output_file_path + 'sofa_cardiovascular.csv', index=False)
    sofa_cns.to_csv(
        inputs.output_file_path + 'sofa_cns.csv', index=False)
    sofa_bilirubin.to_csv(
        inputs.output_file_path + 'sofa_bilirubin.csv', index=False)
    sofa_cardiovascular_meds.to_csv(
        inputs.output_file_path + 'sofa_cardiovascular_meds.csv', index=False)
    sofa_renal.to_csv(
        inputs.output_file_path + 'sofa_renal.csv', index=False)
    sofa_respiration.to_csv(
        inputs.output_file_path + 'sofa_respiration.csv', index=False)
    drugitems_abx.to_csv(
        inputs.output_file_path + 'drugitems_abx.csv', index=False)
    drugitems_abx_all.to_csv(
        inputs.output_file_path + 'drugitems_abx_all.csv', index=False)

###############################################################################
# Here are the SQL queries from AmsterdamUMCdb SOFA script, from which this
# code is developed.

# FiO2/PaO2 SQL
sql_sofa_respiration = """
WITH fio2_table AS (
    SELECT n.admissionid,
        n.measuredat,
        l.valueid,
        l.value AS O2_device,
        CASE
            WHEN n.itemid IN (
                --FiO2 settings on respiratory support
                6699, --FiO2 %: setting on Evita ventilator
                12279, --O2 concentratie --measurement by Servo-i/
                                                            Servo-U ventilator
                12369, --SET %O2: used with BiPap Vision ventilator
                16246 --Zephyros FiO2: Non-invasive ventilation
            ) THEN TRUE
            ELSE FALSE
        END AS ventilatory_support,
        CASE
            WHEN n.itemid IN (
                --FiO2 settings on respiratory support
                6699, --FiO2 %: setting on Evita ventilator
                12279, --O2 concentratie --measurement by Servo-i/
                                                            Servo-U ventilator
                12369, --SET %O2: used with BiPap Vision ventilator
                16246 --Zephyros FiO2: Non-invasive ventilation
            ) THEN
                CASE
                    WHEN NOT n.value IS NULL THEN n.value --use the settings
                    ELSE 0.21
                END
            ELSE -- estimate the FiO2
                CASE
                    WHEN l.valueid IN (
                        2, -- Nasaal
                        7 --O2-bril
                    ) THEN
                        CASE
                            WHEN n.value >= 1 AND n.value < 2 THEN 0.22
                            WHEN n.value >= 2 AND n.value < 3 THEN 0.25
                            WHEN n.value >= 3 AND n.value < 4 THEN 0.27
                            WHEN n.value >= 4 AND n.value < 5 THEN 0.30
                            WHEN n.value >= 5 THEN 0.35
                            ELSE 0.21
                        END
                    WHEN l.valueid IN (
                        1, --Diep Nasaal
                        3, --Kapje
                        8, --Kinnebak
                        9, --Nebulizer
                        4, --Kunstneus
                        18, --Spreekcanule
                        19 --Spreekklepje
                    ) THEN
                        CASE
                            WHEN n.value >= 1 AND n.value < 2 THEN 0.22
                                                        -- not defined by NICE
                            WHEN n.value >= 2 AND n.value < 3 THEN 0.25
                            WHEN n.value >= 3 AND n.value < 4 THEN 0.27
                            WHEN n.value >= 4 AND n.value < 5 THEN 0.30
                            WHEN n.value >= 5 AND n.value < 6 THEN 0.35
                            WHEN n.value >= 6 AND n.value < 7 THEN 0.40
                            WHEN n.value >= 7 AND n.value < 8 THEN 0.45
                            WHEN n.value >= 8 THEN 0.50
                            ELSE 0.21
                        END
                    WHEN l.valueid IN (
                        10, --Waterset
                        11, --Trach.stoma
                        13, --Ambu
                        14, --Guedel
                        15, --DL-tube
                        16, --CPAP
                        17 --Non-Rebreathing masker
                    ) THEN
                        CASE
                            WHEN n.value >= 6 AND n.value < 7 THEN 0.60
                            WHEN n.value >= 7 AND n.value < 8 THEN 0.70
                            WHEN n.value >= 8 AND n.value < 9 THEN 0.80
                            WHEN n.value >= 9 AND n.value < 10 THEN 0.85
                            WHEN n.value >= 10 THEN 0.90
                            ELSE 0.21
                        END
                    WHEN l.valueid IN (
                        12 --B.Lucht
                    ) THEN 0.21
                ELSE 0.21
            END
        END AS fio2
    FROM numericitems n
    LEFT JOIN admissions a ON
        n.admissionid = a.admissionid
    LEFT JOIN listitems l ON
        n.admissionid = l.admissionid AND
        n.measuredat = l.measuredat AND
        l.itemid = 8189 -- Toedieningsweg (Oxygen device)
    WHERE
        n.itemid IN (
            --Oxygen Flow settings without respiratory support
            8845, -- O2 l/min
            10387, --Zuurstof toediening (bloed)
            18587, --Zuurstof toediening

            --FiO2 settings on respiratory support
            6699, --FiO2 %: setting on Evita ventilator
            12279, --O2 concentratie --measurement by Servo-i/
                                                            Servo-U ventilator
            12369, --SET %O2: used with BiPap Vision ventilator
            16246 --Zephyros FiO2: Non-invasive ventilation
        )
    --measurements within 24 hours of ICU stay:
    AND (n.measuredat - a.admittedat) <= 1000*60*60*24 AND (n.measuredat -
                                                            a.admittedat) >= 0
    AND n.value > 0 --ignore stand by values from Evita ventilator
),
oxygenation AS (
    SELECT
        pao2.admissionid,
        CASE pao2.unitid
            WHEN 152 THEN pao2.value * 7.50061683 -- Conversion: kPa to mmHg
            ELSE pao2.value
        END AS pao2,
        f.value AS specimen_source,
        CASE
            WHEN NOT pao2.registeredby LIKE '%Systeem%' THEN TRUE
            ELSE FALSE
        END AS manual_entry,
        (pao2.measuredat - a.admittedat)/(1000*60) AS time,
        fio2_table.fio2,
        fio2_table.ventilatory_support,
        (fio2_table.measuredat - pao2.measuredat)/(60*1000) AS
                                                        FiO2_time_difference,
        ROW_NUMBER() OVER(
            PARTITION BY pao2.admissionid, pao2.measuredat
            ORDER BY ABS(fio2_table.measuredat - pao2.measuredat)
        ) AS priority --give priority to nearest FiO2 measurement
    FROM numericitems pao2
    LEFT JOIN admissions a ON
        pao2.admissionid = a.admissionid
    LEFT JOIN freetextitems f ON
        pao2.admissionid = f.admissionid AND
        pao2.measuredat = f.measuredat AND
        f.itemid = 11646 --Afname (bloed): source of specimen
    LEFT JOIN numericitems paco2 ON
        pao2.admissionid = paco2.admissionid AND
        pao2.measuredat = paco2.measuredat AND
        paco2.itemid IN (
            6846, --PCO2
            9990, --pCO2 (bloed)
            21213 --PCO2 (bloed) - kPa
        )
    LEFT JOIN fio2_table ON
        pao2.admissionid = fio2_table.admissionid AND
        fio2_table.measuredat > pao2.measuredat - 60*60*1000 AND
                        --no earlier than 60 minutes before pao2 measurement
        fio2_table.measuredat < pao2.measuredat + 15*60*1000
                        --no later than 15 minutes after pao2 measurement
    WHERE
        pao2.itemid IN (
            7433, --PO2
            9996, --PO2 (bloed)
            21214 --PO2 (bloed) - kPa
        )
    --measurements within 24 hours of ICU stay
            (use 30 minutes before admission to allow for time differences):
    AND (pao2.measuredat - a.admittedat) <= 1000*60*60*24 AND
        (pao2.measuredat - a.admittedat) >= -(1000*60*30) AND
    (f.value LIKE '%art.%' OR f.value IS NULL)
            -- source is arterial or undefined (assume arterial)
)
SELECT * FROM oxygenation
WHERE priority = 1
"""

# Platelets SQL
sql_sofa_platelets = """
SELECT
    n.admissionid,
    n.itemid,
    n.item,
    n.value,
    n.registeredby,
    CASE
        WHEN NOT n.registeredby LIKE '%Systeem%' THEN TRUE
        ELSE FALSE
    END AS manual_entry,
    (n.measuredat - a.admittedat)/(1000*60) AS time
FROM numericitems n
LEFT JOIN admissions a ON
    n.admissionid = a.admissionid
WHERE n.itemid IN (
    9964, --Thrombo's (bloed)
    6797, --Thrombocyten
    10409, --Thrombo's citr. bloed (bloed)
    14252 --Thrombo CD61 (bloed)
    )
--measurements within 24 hours of ICU stay
            (use 30 minutes before admission to allow for time differences):
AND (n.measuredat - a.admittedat) <= 1000*60*60*24 AND
    (n.measuredat - a.admittedat) >= -(1000*60*30)
"""

# Bilirubin SQL
sql_sofa_bilirubin = """
SELECT
    n.admissionid,
    n.itemid,
    n.item,
    n.value,
    n.registeredby,
    CASE
        WHEN NOT n.registeredby LIKE '%Systeem%' THEN TRUE
        ELSE FALSE
    END AS manual_entry,
    (n.measuredat - a.admittedat)/(1000*60) AS time
FROM numericitems n
LEFT JOIN admissions a ON
    n.admissionid = a.admissionid
    WHERE n.itemid IN (
        6813, --Bili Totaal
        9945 --Bilirubine (bloed)
    )
--measurements within 24 hours of ICU stay
        (use 30 minutes before admission to allow for time differences):
AND (n.measuredat - a.admittedat) <= 1000*60*60*24 AND
        (n.measuredat - a.admittedat) >= -(1000*60*30)
"""

# Cardiovascular SQL
sql_sofa_cardiovascular = """
WITH dosing AS (
    SELECT
        drugitems.admissionid,
        itemid,
        item,
        (start - admissions.admittedat)/(1000*60) AS start_time,
        (stop - admissions.admittedat)/(1000*60) AS stop_time,
        duration,
        rate,
        rateunit,
        dose,
        doseunit,
        doseunitid,
        doserateperkg,
        doserateunitid,
        doserateunit,
        CASE
            WHEN weightgroup LIKE '59' THEN 55
            WHEN weightgroup LIKE '60' THEN 65
            WHEN weightgroup LIKE '70' THEN 75
            WHEN weightgroup LIKE '80' THEN 85
            WHEN weightgroup LIKE '90' THEN 95
            WHEN weightgroup LIKE '100' THEN 105
            WHEN weightgroup LIKE '110' THEN 115
            ELSE 80 --mean weight for all years
        END as patientweight
    FROM drugitems
    LEFT JOIN admissions
    ON drugitems.admissionid = admissions.admissionid
    WHERE ordercategoryid = 65 -- continuous i.v. perfusor
    AND itemid IN (
            7179, -- Dopamine (Inotropin)
            7178, -- Dobutamine (Dobutrex)
            6818, -- Adrenaline (Epinefrine)
            7229  -- Noradrenaline (Norepinefrine)
        )
    AND rate > 0.1
)
SELECT
    admissionid,
    itemid,
    item,
    duration,
    rate,
    rateunit,
    start_time,
    stop_time,
    patientweight,
    dose,
    doserateperkg,
    doseunitid,
    doserateunitid
    CASE
    --recalculate the dose to µg/kg/min ('gamma')
    WHEN doserateperkg == B'0' AND doseunitid = 11 AND doserateunitid = 4
                                                --unit: µg/min -> µg/kg/min
        THEN CASE
            WHEN patientweight > 0
            THEN dose/patientweight
            ELSE dose/80 --mean weight
        END
    WHEN doserateperkg = B'0' AND doseunitid = 10 AND
    doserateunitid = 4 --unit: mg/min  -> µg/kg/min
        THEN CASE
            WHEN patientweight > 0
            THEN dose*1000/patientweight
            ELSE dose*1000/80 --mean weight
        END
    WHEN doserateperkg = B'0' AND doseunitid = 10 AND doserateunitid = 5
                                                --unit: mg/uur  -> µg/kg/min
        THEN CASE
            WHEN patientweight > 0
            THEN dose*1000/patientweight/60
            ELSE dose*1000/80/60 --mean weight
        END
    WHEN doserateperkg = B'1' AND doseunitid = 11 AND doserateunitid = 4
                                    --unit: µg/kg/min (no conversion needed)
        THEN dose
    WHEN doserateperkg = B'1' AND doseunitid = 11 AND doserateunitid = 5
                                            --unit: µg/kg/uur -> µg/kg/min
        THEN dose/60
    END AS gamma
    FROM dosing
    WHERE
    -- medication given within 24 hours of ICU stay:
    start_time <= 24*60 AND stop_time >= 0
    ORDER BY admissionid, start_time
"""

# Mean ABP SQL
sql_mean_abp = """
SELECT
    n.admissionid,
    n.itemid,
    n.item,
    n.value,
    CASE
        WHEN NOT registeredby IS NULL THEN TRUE
        ELSE FALSE
    END as validated,
    (measuredat - a.admittedat)/(1000*60) AS time
FROM numericitems n
LEFT JOIN admissions a ON
n.admissionid = a.admissionid
WHERE itemid IN (
    6642, --ABP gemiddeld
    6679, --Niet invasieve bloeddruk gemiddeld
    8843 --ABP gemiddeld II
)
AND (measuredat - a.admittedat) <= 1000*60*60*24 --measurements within 24 hours
"""

# Glasgow Coma Scale SQL
sql_gcs = """
WITH gcs_components AS (
    SELECT
        eyes.admissionid,
        --eyes.itemid,
        --eyes.item,
        --eyes.value,
        --eyes.valueid,
        CASE eyes.itemid
            WHEN 6732 THEN 5 - eyes.valueid     --Actief openen van de ogen
            WHEN 13077 THEN eyes.valueid        --A_Eye
            WHEN 14470 THEN eyes.valueid - 4    --RA_Eye
            WHEN 16628 THEN eyes.valueid - 4    --MCA_Eye
            WHEN 19635 THEN eyes.valueid - 4    --E_EMV_NICE_24uur
            WHEN 19638 THEN eyes.valueid - 8    --E_EMV_NICE_Opname
        END AS eyes_score,
        --motor.value,
        --motor.valueid,
        CASE motor.itemid
            WHEN 6734 THEN 7 - motor.valueid
                                        --Beste motore reactie van de armen
            WHEN 13072 THEN motor.valueid       --A_Motoriek
            WHEN 14476 THEN motor.valueid - 6   --RA_Motoriek
            WHEN 16634 THEN motor.valueid - 6   --MCA_Motoriek
            WHEN 19636 THEN motor.valueid - 6   --M_EMV_NICE_24uur
            WHEN 19639 THEN motor.valueid - 12  --M_EMV_NICE_Opname
        END AS motor_score,
        --verbal.value,
        --verbal.valueid,
        CASE verbal.itemid
            WHEN 6735 THEN 6 - verbal.valueid   --Beste verbale reactie
            WHEN 13066 THEN verbal.valueid      --A_Verbal
            WHEN 14482 THEN verbal.valueid - 5  --RA_Verbal
            WHEN 16640 THEN verbal.valueid - 5  --MCA_Verbal
            WHEN 19637 THEN verbal.valueid - 9 --V_EMV_NICE_24uur
            WHEN 19640 THEN verbal.valueid - 15 --V_EMV_NICE_Opname
        END AS verbal_score,
        eyes.registeredby,
        (eyes.measuredat - a.admittedat)/(1000*60) AS time
    FROM listitems eyes
    LEFT JOIN admissions a ON
        eyes.admissionid = a.admissionid
    LEFT JOIN listitems motor ON
        eyes.admissionid = motor.admissionid AND
        eyes.measuredat = motor.measuredat AND
        motor.itemid IN (
            6734, --Beste motore reactie van de armen
            13072, --A_Motoriek
            14476, --RA_Motoriek
            16634, --MCA_Motoriek
            19636, --M_EMV_NICE_24uur
            19639 --M_EMV_NICE_Opname
        )
    LEFT JOIN listitems verbal ON
        eyes.admissionid = verbal.admissionid AND
        eyes.measuredat = verbal.measuredat AND
        verbal.itemid IN (
            6735, --Beste verbale reactie
            13066, --A_Verbal
            14482, --RA_Verbal
            16640, --MCA_Verbal
            19637, --V_EMV_NICE_24uur
            19640 --V_EMV_NICE_Opname
        )
    WHERE eyes.itemid IN (
        6732, --Actief openen van de ogen
        13077, --A_Eye
        14470, --RA_Eye
        16628, --MCA_Eye
        19635, --E_EMV_NICE_24uur
        19638 --E_EMV_NICE_Opname
        )
    -- measurements within 24 hours of ICU stay:
    AND (eyes.measuredat - a.admittedat) <= 1000*60*60*24 AND
        (eyes.measuredat - a.admittedat) >= 0
), gcs AS (
    SELECT *,
        eyes_score + motor_score + (
            CASE
                WHEN verbal_score < 1 THEN 1
                ELSE verbal_score
            END
        ) AS gcs_score
    FROM gcs_components
), gcs_prioritized AS (
    SELECT *,
        ROW_NUMBER() OVER(
            PARTITION BY admissionid
            ORDER BY --orders by discipline
                CASE registeredby
                    WHEN 'ICV_Medisch Staflid' THEN 1
                    WHEN 'ICV_Medisch' THEN 2
                    WHEN 'ANES_Anesthesiologie'THEN 3
                    WHEN 'ICV_Physician assistant' THEN 4
                    WHEN 'ICH_Neurochirurgie'THEN 5
                    WHEN 'ICV_IC-Verpleegkundig' THEN 6
                    WHEN 'ICV_MC-Verpleegkundig' THEN 7
                    ELSE 8
                END, gcs_score
        ) AS priority
    FROM gcs
    ORDER BY admissionid, gcs_score ASC
)
SELECT *
FROM gcs_prioritized
WHERE priority = 1
"""

# Urine output SQL
sql_sofa_renal_urine_output = """
SELECT
    n.admissionid,
    n.itemid,
    n.item,
    n.value,
    n.registeredby,
    (n.measuredat - a.admittedat)/(1000*60) AS time
FROM numericitems n
LEFT JOIN admissions a ON
    n.admissionid = a.admissionid
WHERE n.itemid IN (
    8794, --UrineCAD
    8796, --UrineSupraPubis
    8798, --UrineSpontaan
    8800, --UrineIncontinentie
    8803, --UrineUP
    10743, --Nefrodrain li Uit
    10745, --Nefrodrain re Uit
    19921, --UrineSplint Li
    19922 --UrineSplint Re
    )
-- measurements within 24 hours of ICU stay
            (use 30 minutes before admission to allow for time differences):
AND (n.measuredat - a.admittedat) <= 1000*60*60*24 AND
    (n.measuredat - a.admittedat) >= 0
"""

# Creatinine SQL
sql_creatinine = """
WITH baseline AS (
    SELECT n.admissionid,
    MIN(n.value) AS baseline_creatinine
    FROM numericitems n
    LEFT JOIN admissions a ON
        n.admissionid = a.admissionid
    WHERE itemid IN (
        6836, --Kreatinine µmol/l (erroneously documented as µmol)
        9941, --Kreatinine (bloed) µmol/l
        14216 --KREAT enzym. (bloed) µmol/l
    ) AND
    --search 1 year befor admission
    (n.measuredat - a.admittedat)/(60*60*1000) > -(365*24) AND
    (n.measuredat - a.admittedat) < (24*60*60*1000)
    GROUP BY n.admissionid
),
max_creat AS (
    SELECT n.admissionid,
    MAX(n.value) AS max_creatinine_7days
    FROM numericitems n
    LEFT JOIN admissions a ON
        n.admissionid = a.admissionid
    WHERE itemid IN (
        6836, --Kreatinine µmol/l (erroneously documented as µmol)
        9941, --Kreatinine (bloed) µmol/l
        14216 --KREAT enzym. (bloed) µmol/l
    ) AND
    --within 7 days of admission
    (n.measuredat - a.admittedat) > 0 AND
    (n.measuredat - a.admittedat) < (7*24*60*60*1000)
    GROUP BY n.admissionid
)
SELECT
    n.admissionid,
    n.itemid,
    n.item,
    n.value,
    n.registeredby,
    CASE
        WHEN NOT n.registeredby LIKE '%Systeem%' THEN TRUE
        ELSE FALSE
    END AS manual_entry,
    (n.measuredat - a.admittedat)/(1000*60) AS time,
    b.baseline_creatinine,
    m.max_creatinine_7days,
    CASE
        -- AKI definition: 3 fold increase:
        WHEN baseline_creatinine > 0 AND
        m.max_creatinine_7days/baseline_creatinine > 3 THEN TRUE
        -- AKI definition: increase to >= 354 umol/l AND
        at least 44 umol/l increase:
        WHEN max_creatinine_7days >= 354 AND
        max_creatinine_7days - baseline_creatinine >= 44 THEN TRUE
        ELSE FALSE
    END AS acute_renal_failure
FROM numericitems n
LEFT JOIN admissions a ON
    n.admissionid = a.admissionid
LEFT JOIN baseline b ON -- get the baseline kreat (before admission)
    n.admissionid = b.admissionid
LEFT JOIN max_creat m ON --get the highest within a week of admission
    n.admissionid = m.admissionid
WHERE n.itemid IN (
    6836, --Kreatinine µmol/l (erroneously documented as µmol)
    9941, --Kreatinine (bloed) µmol/l
    14216 --KREAT enzym. (bloed) µmol/l
    )
-- measurements within 24 hours of ICU stay
            (use 30 minutes before admission to allow for time differences):
AND (n.measuredat - a.admittedat) <= 1000*60*60*24 AND
    (n.measuredat - a.admittedat) >= -(1000*60*30)
"""
