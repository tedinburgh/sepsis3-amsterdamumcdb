'''
Author: Tom Edinburgh
v3: date 29/06/2022.

This script is for descriptors of sepsis, following Shah et al. (Critical
Care Medicine, 2021)

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
from scipy.stats import kruskal, mannwhitneyu
from lifelines import CoxPHFitter

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
            in the script sepsis3_amsterdamumcdb.py (default: %(default)s)''',
        type=str)
    parser.add_argument(
        '--output_file_path',
        default='../../data/additional_files/tables/',
        help='''File path to the directory that will contain the outputs of
            this script, including all tables (default: %(default)s)''',
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

end_time = 1000*60*60*24*MAX_TIME if MAX_TIME is not None else None
start_time = 1000*60*60*24*MIN_TIME if MIN_TIME is not None else None
# numerics_sofa = numerics_read(
#     numerics_sofa_itemid + numerics_abp_itemid, admissions_df=admissions_df,
#     end=end_time, start=start_time)
# We need a baseline creatinine, so look back further.
# numerics_creatinine = numerics_read(
#     numerics_creatinine_itemid, admissions_df=admissions_df,
#     end=end_time, start=-1000*60*60*24*365)
# Lactate (for septic shock: max lactate of 2mmol/L + cardiovascular SOFA
# score of >=3 which corresponds to use of vasopressors)
lactate_itemids = [10053, 6837, 9580]
# numerics_lactate = numerics_read(
#     lactate_itemids, admissions_df=admissions_df,
#     end=end_time, start=start_time)

additional_itemid = [6640]  # Hartfrequentie # Heart rate
additional_itemid += [12311]  # O2-Saturatie (bloed)
additional_itemid += [11425]  # O2sat (overig)

# The following can be commented out if this script has already been run!*
numerics_additional = numerics_read(
    additional_itemid, admissions_df=admissions_df,
    end=1000*60*60*24*1, start=1000*60*60*24*0)
numerics_additional = numerics_additional.loc[(
   (numerics_additional['time'] == 0))]
numerics_additional.to_csv(
    inputs.additional_file_path + 'numerics_additional.csv')
# *Up to here
numerics_additional = pd.read_csv(
    inputs.additional_file_path + 'numerics_additional.csv')

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
# Load intermediate files

sofa = pd.read_csv(inputs.additional_file_path + 'sofa.csv')
sepsis_output = pd.read_csv(inputs.additional_file_path + 'sepsis3.csv')
# Including subscores
sepsis = pd.read_csv(inputs.additional_file_path + 'sepsis_all.csv')
sofa_platelets = pd.read_csv(
    inputs.additional_file_path + 'sofa_platelets.csv')
creatinine = pd.read_csv(inputs.additional_file_path + 'creatinine.csv')
sofa_cardiovascular = pd.read_csv(
    inputs.additional_file_path + 'sofa_cardiovascular.csv')
sofa_cns = pd.read_csv(inputs.additional_file_path + 'sofa_cns.csv')
sofa_respiration = pd.read_csv(
    inputs.additional_file_path + 'sofa_respiration.csv')
sofa_bilirubin = pd.read_csv(
    inputs.additional_file_path + 'sofa_bilirubin.csv')
sofa_cardiovascular_meds = pd.read_csv(
    inputs.additional_file_path + 'sofa_cardiovascular_meds.csv')
drugitems_abx = pd.read_csv(inputs.additional_file_path + 'drugitems_abx.csv')
drugitems_abx_all = pd.read_csv(
    inputs.additional_file_path + 'drugitems_abx_all.csv')
sofa_renal = pd.read_csv(inputs.additional_file_path + 'sofa_renal.csv')

sepsis_noradrenaline = pd.read_csv(
    inputs.additional_file_path + 'exclude_noradrenaline_6hr/sepsis_all.csv')

###############################################################################
# Mechanical ventilation

mv_evita1_itemid = [9534]  # Type beademing Evita 1
mv_evita1_valueid = [1]  # IPPV
mv_evita1_valueid += [2]  # IPPV_Assist
mv_evita1_valueid += [3]  # CPPV
mv_evita1_valueid += [4]  # CPPV_Assist
mv_evita1_valueid += [5]  # SIMV
mv_evita1_valueid += [6]  # SIMV_ASB
mv_evita1_valueid += [7]  # ASB
mv_evita1_valueid += [8]  # CPAP
mv_evita1_valueid += [9]  # CPAP_ASB
mv_evita1_valueid += [10]  # MMV
mv_evita1_valueid += [11]  # MMV_ASB
mv_evita1_valueid += [12]  # BIPAP
mv_evita1_valueid += [13]  # Pressure Controled

mv_evita4_itemid = [6685]  # Type Beademing Evita 4
mv_evita4_valueid = [1]  # CPPV
mv_evita4_valueid += [3]  # ASB
mv_evita4_valueid += [5]  # CPPV/ASSIST
mv_evita4_valueid += [6]  # SIMV/ASB
mv_evita4_valueid += [8]  # IPPV
mv_evita4_valueid += [9]  # IPPV/ASSIST
mv_evita4_valueid += [10]  # CPAP
mv_evita4_valueid += [11]  # CPAP/ASB
mv_evita4_valueid += [12]  # MMV
mv_evita4_valueid += [13]  # MMV/ASB
mv_evita4_valueid += [14]  # BIPAP
mv_evita4_valueid += [20]  # BIPAP-SIMV/ASB
mv_evita4_valueid += [22]  # BIPAP/ASB

mv_o2_itemid = [8189]  # Toedieningsweg O2
mv_o2_valueid = [16]  # CPAP

# Ventilatie Mode (Set) - Servo-I and Servo-U ventilators
mv_servo_itemid = [12290]
# Ventilatie Mode (Set) (2) Servo-I and Servo-U ventilators
mv_servo_itemid += [1234]
mv_servo_valueid = [2]  # PC
mv_servo_valueid += [3]  # VC
mv_servo_valueid += [4]  # PRVC
mv_servo_valueid += [5]  # VS
mv_servo_valueid += [6]  # SIMV(VC)+PS
mv_servo_valueid += [7]  # SIMV(PC)+PS
mv_servo_valueid += [8]  # PS/CPAP
mv_servo_valueid += [9]  # Bi Vente
mv_servo_valueid += [10]  # PC (No trig)
mv_servo_valueid += [11]  # VC (No trig)
mv_servo_valueid += [12]  # PRVC (No trig)
mv_servo_valueid += [13]  # PS/CPAP (trig)
mv_servo_valueid += [14]  # VC (trig)
mv_servo_valueid += [15]  # PRVC (trig)
mv_servo_valueid += [16]  # PC in NIV
mv_servo_valueid += [17]  # PS/CPAP in NIV
mv_servo_valueid += [18]  # NAVA

mv_bipap_itemid = [12376]  # Mode (Bipap Vision)
mv_bipap_valueid = [1]  # CPAP
mv_bipap_valueid += [2]  # BIPAP

mv_bool = (
    listitems['itemid'].isin(mv_evita1_itemid) &
    listitems['valueid'].isin(mv_evita1_valueid))
mv_bool |= (
    listitems['itemid'].isin(mv_evita4_itemid) &
    listitems['valueid'].isin(mv_evita4_valueid))
mv_bool |= (
    listitems['itemid'].isin(mv_o2_itemid) &
    listitems['valueid'].isin(mv_o2_valueid))
mv_bool |= (
    listitems['itemid'].isin(mv_servo_itemid) &
    listitems['valueid'].isin(mv_servo_valueid))
mv_bool |= (
    listitems['itemid'].isin(mv_bipap_itemid) &
    listitems['valueid'].isin(mv_bipap_valueid))

# Only for day of admission
mv_bool &= (listitems['time'] == 0)

mechanical_ventilation = listitems.loc[mv_bool].groupby('admissionid').agg(
        mv_modes=pd.NamedAgg(
            column='value', aggfunc=lambda x: '; '.join(x.unique()))
    ).reset_index()

mechanical_ventilation['mv_bool'] = True

###############################################################################
# Start to compile descriptors in one dataframe

descriptors = combined_diagnoses.copy()

descriptors['location_IC'] = descriptors['location'].str.contains('IC')
agegroups = descriptors['agegroup'].sort_values().unique()
for agegroup in agegroups:
    descriptors['agegroup_' + agegroup[:2]] = (
        descriptors['agegroup'] == agegroup)

descriptors['female'] = (descriptors['gender'] == 'Vrouw')
descriptors.loc[descriptors['gender'].isna(), 'female'] = np.nan

# Define ICU death as within 24hr of discharge (or if explicitly recorded)
descriptors['overleden_IC'] = (
    (descriptors['destination'] == 'Overleden') |
    (descriptors['dateofdeath'] -
        descriptors['dischargedat'] < (1000*60*60*24)))

descriptors['repeat_admission'] = (descriptors['admissioncount'] > 1)

descriptors['elective_surgical'] = (
    (descriptors['urgency'] == 0) & (descriptors['surgical'] == 1))
descriptors['emergency_surgical'] = (
    (descriptors['urgency'] == 1) & (descriptors['surgical'] == 1))
descriptors['emergency_medical'] = (
    (descriptors['urgency'] == 1) & (descriptors['surgical'] == 0))
# Check that these are mutually exclusive
admission_type = ['elective_surgical', 'emergency_surgical']
admission_type += ['emergency_medical']
if not descriptors[admission_type].sum(axis=1).all():
    print('Check admission types!')

# Update length of stay, to greater precision
descriptors['lengthofstay'] = (
    descriptors['dischargedat'] - descriptors['admittedat']) / (1000*60*60)

keep_columns = ['patientid', 'admissionid', 'admissioncount', 'admittedat']
keep_columns += ['dischargedat', 'lengthofstay', 'female', 'specialty']
keep_columns += ['agegroup_' + agegroup[:2] for agegroup in agegroups]
keep_columns += ['location_IC', 'overleden_IC', 'elective_surgical']
keep_columns += ['emergency_surgical', 'emergency_medical', 'repeat_admission']
keep_columns += ['sepsis', 'cardiosurg', 'respfailure', 'neurosurg', 'gisurg']
keep_columns += ['cardiacarrest', 'vascsurg', 'trauma', 'neuro', 'cardio']

descriptors = descriptors[keep_columns]
descriptors.rename(columns={'sepsis': 'old_sepsis'}, inplace=True)

admission_cat = ['cardiosurg', 'respfailure', 'neurosurg', 'gisurg']
admission_cat += ['cardiacarrest', 'vascsurg', 'trauma', 'neuro', 'cardio']
descriptors['admission_other'] = ~descriptors[admission_cat].any(axis=1)

sepsis_add_columns = ['admissionid', 'sepsis_episode', 'septic_shock']
sepsis_add_columns += ['antibiotic_escalation', 'sepsis_episode_cultures']
sepsis_add_columns += ['suspected', 'discard', 'sofa_total_score']
sepsis_add_columns += ['n_sofa_scores']
sepsis_add = sepsis.loc[(sepsis['time'] == 0), sepsis_add_columns]
descriptors = pd.merge(descriptors, sepsis_add, on='admissionid', how='left')
descriptors[sepsis_add_columns[1:-2]] = (
    descriptors[sepsis_add_columns[1:-2]].fillna(False))
descriptors['n_sofa_scores'] = descriptors['n_sofa_scores'].fillna(0)

descriptors.rename(columns={
    'sepsis_episode': 'sepsis',
    'new_sepsis_episode': 'new_sepsis',
    'sepsis_episode_cultures': 'sepsis_cultures'}, inplace=True)

drugitems_abx_all_add = drugitems_abx_all.loc[drugitems_abx_all['time'] == 0]
drugitems_abx_all_add = drugitems_abx_all_add[['admissionid', 'prophylactic']]
drugitems_abx_all_add = drugitems_abx_all_add.groupby('admissionid').agg(
        prophylactic=pd.NamedAgg(column='prophylactic', aggfunc='any')
    ).reset_index()
drugitems_abx_all_add['antibiotics_admission'] = True
descriptors = pd.merge(
    descriptors, drugitems_abx_all_add, on='admissionid', how='left')
descriptors['antibiotics_admission'] = (
    descriptors['antibiotics_admission'].fillna(False))
descriptors['prophylactic'] = descriptors['prophylactic'].fillna(False)

descriptors['sepsis_wo_shock'] = (
    ~descriptors['septic_shock'] & descriptors['sepsis'])
descriptors['antibiotics_wo_sepsis'] = (
    ~descriptors['sepsis'] & descriptors['antibiotics_admission'])
descriptors['antibiotic_escalation_wo_sepsis'] = (
    ~descriptors['sepsis'] & descriptors['antibiotic_escalation'])
descriptors['no_antibiotics'] = (
    ~descriptors['sepsis'] & ~descriptors['antibiotics_admission'])

descriptors['septic_shock_cultures'] = (
    descriptors['septic_shock'] & descriptors['sepsis_cultures'])
descriptors['sepsis_wo_shock_cultures'] = (
    ~descriptors['septic_shock_cultures'] & descriptors['sepsis_cultures'])
descriptors['antibiotics_wo_sepsis_cultures'] = (
    ~descriptors['sepsis_cultures'] & descriptors['antibiotics_admission'])
descriptors['antibiotic_escalation_wo_sepsis_cultures'] = (
    ~descriptors['sepsis_cultures'] & descriptors['antibiotic_escalation'])
descriptors['no_antibiotics_cultures'] = (
    ~descriptors['sepsis_cultures'] & ~descriptors['antibiotics_admission'])

sepsis_norad_add_columns = ['admissionid', 'sepsis_episode', 'septic_shock']
sepsis_noradrenaline_add = sepsis_noradrenaline.loc[
    (sepsis_noradrenaline['time'] == 0), sepsis_norad_add_columns]
sepsis_noradrenaline_add.rename(
    columns={'sepsis_episode': 'sepsis'}, inplace=True)
descriptors = pd.merge(
    descriptors, sepsis_noradrenaline_add, on='admissionid', how='left',
    suffixes=['', '_noradrenaline'])

descriptors['sepsis_noradrenaline'] = (
    descriptors['sepsis_noradrenaline'].fillna(False))
descriptors['septic_shock_noradrenaline'] = (
    descriptors['septic_shock_noradrenaline'].fillna(False))
descriptors['sepsis_wo_shock_noradrenaline'] = (
    ~descriptors['septic_shock_noradrenaline'] &
    descriptors['sepsis_noradrenaline'])
descriptors['antibiotics_wo_sepsis_noradrenaline'] = (
    ~descriptors['sepsis_noradrenaline'] &
    descriptors['antibiotics_admission'])
descriptors['no_antibiotics_noradrenaline'] = (
    ~descriptors['sepsis_noradrenaline'] &
    ~descriptors['antibiotics_admission'])

# numerics_additional contains values only in first 24hr after admission
heart_rate_add = numerics_additional[numerics_additional['itemid'] == 6640]
heart_rate_add = heart_rate_add.groupby('admissionid').agg(
        max_heart_rate=pd.NamedAgg(column='value', aggfunc='max')
    ).reset_index()
descriptors = pd.merge(
    descriptors, heart_rate_add, on='admissionid', how='left')

mean_abp_add = sofa_cardiovascular[(
    sofa_cardiovascular['itemid'].isin([6642, 6679, 8843]))]
mean_abp_add = mean_abp_add.groupby('admissionid').agg(
        min_mean_abp=pd.NamedAgg(column='lowest_mean_abp', aggfunc='min')
    ).reset_index()
descriptors = pd.merge(descriptors, mean_abp_add, on='admissionid', how='left')

# Again, we only want values from first 24hr after admission
respiration_add = sofa_respiration[sofa_respiration['time'] == 0]
respiration_add = respiration_add.groupby('admissionid').agg(
        max_fio2=pd.NamedAgg(column='fio2', aggfunc='max'),
        min_pao2=pd.NamedAgg(column='pao2', aggfunc='min'),
        min_pf_ratio=pd.NamedAgg(column='pf_ratio', aggfunc='min')
    ).reset_index()
descriptors = pd.merge(
    descriptors, respiration_add,
    on='admissionid', how='left')

spo2_add = numerics_additional[
    (numerics_additional['time'] == 0) &
    numerics_additional['itemid'].isin([12311, 11425])]
spo2_add = spo2_add.groupby('admissionid').agg(
        min_spo2=pd.NamedAgg(column='value', aggfunc='min')
    ).reset_index()
descriptors = pd.merge(descriptors, spo2_add, on='admissionid', how='left')

creatinine_add = creatinine.loc[creatinine['time'] == 0]
creatinine_add = creatinine_add.groupby('admissionid').agg(
        max_creatinine=pd.NamedAgg(column='max_creatinine', aggfunc='max')
    ).reset_index()
descriptors = pd.merge(
    descriptors, creatinine_add,
    on='admissionid', how='left')

gcs_add = sofa_cns.loc[sofa_cns['time'] == 0]
gcs_add = gcs_add.groupby('admissionid').agg(
        min_gcs=pd.NamedAgg(column='min_gcs', aggfunc='min')
    ).reset_index()
descriptors = pd.merge(descriptors, gcs_add, on='admissionid', how='left')

platelets_add = sofa_platelets[sofa_platelets['time'] == 0]
platelets_add = platelets_add.groupby('admissionid').agg(
        min_platelets=pd.NamedAgg(column='value', aggfunc='min')
    ).reset_index()
descriptors = pd.merge(
    descriptors, platelets_add, on='admissionid', how='left')

bilirubin_add = sofa_bilirubin[sofa_bilirubin['time'] == 0]
bilirubin_add = bilirubin_add.groupby('admissionid').agg(
        max_bilirubin=pd.NamedAgg(column='value', aggfunc='max')
    ).reset_index()
descriptors = pd.merge(
    descriptors, bilirubin_add, on='admissionid', how='left')

vasopressors = sofa_cardiovascular_meds[sofa_cardiovascular_meds['time'] == 0]
descriptors['vasopressors'] = descriptors['admissionid'].isin(
    vasopressors['admissionid'])

# Note, this is for non-prophylactic antibiotics!
antibiotics_add = drugitems_abx[
    (drugitems_abx['time'] >= 0) & (drugitems_abx['intravenous'])]
antibiotics_add = antibiotics_add[['admissionid', 'time']]
antibiotics_add = antibiotics_add.sort_values(by=['admissionid', 'time'])
antibiotics_add = antibiotics_add[~antibiotics_add.duplicated()]
# The idea here is to find the first point where the difference in days between
# two entries for an admission is greater than 1. Doing this is a bit
# convoluted: find the time difference, remove all values not equal to 1,
# insert a value of 1 if there is an entry on day 0. Cumulative summing this
# will give the number of days from day 0 until the first day without
# antibiotics, which will be NaN (using skipna=False). All times after this
# for that admission will also be NaN, so take the maximum to find the max
# number of consecutive days of antibiotics, started from admission.
antibiotics_add['time_diff'] = antibiotics_add['time'].diff(1)
antibiotics_add.loc[antibiotics_add['time_diff'] != 1, 'time_diff'] = np.nan
antibiotics_add.loc[(
        (antibiotics_add['admissionid'].duplicated() == 0) &
        (antibiotics_add['time'] == 0)),
    'time_diff'] = 1
antibiotics_add = antibiotics_add.groupby('admissionid').agg(
        n_days_total_abx=pd.NamedAgg(column='time', aggfunc='nunique'),
        n_days_consecutive_start_abx=pd.NamedAgg(
            column='time_diff', aggfunc=lambda x: x.cumsum(skipna=False).max())
    ).reset_index()
antibiotics_add['n_days_consecutive_start_abx'] = (
    antibiotics_add['n_days_consecutive_start_abx'].fillna(0).astype(int))
descriptors = pd.merge(
    descriptors, antibiotics_add, on='admissionid', how='left')
descriptors['antibiotics_min_4d_start'] = (
    (descriptors['n_days_consecutive_start_abx'] >= 4)).fillna(False)
descriptors['antibiotics_min_4d_departure'] = (
    (descriptors['n_days_total_abx'] >= 4) |
    (descriptors['n_days_total_abx'] == (
            1 + (descriptors['lengthofstay'] // 24)))).fillna(False)
descriptors['antibiotics_min_4d_death'] = (
    (descriptors['n_days_total_abx'] >= 4) |
    ((descriptors['n_days_total_abx'] == (  # number of days in ICU
            1 + (descriptors['lengthofstay'] // 24))) &
        (descriptors['overleden_IC'] == 1))).fillna(False)
descriptors.loc[
    descriptors['lengthofstay'] < 4*24, 'antibiotics_min_4d_death'] = False
descriptors['discharge_less_4d'] = (
    (descriptors['lengthofstay'] < 4*24) &
    (descriptors['overleden_IC'] == 0)).fillna(False)
descriptors['antibiotics_less_4d'] = (
    (descriptors['antibiotics_min_4d_death'] == 0) &
    (descriptors['discharge_less_4d'] == 0))
descriptors.drop(
    columns=['n_days_consecutive_start_abx', 'n_days_total_abx'], inplace=True)

mv_bool = mechanical_ventilation['mv_bool']
descriptors['mechanical_ventilation'] = descriptors['admissionid'].isin(
    mechanical_ventilation.loc[mv_bool, 'admissionid'])

# Specialty
descriptors['specialty_cardiosurg'] = (
    descriptors['specialty'] == 'Cardiochirurgie')
descriptors['specialty_neurosurg'] = (
    descriptors['specialty'] == 'Neurochirurgie')
descriptors['specialty_vascsurg'] = (
    descriptors['specialty'] == 'Vaatchirurgie')
descriptors['specialty_cardio'] = (descriptors['specialty'] == 'Cardiologie')
descriptors['specialty_gisurg'] = (
    descriptors['specialty'] == 'Heelkunde Gastro-enterologie')
descriptors['specialty_trauma'] = (descriptors['specialty'] == 'Traumatologie')
descriptors['specialty_neuro'] = (descriptors['specialty'] == 'Neurologie')
descriptors['specialty_respiratory'] = (
    (descriptors['specialty'] == 'Longziekte') |
    (descriptors['specialty'] == 'Heelkunde Longen/Oncologie'))
descriptors['specialty_oncology'] = (
    (descriptors['specialty'] == 'Oncologie Inwendig') |
    (descriptors['specialty'] == 'Heelkunde Oncologie'))
descriptors['specialty_internal'] = (
    descriptors['specialty'] == 'Inwendig')
descriptors['specialty_hematology'] = (
    descriptors['specialty'] == 'Hematologie')
descriptors['specialty_nephrology'] = (
    descriptors['specialty'] == 'Nefrologie')
descriptors['specialty_gynecology'] = (
    descriptors['specialty'] == 'Gynaecologie')
descriptors['specialty_urology'] = (
    descriptors['specialty'] == 'Urologie')
descriptors['specialty_ICU'] = (
    descriptors['specialty'] == 'Intensive Care Volwassenen')
descriptors['othersurg'] = (
    (descriptors.loc[:, (
        descriptors.columns.str.startswith('specialty_'))].sum(axis=1) == 0) &
    (descriptors['emergency_surgical'] | descriptors['elective_surgical']))
descriptors['othermedical'] = (
    (descriptors.loc[:, (
        descriptors.columns.str.startswith('specialty_'))].sum(axis=1) == 0) &
    (descriptors['emergency_medical']))

# The df 'descriptors' has inherited 'discard' from the df 'sepsis', add to it
descriptors['discard'] |= (descriptors['lengthofstay'] < 1)
descriptors['discard'] |= (descriptors['n_sofa_scores'] < 3)
descriptors['discard'] |= (
    descriptors['admissionid'].isin(sepsis['admissionid']) == 0)

###############################################################################
# Save descriptors dataframe
descriptors.to_csv(
    inputs.additional_file_path + 'descriptors.csv', index=False)

# Exclusion
descriptors_all = descriptors.copy()
descriptors = descriptors.loc[(descriptors['discard'] == 0)]

###############################################################################
# Tables


def n_fun(x):
    n = x.sum()
    pct = x.mean() * 100
    if n == 0:
        pct = 0
    return '{:.0f}'.format(n) + ' ({:.1f})'.format(pct)


def median_fun(x):
    if len(x) == 0:
        output = 'nan'
    else:
        m, lq, uq = x.quantile([0.5, 0.25, 0.75]).values
        mmin = min((abs(m), abs(lq), abs(uq)))
        format_str = '{:.' + str(max((0, int(np.ceil(1 - np.log10(mmin))))))
        format_str += 'f}'
        m = format_str.format(m)
        lq = format_str.format(lq)
        uq = format_str.format(uq)
        output = m + ' (' + lq + '-' + uq + ')'
    return output


def mean_fun(x):
    m = x.mean()
    sd = x.std()
    mmin = min((abs(m), abs(sd)))
    format_str = '{:.' + str(max((0, int(np.ceil(1 - np.log10(mmin)))))) + 'f}'
    m = format_str.format(m)
    sd = format_str.format(sd)
    return m + ' (' + sd + ')'


table_setup = pd.DataFrame(columns=['column', 'function', 'category'])
n_index = ['Number of admissions', 'Number of patients']
table_setup.loc['Number of admissions', 'column'] = 'admissionid'
table_setup.loc['Number of patients', 'column'] = 'patientid'
table_setup.loc[n_index, 'function'] = lambda x: x.nunique()
table_setup.loc[n_index, 'category'] = 'n'

demographic_index = ['Women, n (%)', '18-39', '40-49']
demographic_index += ['50-59', '60-69', '70-79', '80+']
table_setup.loc['Women, n (%)', 'column'] = 'female'
table_setup.loc['18-39', 'column'] = 'agegroup_18'
table_setup.loc['40-49', 'column'] = 'agegroup_40'
table_setup.loc['50-59', 'column'] = 'agegroup_50'
table_setup.loc['60-69', 'column'] = 'agegroup_60'
table_setup.loc['70-79', 'column'] = 'agegroup_70'
table_setup.loc['80+', 'column'] = 'agegroup_80'
table_setup.loc[demographic_index, 'function'] = n_fun
table_setup.loc[demographic_index, 'category'] = 'demographic'
demographic_index.insert(1, 'Age group, n (%)')

admission_cat_index = ['Elective surgical', 'Emergency surgical']
admission_cat_index += ['Emergency medical']
# admission_cat_index += ['Unknown']
table_setup.loc['Elective surgical', 'column'] = 'elective_surgical'
table_setup.loc['Emergency surgical', 'column'] = 'emergency_surgical'
table_setup.loc['Emergency medical', 'column'] = 'emergency_medical'
# table_setup.loc['Repeat admission', 'column'] = 'repeat_admission'
# table_setup.loc['Unknown', 'column'] = 'unknown_surgical_medical'
table_setup.loc[admission_cat_index, 'function'] = n_fun
table_setup.loc[admission_cat_index, 'category'] = 'admission_category'

admission_reason_index = ['Cardiothoracic surgery', 'Lung disease/surgery']
admission_reason_index += ['Neurosurgery', 'Trauma', 'Cardiology']
admission_reason_index += ['Vascular surgery', 'Gastrointestinal surgery']
admission_reason_index += ['Neurology', 'Internal', 'Hematology']
admission_reason_index += ['Gynaecology', 'Urology', 'From other ICU']
admission_reason_index += ['Other surgery', 'Other medical']
table_setup.loc['Cardiothoracic surgery', 'column'] = 'specialty_cardiosurg'
table_setup.loc['Lung disease/surgery', 'column'] = 'specialty_respiratory'
table_setup.loc['Neurosurgery', 'column'] = 'specialty_neurosurg'
table_setup.loc['Trauma', 'column'] = 'specialty_trauma'
table_setup.loc['Gastrointestinal surgery', 'column'] = 'specialty_gisurg'
table_setup.loc['Vascular surgery', 'column'] = 'specialty_vascsurg'
table_setup.loc['Cardiology', 'column'] = 'specialty_cardio'
table_setup.loc['Neurology', 'column'] = 'specialty_neuro'
table_setup.loc['Internal', 'column'] = 'specialty_internal'
table_setup.loc['Hematology', 'column'] = 'specialty_hematology'
table_setup.loc['Gynaecology', 'column'] = 'specialty_gynecology'
table_setup.loc['Urology', 'column'] = 'specialty_urology'
table_setup.loc['From other ICU', 'column'] = 'specialty_ICU'
table_setup.loc['Other surgery', 'column'] = 'othersurg'
table_setup.loc['Other medical', 'column'] = 'othermedical'
table_setup.loc[admission_reason_index, 'function'] = n_fun
table_setup.loc[admission_reason_index, 'category'] = 'admission_reason'
# admission_reason_index = ['Cardiothoracic', 'Respiratory failure']
# admission_reason_index += ['Neurosurgery', 'Trauma', 'Gastrointestinal']
# admission_reason_index += ['Vascular surgery', 'Cardiac arrest']
# admission_reason_index += ['Neurological disorder', 'Cardiac disorder']
# admission_reason_index += ['Other']
# table_setup.loc['Cardiothoracic', 'column'] = 'cardiosurg'
# table_setup.loc['Respiratory failure', 'column'] = 'respfailure'
# table_setup.loc['Neurosurgery', 'column'] = 'neurosurg'
# table_setup.loc['Trauma', 'column'] = 'trauma'
# table_setup.loc['Gastrointestinal surgery', 'column'] = 'gisurg'
# table_setup.loc['Vascular surgery', 'column'] = 'vascsurg'
# table_setup.loc['Cardiac arrest', 'column'] = 'cardiacarrest'
# table_setup.loc['Neurological disorder', 'column'] = 'neuro'
# table_setup.loc['Cardiac disorder', 'column'] = 'cardio'
# table_setup.loc['Other', 'column'] = 'admission_other'
# table_setup.loc[admission_reason_index, 'function'] = n_fun
# table_setup.loc[admission_reason_index, 'category'] = 'admission_reason'

physiology_index = ['Maximum heart rate']
physiology_index += ['Minimum mean arterial pressure, mmHg', 'Maximum FiO2']
physiology_index += ['Minimum SpO2', 'Minimum PaO2, mmHg']
physiology_index += ['Minimum PaO2:FiO2 ratio', 'Minimum GCS']
physiology_index += ['Maximum creatinine, µmol/L', 'Minimum platelets']
physiology_index += ['Maximum bilirubin, µmol/L', 'Maximum SOFA score']
physiology_index += ['Use of vasopressors, n (%)']
physiology_index += ['Mechanical ventilation, n (%)']
table_setup.loc['Maximum heart rate', 'column'] = 'max_heart_rate'
table_setup.loc[
    'Minimum mean arterial pressure, mmHg', 'column'] = 'min_mean_abp'
table_setup.loc['Maximum FiO2', 'column'] = 'max_fio2'
table_setup.loc['Minimum SpO2', 'column'] = 'min_spo2'
table_setup.loc['Minimum PaO2, mmHg', 'column'] = 'min_pao2'
table_setup.loc['Minimum PaO2:FiO2 ratio', 'column'] = 'min_pf_ratio'
table_setup.loc['Minimum GCS', 'column'] = 'min_gcs'
table_setup.loc['Maximum creatinine, µmol/L', 'column'] = 'max_creatinine'
table_setup.loc['Minimum platelets', 'column'] = 'min_platelets'
table_setup.loc['Maximum bilirubin, µmol/L', 'column'] = 'max_bilirubin'
table_setup.loc['Maximum SOFA score', 'column'] = 'sofa_total_score'
table_setup.loc['Use of vasopressors, n (%)', 'column'] = 'vasopressors'
table_setup.loc[
    'Mechanical ventilation, n (%)', 'column'] = 'mechanical_ventilation'
table_setup.loc[physiology_index, 'function'] = median_fun
table_setup.loc['Use of vasopressors, n (%)', 'function'] = n_fun
table_setup.loc['Mechanical ventilation, n (%)', 'function'] = n_fun
table_setup.loc[physiology_index, 'category'] = 'physiology'

outcomes_index = ['Antibiotic escalation, first 24hr, n (%)']
outcomes_index += ['IV antibiotics for first 4d (*), n (%)']
outcomes_index += ['IV antibiotics for at least 4d (**), n (%)']
outcomes_index += ['ICU length of stay, h, median (IQR)']
outcomes_index += ['ICU mortality, n (%)']
table_setup.loc[
    'Antibiotic escalation, first 24hr, n (%)',
    'column'] = 'antibiotic_escalation'
table_setup.loc[
    'IV antibiotics for first 4d (*), n (%)',
    'column'] = 'antibiotics_min_4d_start'
table_setup.loc[
    'IV antibiotics for at least 4d (**), n (%)',
    'column'] = 'antibiotics_min_4d_departure'
table_setup.loc[
    'ICU length of stay, h, median (IQR)', 'column'] = 'lengthofstay'
table_setup.loc['ICU mortality, n (%)', 'column'] = 'overleden_IC'
table_setup.loc[outcomes_index, 'function'] = n_fun
table_setup.loc['ICU length of stay, h, median (IQR)', 'function'] = median_fun
table_setup.loc[outcomes_index, 'category'] = 'outcome'

table_index = n_index + [''] + demographic_index
table_index += ['', 'Admission category, n (%)'] + admission_cat_index
table_index += ['', 'Specialty, n (%)'] + admission_reason_index
table_index += ['', 'First 24hr physiology, median (IQR)'] + physiology_index
table_index += ['', 'Outcomes'] + outcomes_index

# ICU stays
tableI_columns = ['Overall', 'Septic shock', 'Sepsis without shock']
tableI_columns += ['Antibiotics without sepsis', 'Not on antibiotics']

tableI = pd.DataFrame('', columns=tableI_columns, index=table_index)
missingnessI = pd.DataFrame(
    columns=['Missing, n', 'Missing, %'], index=table_index)
for row in table_index:
    if row in list(table_setup.index):
        temp_fun = table_setup.loc[row, 'function']
        temp_col = table_setup.loc[row, 'column']
        tableI.loc[row, 'Overall'] = temp_fun(
            descriptors.loc[descriptors['location_IC'], temp_col])
        # tableI.loc[row, 'Sepsis'] = temp_fun(
        #     descriptors.loc[(
        #             descriptors['sepsis'] &
        #             descriptors['location_IC']),
        #         temp_col])
        tableI.loc[row, 'Septic shock'] = temp_fun(
            descriptors.loc[(
                    descriptors['septic_shock'] &
                    descriptors['location_IC']),
                temp_col])
        tableI.loc[row, 'Sepsis without shock'] = temp_fun(
            descriptors.loc[(
                    descriptors['sepsis_wo_shock'] &
                    descriptors['location_IC']),
                temp_col])
        tableI.loc[row, 'Antibiotics without sepsis'] = temp_fun(
            descriptors.loc[(
                    descriptors['antibiotics_wo_sepsis'] &
                    descriptors['location_IC']),
                temp_col])
        tableI.loc[row, 'Not on antibiotics'] = temp_fun(
            descriptors.loc[(
                    descriptors['no_antibiotics'] &
                    descriptors['location_IC']),
                temp_col])
        missingnessI.loc[row, 'Missing, n'] = descriptors.loc[
            (descriptors['location_IC']), temp_col].isna().sum()
        missingnessI.loc[row, 'Missing, %'] = descriptors.loc[
            (descriptors['location_IC']), temp_col].isna().mean() * 100

# All admissions
tableIa = pd.DataFrame('', columns=tableI_columns, index=table_index)
missingnessIa = pd.DataFrame(
    columns=['Missing, n', 'Missing, %'], index=table_index)
for row in table_index:
    if row in list(table_setup.index):
        temp_fun = table_setup.loc[row, 'function']
        temp_col = table_setup.loc[row, 'column']
        tableIa.loc[row, 'Overall'] = temp_fun(descriptors.loc[:, temp_col])
        # tableIa.loc[row, 'Sepsis'] = temp_fun(
        #     descriptors.loc[descriptors['sepsis'], temp_col])
        tableIa.loc[row, 'Septic shock'] = temp_fun(
            descriptors.loc[descriptors['septic_shock'], temp_col])
        tableIa.loc[row, 'Sepsis without shock'] = temp_fun(
            descriptors.loc[descriptors['sepsis_wo_shock'], temp_col])
        tableIa.loc[row, 'Antibiotics without sepsis'] = temp_fun(
            descriptors.loc[descriptors['antibiotics_wo_sepsis'], temp_col])
        tableIa.loc[row, 'Not on antibiotics'] = temp_fun(
            descriptors.loc[descriptors['no_antibiotics'], temp_col])
        missingnessIa.loc[row, 'Missing, n'] = descriptors[
            temp_col].isna().sum()
        missingnessIa.loc[row, 'Missing, %'] = descriptors[
            temp_col].isna().mean() * 100

# MCU stays
tableIb = pd.DataFrame('', columns=tableI_columns, index=table_index)
missingnessIb = pd.DataFrame(
    columns=['Missing, n', 'Missing, %'], index=table_index)
for row in table_index:
    if row in list(table_setup.index):
        temp_fun = table_setup.loc[row, 'function']
        temp_col = table_setup.loc[row, 'column']
        tableIb.loc[row, 'Overall'] = temp_fun(
            descriptors.loc[~descriptors['location_IC'], temp_col])
        # tableIb.loc[row, 'Sepsis'] = temp_fun(
        #     descriptors.loc[(
        #             descriptors['sepsis'] &
        #             ~descriptors['location_IC']),
        #         temp_col])
        tableIb.loc[row, 'Septic shock'] = temp_fun(
            descriptors.loc[(
                    descriptors['septic_shock'] &
                    ~descriptors['location_IC']),
                temp_col])
        tableIb.loc[row, 'Sepsis without shock'] = temp_fun(
            descriptors.loc[(
                    descriptors['sepsis_wo_shock'] &
                    ~descriptors['location_IC']),
                temp_col])
        tableIb.loc[row, 'Antibiotics without sepsis'] = temp_fun(
            descriptors.loc[(
                    descriptors['antibiotics_wo_sepsis'] &
                    ~descriptors['location_IC']),
                temp_col])
        tableIb.loc[row, 'Not on antibiotics'] = temp_fun(
            descriptors.loc[(
                    descriptors['no_antibiotics'] &
                    ~descriptors['location_IC']),
                temp_col])
        missingnessIb.loc[row, 'Missing, n'] = descriptors.loc[
            (~descriptors['location_IC']), temp_col].isna().sum()
        missingnessIb.loc[row, 'Missing, %'] = descriptors.loc[
            (~descriptors['location_IC']), temp_col].isna().mean() * 100

# Sensitivity analysis: excluding noradrenaline
tableIc = pd.DataFrame('', columns=tableI_columns, index=table_index)
missingnessIc = pd.DataFrame(
    columns=['Missing, n', 'Missing, %'], index=table_index)
for row in table_index:
    if row in list(table_setup.index):
        temp_fun = table_setup.loc[row, 'function']
        temp_col = table_setup.loc[row, 'column']
        tableIc.loc[row, 'Overall'] = temp_fun(
            descriptors.loc[descriptors['location_IC'], temp_col])
        # tableIc.loc[row, 'Sepsis'] = temp_fun(
        #     descriptors.loc[(
        #             descriptors['sepsis_noradrenaline'] &
        #             descriptors['location_IC']),
        #         temp_col])
        tableIc.loc[row, 'Septic shock'] = temp_fun(
            descriptors.loc[(
                    descriptors['septic_shock_noradrenaline'] &
                    descriptors['location_IC']),
                temp_col])
        tableIc.loc[row, 'Sepsis without shock'] = temp_fun(
            descriptors.loc[(
                    descriptors['sepsis_wo_shock_noradrenaline'] &
                    descriptors['location_IC']),
                temp_col])
        tableIc.loc[row, 'Antibiotics without sepsis'] = temp_fun(
            descriptors.loc[(
                    descriptors['antibiotics_wo_sepsis_noradrenaline'] &
                    descriptors['location_IC']),
                temp_col])
        tableIc.loc[row, 'Not on antibiotics'] = temp_fun(
            descriptors.loc[(
                    descriptors['no_antibiotics_noradrenaline'] &
                    descriptors['location_IC']),
                temp_col])
        missingnessIc.loc[row, 'Missing, n'] = descriptors.loc[
            (descriptors['location_IC']), temp_col].isna().sum()
        missingnessIc.loc[row, 'Missing, %'] = descriptors.loc[
            (descriptors['location_IC']), temp_col].isna().mean() * 100

# Sensitivity analysis: cultures and antibiotics for suspected infection
tableId = pd.DataFrame('', columns=tableI_columns, index=table_index)
missingnessId = pd.DataFrame(
    columns=['Missing, n', 'Missing, %'], index=table_index)
for row in table_index:
    if row in list(table_setup.index):
        temp_fun = table_setup.loc[row, 'function']
        temp_col = table_setup.loc[row, 'column']
        tableId.loc[row, 'Overall'] = temp_fun(
            descriptors.loc[descriptors['location_IC'], temp_col])
        # tableId.loc[row, 'Sepsis'] = temp_fun(
        #     descriptors.loc[(
        #             descriptors['sepsis_cultures'] &
        #             descriptors['location_IC']),
        #         temp_col])
        tableId.loc[row, 'Septic shock'] = temp_fun(
            descriptors.loc[(
                    descriptors['septic_shock_cultures'] &
                    descriptors['location_IC']),
                temp_col])
        tableId.loc[row, 'Sepsis without shock'] = temp_fun(
            descriptors.loc[(
                    descriptors['sepsis_wo_shock_cultures'] &
                    descriptors['location_IC']),
                temp_col])
        tableId.loc[row, 'Antibiotics without sepsis'] = temp_fun(
            descriptors.loc[(
                    descriptors['antibiotics_wo_sepsis_cultures'] &
                    descriptors['location_IC']),
                temp_col])
        tableId.loc[row, 'Not on antibiotics'] = temp_fun(
            descriptors.loc[(
                    descriptors['no_antibiotics_cultures'] &
                    descriptors['location_IC']),
                temp_col])
        missingnessId.loc[row, 'Missing, n'] = descriptors.loc[
            (descriptors['location_IC']), temp_col].isna().sum()
        missingnessId.loc[row, 'Missing, %'] = descriptors.loc[
            (descriptors['location_IC']), temp_col].isna().mean() * 100

# ICU admissions
tableSII_columns = ['Antibiotics for at least 4 days or until death']
tableSII_columns += ['Discharged from ICU after < 4 days']
tableSII_columns += ['Antibiotics for < 4 days']
tableSII_columns += ['Overall (all patients with sepsis on admission)']

tableSII = pd.DataFrame('', columns=tableSII_columns, index=table_index)
missingnessSII = pd.DataFrame(
    columns=['Missing, n', 'Missing, %'], index=table_index)
for row in table_index:
    if row in list(table_setup.index):
        temp_fun = table_setup.loc[row, 'function']
        temp_col = table_setup.loc[row, 'column']
        temp_table_col = 'Overall (all patients with sepsis on admission)'
        tableSII.loc[row, temp_table_col] = temp_fun(
            descriptors.loc[(
                    descriptors['location_IC'] &
                    descriptors['sepsis']),
                temp_col])
        temp_table_col = 'Antibiotics for at least 4 days or until death'
        tableSII.loc[row, temp_table_col] = temp_fun(
            descriptors.loc[(
                    descriptors['location_IC'] &
                    descriptors['sepsis'] &
                    descriptors['antibiotics_min_4d_death']),
                temp_col])
        temp_table_col = 'Discharged from ICU after < 4 days'
        tableSII.loc[row, temp_table_col] = temp_fun(
            descriptors.loc[(
                    descriptors['location_IC'] &
                    descriptors['sepsis'] &
                    descriptors['discharge_less_4d']),
                temp_col])
        tableSII.loc[row, 'Antibiotics for < 4 days'] = temp_fun(
            descriptors.loc[(
                    descriptors['location_IC'] &
                    descriptors['sepsis'] &
                    (descriptors['antibiotics_less_4d'] |
                        descriptors['discharge_less_4d'])),
                temp_col])
        missingnessSII.loc[row, 'Missing, n'] = (
            descriptors.loc[
                descriptors['location_IC'] & descriptors['sepsis'], temp_col]
            ).isna().sum()
        missingnessSII.loc[row, 'Missing, %'] = (
            descriptors.loc[
                descriptors['location_IC'] & descriptors['sepsis'], temp_col]
            ).isna().mean() * 100

# All admissions
tableSII_index = table_index.copy()
tableSII_index.remove('Antibiotic escalation, first 24hr, n (%)')
tableSIIa = pd.DataFrame('', columns=tableSII_columns, index=table_index)
missingnessSIIa = pd.DataFrame(
    columns=['Missing, n', 'Missing, %'], index=table_index)
for row in table_index:
    if row in list(table_setup.index):
        temp_fun = table_setup.loc[row, 'function']
        temp_col = table_setup.loc[row, 'column']
        temp_table_col = 'Overall (all patients with sepsis on admission)'
        tableSIIa.loc[row, temp_table_col] = temp_fun(
            descriptors.loc[descriptors['sepsis'], temp_col])
        temp_table_col = 'Antibiotics for at least 4 days or until death'
        tableSIIa.loc[row, temp_table_col] = temp_fun(
            descriptors.loc[(
                    descriptors['sepsis'] &
                    descriptors['antibiotics_min_4d_death']),
                temp_col])
        temp_table_col = 'Discharged from ICU after < 4 days'
        tableSIIa.loc[row, temp_table_col] = temp_fun(
            descriptors.loc[(
                    descriptors['sepsis'] &
                    descriptors['discharge_less_4d']),
                temp_col])
        tableSIIa.loc[row, 'Antibiotics for < 4 days'] = temp_fun(
            descriptors.loc[(
                    descriptors['sepsis'] &
                    descriptors['antibiotics_less_4d']),
                temp_col])
        missingnessSIIa.loc[row, 'Missing, n'] = descriptors.loc[
            descriptors['sepsis'], temp_col].isna().sum()
        missingnessSIIa.loc[row, 'Missing, %'] = descriptors.loc[
            descriptors['sepsis'], temp_col].isna().mean() * 100

# MCU admissions
tableSIIb = pd.DataFrame('', columns=tableSII_columns, index=table_index)
missingnessSIIb = pd.DataFrame(
    columns=['Missing, n', 'Missing, %'], index=table_index)
for row in table_index:
    if row in list(table_setup.index):
        temp_fun = table_setup.loc[row, 'function']
        temp_col = table_setup.loc[row, 'column']
        temp_table_col = 'Overall (all patients with sepsis on admission)'
        tableSIIb.loc[row, temp_table_col] = temp_fun(
            descriptors.loc[(
                    ~descriptors['location_IC'] &
                    descriptors['sepsis']),
                temp_col])
        temp_table_col = 'Antibiotics for at least 4 days or until death'
        tableSIIb.loc[row, temp_table_col] = temp_fun(
            descriptors.loc[(
                    ~descriptors['location_IC'] &
                    descriptors['sepsis'] &
                    descriptors['antibiotics_min_4d_death']),
                temp_col])
        temp_table_col = 'Discharged from ICU after < 4 days'
        tableSIIb.loc[row, temp_table_col] = temp_fun(
            descriptors.loc[(
                    ~descriptors['location_IC'] &
                    descriptors['sepsis'] &
                    descriptors['discharge_less_4d']),
                temp_col])
        tableSIIb.loc[row, 'Antibiotics for < 4 days'] = temp_fun(
            descriptors.loc[(
                    ~descriptors['location_IC'] &
                    descriptors['sepsis'] &
                    (descriptors['antibiotics_less_4d'] |
                        descriptors['discharge_less_4d'])),
                temp_col])
        missingnessSIIb.loc[row, 'Missing, n'] = (
            descriptors.loc[
                ~descriptors['location_IC'] & descriptors['sepsis'], temp_col]
            ).isna().sum()
        missingnessSIIb.loc[row, 'Missing, %'] = (
            descriptors.loc[
                ~descriptors['location_IC'] & descriptors['sepsis'], temp_col]
            ).isna().mean() * 100

# Save tables
tableI['Missing data, n (%)'] = (
    missingnessI['Missing, n'].astype(str) + ' (' +
    missingnessI['Missing, %'].astype(float).round(1).astype(str) + ')')
tableIb['Missing data, n (%)'] = (
    missingnessIb['Missing, n'].astype(str) + ' (' +
    missingnessIb['Missing, %'].astype(float).round(1).astype(str) + ')')
tableIc['Missing data, n (%)'] = (
    missingnessIc['Missing, n'].astype(str) + ' (' +
    missingnessIc['Missing, %'].astype(float).round(1).astype(str) + ')')
tableId['Missing data, n (%)'] = (
    missingnessId['Missing, n'].astype(str) + ' (' +
    missingnessId['Missing, %'].astype(float).round(1).astype(str) + ')')
tableSII['Missing data, n (%)'] = (
    missingnessSII['Missing, n'].astype(str) + ' (' +
    missingnessSII['Missing, %'].astype(float).round(1).astype(str) + ')')
tableSIIb['Missing data, n (%)'] = (
    missingnessSIIb['Missing, n'].astype(str) + ' (' +
    missingnessSIIb['Missing, %'].astype(float).round(1).astype(str) + ')')

tableI.to_csv(inputs.output_file_path + 'icu/table_characteristics_icu.csv')
tableIb.to_csv(inputs.output_file_path + 'mcu/table_characteristics_mcu.csv')
tableIc.to_csv(
    inputs.output_file_path +
    'icu/table_characteristics_icu_exclude_noradrenaline_6hr.csv')
tableId.to_csv(
    inputs.output_file_path + 'icu/table_characteristics_icu_cultures.csv')
tableSII.to_csv(
    inputs.output_file_path + 'icu/table_characteristics_abx_icu.csv')
tableSIIb.to_csv(
    inputs.output_file_path + 'mcu/table_characteristics_abx_mcu.csv')

###############################################################################
# Antibiotics usage (all)

drugitems_abx_all['item'] = drugitems_abx_all['item'].apply(
    lambda x: x.replace(')', ')~~').split('~~')[0])
drugitems_abx_all['item'] = drugitems_abx_all['item'].apply(
    lambda x: x.replace(' ( ', ' ('))
drugitems_abx_all = drugitems_abx_all.loc[
    (drugitems_abx_all[
        ['admissionid', 'start_time', 'itemid']].duplicated() == 0)]
drugitems_abx_all = drugitems_abx_all.loc[
    drugitems_abx_all['admissionid'].isin(
        descriptors.loc[(descriptors['discard'] == 0), 'admissionid'])]
drugitems_abx_all.reset_index(inplace=True)

drugitems_abx_all.sort_values(
    by=['admissionid', 'item', 'start_time'], inplace=True)
drugitems_abx_all['time_diff'] = drugitems_abx_all['start_time'].diff(1)
drugitems_abx_all.loc[(
    (drugitems_abx_all['time_diff'] != 1)), 'time_diff'] = np.nan
drugitems_abx_all.loc[(
        (drugitems_abx_all[['admissionid', 'item']].duplicated() == 0)),
    'time_diff'] = 0
drugitems_abx_all['new_course'] = (
    (drugitems_abx_all['time_diff'].isna() == 0) &
    (drugitems_abx_all['time_diff'] != 1))
drugitems_abx_all.sort_values(by='index', inplace=True)

abx_usage_columns = ['admissionid', 'start_time', 'item', 'rank', 'new_course']
abx_usage = drugitems_abx_all[abx_usage_columns]

dutch_labels = abx_usage['item'].apply(
    lambda x: x.split('(')[0].strip())
english_labels = dutch_labels.copy()
str_replace_dict = {
    'Belcomycine': 'Colistin', 'Colistine Inhalatie': 'Colistin',
    'Metronidazol-Flagyl': 'Metronidazol',
    'Co-trimoxazol forte': 'Co-Trimoxazol', 'oïne': 'oin',
    'xon': 'xone', 'ine': 'in', 'dazol': 'dazole', 'xim': 'xime',
    'dim': 'dime', 'sulfaat': 'sulphate'}
for key in str_replace_dict.keys():
    value = str_replace_dict[key]
    english_labels = english_labels.str.replace(key, value)
abx_usage['item_en'] = english_labels
abx_usage.loc[(
        (abx_usage['item_en'] == 'Amoxicillin/Clavulaanzuur')),
    'item_en'] = 'Co-Amoxiclav'

abx_usage = pd.merge(
    abx_usage, descriptors[['admissionid', 'location_IC']],
    on='admissionid', how='left')

# ICU
abx_usage_icu = abx_usage.loc[(abx_usage['location_IC'] == 1)]
abx_total_usage_icu = abx_usage_icu.groupby('item_en').agg(
        total_courses=pd.NamedAgg(column='new_course', aggfunc='sum'),
        total_days=pd.NamedAgg(column='start_time', aggfunc='size'),
        rank=pd.NamedAgg(column='rank', aggfunc=lambda x: x.iloc[0]),
    ).reset_index()

abx_total_usage_icu = abx_total_usage_icu.sort_values(
    ['rank', 'total_courses'], ascending=[False, False]).reset_index(drop=True)
abx_total_usage_icu['item_en'] = abx_total_usage_icu['item_en'].apply(
    lambda x: '  ' + x.split('(')[0].strip())

totals = abx_total_usage_icu['total_courses'].sum()
abx_total_usage_icu['percent_courses'] = (
    abx_total_usage_icu['total_courses'] / totals)
abx_total_usage_icu['percent_courses'] = (
    abx_total_usage_icu['percent_courses'].apply(
        lambda x: '{:.2f}'.format(x * 100)))

for col in abx_total_usage_icu.columns:
    if col.startswith('total'):
        abx_total_usage_icu[col] = abx_total_usage_icu[col].apply(
            lambda x: '{:.0f}'.format(x))

for rank in abx_total_usage_icu['rank'].sort_values().unique():
    idx = (abx_total_usage_icu['rank'] == rank).idxmax()
    n_columns = len(abx_total_usage_icu.columns)
    insert_row = pd.Series(
        ['Rank ' + str(rank) + ':'] + [''] * (n_columns - 1),
        index=abx_total_usage_icu.columns,
        name='')
    abx_total_usage_icu = abx_total_usage_icu.iloc[:idx, :].append(
        insert_row).append(abx_total_usage_icu.iloc[idx:, :])

abx_total_usage_icu.reset_index(drop=True, inplace=True)
abx_total_columns = {
    'item_en': 'Antibiotic',
    'total_courses': 'Total number of courses',
    'percent_courses': 'Percentage of courses (%)',
    'total_days': 'Total number of days'}
abx_total_usage_icu = abx_total_usage_icu[abx_total_columns.keys()]
abx_total_usage_icu = abx_total_usage_icu.rename(columns=abx_total_columns)
abx_total_usage_icu.set_index('Antibiotic', inplace=True)

# MCU
abx_usage_mcu = abx_usage.loc[(abx_usage['location_IC'] == 0)]
abx_total_usage_mcu = abx_usage_mcu.groupby('item_en').agg(
        total_courses=pd.NamedAgg(column='new_course', aggfunc='sum'),
        total_days=pd.NamedAgg(column='start_time', aggfunc='size'),
        rank=pd.NamedAgg(column='rank', aggfunc=lambda x: x.iloc[0]),
    ).reset_index()

abx_total_usage_mcu = abx_total_usage_mcu.sort_values(
    ['rank', 'total_courses'], ascending=[False, False]).reset_index(drop=True)
abx_total_usage_mcu['item_en'] = abx_total_usage_mcu['item_en'].apply(
    lambda x: '  ' + x.split('(')[0].strip())

totals = abx_total_usage_mcu['total_courses'].sum()
abx_total_usage_mcu['percent_courses'] = (
    abx_total_usage_mcu['total_courses'] / totals)
abx_total_usage_mcu['percent_courses'] = (
    abx_total_usage_mcu['percent_courses'].apply(
        lambda x: '{:.2f}'.format(x * 100)))

for col in abx_total_usage_mcu.columns:
    if col.startswith('total'):
        abx_total_usage_mcu[col] = abx_total_usage_mcu[col].apply(
            lambda x: '{:.0f}'.format(x))

for rank in abx_total_usage_mcu['rank'].sort_values().unique():
    idx = (abx_total_usage_mcu['rank'] == rank).idxmax()
    n_columns = len(abx_total_usage_mcu.columns)
    insert_row = pd.Series(
        ['Rank ' + str(rank) + ':'] + [''] * (n_columns - 1),
        index=abx_total_usage_mcu.columns,
        name='')
    abx_total_usage_mcu = abx_total_usage_mcu.iloc[:idx, :].append(
        insert_row).append(abx_total_usage_mcu.iloc[idx:, :])

abx_total_usage_mcu.reset_index(drop=True, inplace=True)
abx_total_usage_mcu = abx_total_usage_mcu[abx_total_columns.keys()]
abx_total_usage_mcu = abx_total_usage_mcu.rename(columns=abx_total_columns)
abx_total_usage_mcu.set_index('Antibiotic', inplace=True)

abx_total_usage_icu.to_csv(
    inputs.output_file_path + 'icu/table_abx_total_usage_icu.csv')
abx_total_usage_mcu.to_csv(
    inputs.output_file_path + 'mcu/table_abx_total_usage_mcu.csv')

abx_usage.to_csv(inputs.additional_file_path + 'abx_usage.csv', index=False)

###############################################################################
# Antibiotics usage (non-prophylactic only)

abx_usage_columns = ['admissionid', 'start_time', 'item', 'rank', 'new_course']
abx_usage = drugitems_abx_all.loc[(
    (drugitems_abx_all['prophylactic'] == 0)), abx_usage_columns]

abx_usage = pd.merge(
    abx_usage, descriptors[['admissionid', 'location_IC']],
    on='admissionid', how='left')

dutch_labels = abx_usage['item'].apply(
    lambda x: x.split('(')[0].strip())
english_labels = dutch_labels.copy()
str_replace_dict = {
    'Belcomycine': 'Colistin', 'Colistine Inhalatie': 'Colistin',
    'Metronidazol-Flagyl': 'Metronidazol',
    'Co-trimoxazol forte': 'Co-Trimoxazol', 'oïne': 'oin',
    'xon': 'xone', 'ine': 'in', 'dazol': 'dazole', 'xim': 'xime',
    'dim': 'dime', 'sulfaat': 'sulphate'}
for key in str_replace_dict.keys():
    value = str_replace_dict[key]
    english_labels = english_labels.str.replace(key, value)
abx_usage['item_en'] = english_labels
abx_usage.loc[(
        (abx_usage['item_en'] == 'Amoxicillin/Clavulaanzuur')),
    'item_en'] = 'Co-Amoxiclav'

# ICU
abx_usage_icu = abx_usage.loc[(abx_usage['location_IC'] == 1)]
abx_total_usage_icu = abx_usage_icu.groupby('item_en').agg(
        total_courses=pd.NamedAgg(column='new_course', aggfunc='sum'),
        total_days=pd.NamedAgg(column='start_time', aggfunc='size'),
        rank=pd.NamedAgg(column='rank', aggfunc=lambda x: x.iloc[0]),
    ).reset_index()

abx_total_usage_icu = abx_total_usage_icu.sort_values(
    ['rank', 'total_courses'], ascending=[False, False]).reset_index(drop=True)
abx_total_usage_icu['item_en'] = abx_total_usage_icu['item_en'].apply(
    lambda x: '  ' + x.split('(')[0].strip())

totals = abx_total_usage_icu['total_courses'].sum()
abx_total_usage_icu['percent_courses'] = (
    abx_total_usage_icu['total_courses'] / totals)
abx_total_usage_icu['percent_courses'] = (
    abx_total_usage_icu['percent_courses'].apply(
        lambda x: '{:.2f}'.format(x * 100)))

for col in abx_total_usage_icu.columns:
    if col.startswith('total'):
        abx_total_usage_icu[col] = abx_total_usage_icu[col].apply(
            lambda x: '{:.0f}'.format(x))

for rank in abx_total_usage_icu['rank'].sort_values().unique():
    idx = (abx_total_usage_icu['rank'] == rank).idxmax()
    n_columns = len(abx_total_usage_icu.columns)
    insert_row = pd.Series(
        ['Rank ' + str(rank) + ':'] + [''] * (n_columns - 1),
        index=abx_total_usage_icu.columns,
        name='')
    abx_total_usage_icu = abx_total_usage_icu.iloc[:idx, :].append(
        insert_row).append(abx_total_usage_icu.iloc[idx:, :])

abx_total_usage_icu.reset_index(drop=True, inplace=True)
abx_total_columns = {
    'item_en': 'Antibiotic',
    'total_courses': 'Total number of courses',
    'percent_courses': 'Percentage of courses (%)',
    'total_days': 'Total number of antibiotic-days'}
abx_total_usage_icu = abx_total_usage_icu[abx_total_columns.keys()]
abx_total_usage_icu = abx_total_usage_icu.rename(columns=abx_total_columns)
abx_total_usage_icu.set_index('Antibiotic', inplace=True)

# MCU
abx_usage_mcu = abx_usage.loc[(abx_usage['location_IC'] == 0)]
abx_total_usage_mcu = abx_usage_mcu.groupby('item_en').agg(
        total_courses=pd.NamedAgg(column='new_course', aggfunc='sum'),
        total_days=pd.NamedAgg(column='start_time', aggfunc='size'),
        rank=pd.NamedAgg(column='rank', aggfunc=lambda x: x.iloc[0]),
    ).reset_index()

abx_total_usage_mcu = abx_total_usage_mcu.sort_values(
    ['rank', 'total_courses'], ascending=[False, False]).reset_index(drop=True)
abx_total_usage_mcu['item_en'] = abx_total_usage_mcu['item_en'].apply(
    lambda x: '  ' + x.split('(')[0].strip())

totals = abx_total_usage_mcu['total_courses'].sum()
abx_total_usage_mcu['percent_courses'] = (
    abx_total_usage_mcu['total_courses'] / totals)
abx_total_usage_mcu['percent_courses'] = (
    abx_total_usage_mcu['percent_courses'].apply(
        lambda x: '{:.2f}'.format(x * 100)))

for col in abx_total_usage_mcu.columns:
    if col.startswith('total'):
        abx_total_usage_mcu[col] = abx_total_usage_mcu[col].apply(
            lambda x: '{:.0f}'.format(x))

for rank in abx_total_usage_mcu['rank'].sort_values().unique():
    idx = (abx_total_usage_mcu['rank'] == rank).idxmax()
    n_columns = len(abx_total_usage_mcu.columns)
    insert_row = pd.Series(
        ['Rank ' + str(rank) + ':'] + [''] * (n_columns - 1),
        index=abx_total_usage_mcu.columns,
        name='')
    abx_total_usage_mcu = abx_total_usage_mcu.iloc[:idx, :].append(
        insert_row).append(abx_total_usage_mcu.iloc[idx:, :])

abx_total_usage_mcu.reset_index(drop=True, inplace=True)
abx_total_usage_mcu = abx_total_usage_mcu[abx_total_columns.keys()]
abx_total_usage_mcu = abx_total_usage_mcu.rename(columns=abx_total_columns)
abx_total_usage_mcu.set_index('Antibiotic', inplace=True)

abx_total_usage_icu.to_csv(
    inputs.output_file_path + 'icu/table_abx_non_prophylactic_usage_icu.csv')
abx_total_usage_mcu.to_csv(
    inputs.output_file_path + 'mcu/table_abx_non_prophylactic_usage_mcu.csv')

abx_usage.to_csv(
    inputs.additional_file_path + 'abx_non_prophylactic_usage.csv',
    index=False)

###############################################################################
# Numbers for the text

descriptors = descriptors_all.copy()

text_n = pd.DataFrame(columns=['Description', 'n', 'denom', 'percent'])
text_n.loc[0] = [
    'Starting admissions', descriptors.shape[0], np.nan, np.nan]
text_n.loc[1] = [
    'LoS <1hr', (descriptors['lengthofstay'] < 1).sum(),
    descriptors.shape[0], (descriptors['lengthofstay'] < 1).mean()*100]
descriptors = descriptors.loc[(descriptors['lengthofstay'] >= 1)]
text_n.loc[2] = [
    'Fewer than 3 SOFA', (descriptors['n_sofa_scores'] < 3).sum(),
    descriptors.shape[0], (descriptors['n_sofa_scores'] < 3).mean()*100]
text_n.loc[3] = [
    'Remaining ICU/MCU', (descriptors['n_sofa_scores'] >= 3).sum(),
    descriptors.shape[0],
    (descriptors['n_sofa_scores'] >= 3).mean()*100]
descriptors = descriptors.loc[(descriptors['n_sofa_scores'] >= 3)]
text_n.loc[4] = [
    'MCU, excluded', (descriptors['location_IC'] == 0).sum(),
    descriptors.shape[0], (descriptors['location_IC'] == 0).mean()*100]
descriptors = descriptors.loc[(descriptors['location_IC'] == 1)]
text_n.loc[5] = [
    'Remaining ICU admissions', descriptors.shape[0], np.nan, np.nan]
text_n.loc[6] = [
    'Remaining ICU patients', descriptors['patientid'].nunique(),
    np.nan, np.nan]
sepsis = sepsis.loc[sepsis['admissionid'].isin(descriptors['admissionid'])]
sepsis['max_sofa_increase'] = sepsis[
    ['sofa_diff0', 'sofa_diff1', 'sofa_diff2']].max(axis=1)
text_n.loc[7] = [
    'SOFA increase >=2', (sepsis['max_sofa_increase'] >= 2).sum(),
    np.nan, np.nan]
text_n.loc[8] = [
    'Associated with antibiotic escalation',
    ((sepsis['max_sofa_increase'] >= 2) &
        (sepsis['antibiotic_escalation'] == 1)).sum(),
    (sepsis['max_sofa_increase'] >= 2).sum(),
    (sepsis.loc[(
            (sepsis['max_sofa_increase'] >= 2)),
        'antibiotic_escalation'] == 1).mean()*100]
text_n.loc[9] = [
    'Sepsis episodes (including during SDD when cefotaxime continued >4d)',
    (sepsis['sepsis_episode'] == 1).sum(), np.nan, np.nan]
text_n.loc[10] = [
    'New sepsis episode', (sepsis['new_sepsis_episode'] == 1).sum(),
    (sepsis['sepsis_episode'] == 1).sum(),
    (sepsis.loc[(
            (sepsis['sepsis_episode'] == 1)),
        'new_sepsis_episode'] == 1).mean()*100]
text_n.loc[11] = [
    'Septic shock', (sepsis['septic_shock'] == 1).sum(),
    (sepsis['new_sepsis_episode'] == 1).sum(),
    (sepsis.loc[(
            (sepsis['new_sepsis_episode'] == 1)),
        'septic_shock'] == 1).mean()*100]
text_n.loc[12] = [
    'Sepsis on admission',
    (sepsis.loc[(sepsis['new_sepsis_episode'] == 1), 'time'] == 0).sum(),
    (sepsis['new_sepsis_episode'] == 1).sum(),
    (sepsis.loc[(
            (sepsis['new_sepsis_episode'] == 1)),
        'time'] == 0).mean()*100]
text_n.loc[13] = [
    'Septic shock on admission',
    (sepsis.loc[(sepsis['time'] == 0), 'septic_shock'] == 1).sum(),
    (sepsis.loc[(sepsis['time'] == 0), 'new_sepsis_episode'] == 1).sum(),
    (sepsis.loc[(
            (sepsis['new_sepsis_episode'] == 1) & (sepsis['time'] == 0)),
        'septic_shock'] == 1).mean()*100]
text_n.loc[14] = [
    'Medical x sepsis',
    ((descriptors['emergency_medical'] == 1) &
        (descriptors['sepsis'] == 1)).sum(),
    (descriptors['emergency_medical'] == 1).sum(),
    (descriptors.loc[(
            (descriptors['emergency_medical'] == 1)),
        'sepsis'] == 1).mean()*100]
text_n.loc[15] = [
    'Medical x septic shock',
    ((descriptors['emergency_medical'] == 1) &
        (descriptors['septic_shock'] == 1)).sum(),
    ((descriptors['emergency_medical'] == 1) &
        (descriptors['sepsis'] == 1)).sum(),
    (descriptors.loc[(
            (descriptors['emergency_medical'] == 1) &
            (descriptors['sepsis'] == 1)),
        'septic_shock'] == 1).mean()*100]
text_n.loc[16] = [
    'Emergency surgery x sepsis',
    ((descriptors['emergency_surgical'] == 1) &
        (descriptors['sepsis'] == 1)).sum(),
    (descriptors['emergency_surgical'] == 1).sum(),
    (descriptors.loc[(
            (descriptors['emergency_surgical'] == 1)),
        'sepsis'] == 1).mean()*100]
text_n.loc[17] = [
    'Emergency surgery x septic shock',
    ((descriptors['emergency_surgical'] == 1) &
        (descriptors['septic_shock'] == 1)).sum(),
    ((descriptors['emergency_surgical'] == 1) &
        (descriptors['sepsis'] == 1)).sum(),
    (descriptors.loc[(
            (descriptors['emergency_surgical'] == 1) &
            (descriptors['sepsis'] == 1)),
        'septic_shock'] == 1).mean()*100]
text_n.loc[18] = [
    'Sepsis without shock x IV antibiotics >4d',
    ((descriptors['antibiotics_min_4d_departure'] == 1) &
        (descriptors['sepsis_wo_shock'] == 1)).sum(),
    (descriptors['sepsis_wo_shock'] == 1).sum(),
    (descriptors.loc[(
            (descriptors['sepsis_wo_shock'] == 1)),
        'antibiotics_min_4d_departure'] == 1).mean()*100]
text_n.loc[19] = [
    'Septic shock x IV antibiotics >4d',
    ((descriptors['antibiotics_min_4d_departure'] == 1) &
        (descriptors['septic_shock'] == 1)).sum(),
    (descriptors['septic_shock'] == 1).sum(),
    (descriptors.loc[(
            (descriptors['septic_shock'] == 1)),
        'antibiotics_min_4d_departure'] == 1).mean()*100]
n_sepsis = sepsis[
    ['admissionid', 'new_sepsis_episode']].groupby('admissionid').sum()
text_n.loc[20] = [
    '2 separate sepsis episodes',
    (n_sepsis.values > 1).sum(),
    n_sepsis.size, (n_sepsis.values > 1).mean()*100]
text_n.loc[21] = [
    '>2 separate sepsis episodes',
    (n_sepsis.values > 2).sum(),
    n_sepsis.size, (n_sepsis.values > 1).mean()*100]
n_septic_shock = sepsis[
    ['admissionid', 'new_septic_shock']].groupby('admissionid').sum()
text_n.loc[22] = [
    '>1 separate septic shock episodes',
    (n_septic_shock.values > 1).sum(),
    n_septic_shock.size, (n_septic_shock.values > 1).mean()*100]
abx_usage_icu = abx_usage_icu.loc[
    abx_usage_icu['admissionid'].isin(
        descriptors.loc[(descriptors['discard'] == 0), 'admissionid'])]
n_regimens = abx_usage_icu.groupby('admissionid').agg(
    n_items=pd.NamedAgg(column='item', aggfunc=lambda x: x.nunique()))
text_n.loc[23] = [
    '0 abx regimen (during all ICU stay)',
    descriptors.shape[0] - n_regimens.shape[0], np.nan, np.nan]
text_n.loc[24] = ['1 abx regimen', (n_regimens == 1).sum()[0], np.nan, np.nan]
text_n.loc[25] = ['>1 abx regimen', (n_regimens > 1).sum()[0], np.nan, np.nan]
text_n.loc[26] = [
    'Abx courses', abx_usage_icu['new_course'].sum(), np.nan, np.nan]
text_n.loc[27] = [
    'Abx days', abx_usage_icu.shape[0], np.nan, np.nan]
text_n.loc[28] = [
    'Rank 1 courses',
    abx_usage_icu.loc[(abx_usage_icu['rank'] == 1), 'new_course'].sum(),
    abx_usage_icu['new_course'].sum(),
    (abx_usage_icu.loc[(
            (abx_usage_icu['new_course'] == 1)),
        'rank'] == 1).mean()*100]
text_n.loc[29] = [
    'Rank 2 courses',
    abx_usage_icu.loc[(abx_usage_icu['rank'] == 2), 'new_course'].sum(),
    abx_usage_icu['new_course'].sum(),
    (abx_usage_icu.loc[(
            (abx_usage_icu['new_course'] == 1)),
        'rank'] == 2).mean()*100]
text_n.loc[30] = [
    'Rank 3 courses',
    abx_usage_icu.loc[(abx_usage_icu['rank'] == 3), 'new_course'].sum(),
    abx_usage_icu['new_course'].sum(),
    (abx_usage_icu.loc[(
            (abx_usage_icu['new_course'] == 1)),
        'rank'] == 3).mean()*100]
text_n.loc[31] = [
    'Rank 4 courses',
    abx_usage_icu.loc[(abx_usage_icu['rank'] == 4), 'new_course'].sum(),
    abx_usage_icu['new_course'].sum(),
    (abx_usage_icu.loc[(
            (abx_usage_icu['new_course'] == 1)),
        'rank'] == 4).mean()*100]
kruskal_p = kruskal(
    descriptors.loc[descriptors['sepsis_wo_shock'], 'lengthofstay'],
    descriptors.loc[descriptors['septic_shock'], 'lengthofstay'],
    descriptors.loc[descriptors['antibiotics_wo_sepsis'], 'lengthofstay'],
    descriptors.loc[descriptors['no_antibiotics'], 'lengthofstay'])
los_v1 = descriptors.loc[(
    (descriptors['antibiotics_less_4d'] == 1) |
    (descriptors['discharge_less_4d'] == 1)), 'lengthofstay']
los_v2 = descriptors.loc[(
    (descriptors['antibiotics_min_4d_death'] == 1)), 'lengthofstay']
los_v3 = descriptors.loc[(
    (descriptors['antibiotics_less_4d'] == 1)), 'lengthofstay']
los_v4 = descriptors.loc[(
    (descriptors['antibiotics_min_4d_departure'] == 1)), 'lengthofstay']
text_n.loc[32] = [
    'LoS, abx <=4d or discharge <=4d (a)',
    ('%.0f' % los_v1.median() + ' (%.0f' % los_v1.quantile(0.25) +
        '-%.0f' % los_v1.quantile(0.75) + ')'),
    np.nan, np.nan]
text_n.loc[33] = [
    'LoS, abx >4d or until death (b)',
    ('%.0f' % los_v2.median() + ' (%.0f' % los_v2.quantile(0.25) +
        '-%.0f' % los_v2.quantile(0.75) + ')'),
    np.nan, np.nan]
text_n.loc[34] = [
    'LoS, abx <=4d (c)',
    ('%.0f' % los_v3.median() + ' (%.0f' % los_v3.quantile(0.25) +
        '-%.0f' % los_v3.quantile(0.75) + ')'),
    np.nan, np.nan]
text_n.loc[35] = [
    'LoS, abx <=4d or until death/discharge (d)',
    ('%.0f' % los_v4.median() + ' (%.0f' % los_v4.quantile(0.25) +
        '-%.0f' % los_v4.quantile(0.75) + ')'),
    np.nan, np.nan]
text_n.loc[36] = [
    'KW p-val for ICU LoS (test statistic)',
    '%.3f' % kruskal_p.pvalue + ' (' + '%.1f' % kruskal_p.statistic + ')',
    np.nan, np.nan]
mannwhitneyu_p = mannwhitneyu(los_v1, los_v2)
text_n.loc[37] = [
    'MWU p-val for LoS (a) vs LoS (b) (test statistics)',
    '%.3f' % mannwhitneyu_p.pvalue + ' (' +
    '%.1f' % mannwhitneyu_p.statistic + ')',
    np.nan, np.nan]
mannwhitneyu_p = mannwhitneyu(los_v3, los_v4)
text_n.loc[38] = [
    'MWU p-val for LoS (c) vs LoS (d) (test statistics)',
    '%.3f' % mannwhitneyu_p.pvalue + ' (' +
    '%.1f' % mannwhitneyu_p.statistic + ')',
    np.nan, np.nan]

text_n.to_csv(inputs.output_file_path + 'text_n.csv')

###############################################################################
# Hazard ratios

combined_diagnoses['traumatic_brain_injury'] = (
    (combined_diagnoses['trauma'] == 1) &
    (combined_diagnoses['diagnosis'].str.contains(
        r'Neurochirurgisch|Head|Hoofdtrauma|Multi')))
combined_diagnoses['shock'] = (
    (combined_diagnoses['cardio'] == 1) &
    (combined_diagnoses['diagnosis'].str.contains(r'shock|Shock')))
combined_diagnoses['heart_failure'] = (
    ((combined_diagnoses['cardio'] == 1) &
        (combined_diagnoses['diagnosis'].str.contains(r'CHF|cordis'))) |
    (combined_diagnoses['cardiacarrest'] == 1))
a_str = r'Hemorrhage/hematoma, intracranial|stroke|Intracerebraal haema|S.A.B'
combined_diagnoses['aneurysm'] = (
    ((combined_diagnoses['neuro'] == 1) &
        (combined_diagnoses['diagnosis'].str.contains(a_str))) |
    ((combined_diagnoses['neurosurg'] == 1) &
        (combined_diagnoses['diagnosis'].str.contains(
            r'Hemorrhage/hematoma-intracranial'))))
combined_diagnoses['respiratory_failure'] = (
    (combined_diagnoses['respfailure'] == 1) &
    (combined_diagnoses['diagnosis'].str.contains(
        r'Respiratoire insufficiëntie|Pneumoni|ARDS')))
combined_diagnoses['trauma_other'] = (
    (combined_diagnoses['trauma'] == 1) &
    (combined_diagnoses['traumatic_brain_injury'] == 0))

diagnoses_add = ['traumatic_brain_injury', 'shock', 'heart_failure']
diagnoses_add += ['aneurysm', 'respiratory_failure', 'trauma_other']
descriptors = pd.merge(
    descriptors, combined_diagnoses[['admissionid'] + diagnoses_add],
    on='admissionid', how='left')

cph = CoxPHFitter()
cox_columns = ['lengthofstay', 'overleden_IC', 'emergency_medical']
cox_columns += ['emergency_surgical', 'agegroup_18', 'agegroup_40']
cox_columns += ['agegroup_50', 'agegroup_60', 'agegroup_70', 'female']
cox_columns += ['sepsis_wo_shock', 'septic_shock']
cox_df = descriptors.loc[(descriptors['female'].isna() == 0), cox_columns]
cph.fit(cox_df, 'lengthofstay', 'overleden_IC')

hr_v1 = cph.summary[
    ['exp(coef)', 'exp(coef) lower 95%', 'exp(coef) upper 95%', 'p']]

cph = CoxPHFitter()
cox_columns = ['lengthofstay', 'overleden_IC', 'agegroup_18', 'agegroup_40']
cox_columns += ['agegroup_50', 'agegroup_60', 'agegroup_70', 'female']
cox_columns += ['sepsis_wo_shock', 'septic_shock']
cox_df = descriptors.loc[(descriptors['female'].isna() == 0), cox_columns]
cph.fit(cox_df, 'lengthofstay', 'overleden_IC')

hr_v2 = cph.summary[
    ['exp(coef)', 'exp(coef) lower 95%', 'exp(coef) upper 95%', 'p']]
hr = pd.merge(
    hr_v1, hr_v2, on='covariate', how='outer', suffixes=['_v1', '_v2'])

cph = CoxPHFitter()
cox_columns = ['lengthofstay', 'overleden_IC', 'agegroup_18', 'agegroup_40']
cox_columns += ['agegroup_50', 'agegroup_60', 'agegroup_70', 'female']
cox_columns += ['sepsis_wo_shock', 'septic_shock']
cox_df = descriptors.loc[(
    ((descriptors['emergency_surgical'] == 1) |
        (descriptors['emergency_medical'] == 1)) &
    (descriptors['female'].isna() == 0)), cox_columns]
cph.fit(cox_df, 'lengthofstay', 'overleden_IC')

hr_v3 = cph.summary[
    ['exp(coef)', 'exp(coef) lower 95%', 'exp(coef) upper 95%', 'p']]
hr = pd.merge(hr, hr_v3, on='covariate', how='outer')

cph = CoxPHFitter()
descriptors['sofa_4-7'] = (
    (descriptors['sofa_total_score'] >= 4) &
    (descriptors['sofa_total_score'] < 8))
descriptors['sofa_8-11'] = (
    (descriptors['sofa_total_score'] >= 8) &
    (descriptors['sofa_total_score'] < 12))
descriptors['sofa_12-15'] = (
    (descriptors['sofa_total_score'] >= 12) &
    (descriptors['sofa_total_score'] < 16))
descriptors['sofa_16-'] = (descriptors['sofa_total_score'] >= 16)
cox_columns = ['lengthofstay', 'overleden_IC', 'agegroup_18', 'agegroup_40']
cox_columns += ['agegroup_50', 'agegroup_60', 'agegroup_70', 'female']
cox_columns += ['sofa_4-7', 'sofa_8-11', 'sofa_12-15', 'sofa_16-']
cox_df = descriptors.loc[(
    (descriptors['sepsis_wo_shock'] == 1) &
    (descriptors['female'].isna() == 0)), cox_columns]
cph.fit(cox_df, 'lengthofstay', 'overleden_IC')

hr_v4 = cph.summary[
    ['exp(coef)', 'exp(coef) lower 95%', 'exp(coef) upper 95%', 'p']]
hr = pd.merge(hr, hr_v4, on='covariate', how='outer', suffixes=['_v3', '_v4'])

cph = CoxPHFitter()
cox_columns = ['lengthofstay', 'overleden_IC', 'agegroup_18', 'agegroup_40']
cox_columns += ['agegroup_50', 'agegroup_60', 'agegroup_70', 'female']
cox_columns += ['sofa_8-11', 'sofa_12-15', 'sofa_16-']
cox_df = descriptors.loc[(
    (descriptors['septic_shock'] == 1) &
    (descriptors['female'].isna() == 0)), cox_columns]
cph.fit(cox_df, 'lengthofstay', 'overleden_IC')

hr_v5 = cph.summary[
    ['exp(coef)', 'exp(coef) lower 95%', 'exp(coef) upper 95%', 'p']]
hr = pd.merge(hr, hr_v5, on='covariate', how='outer')

cph = CoxPHFitter()
cox_columns = ['lengthofstay', 'overleden_IC', 'agegroup_18', 'agegroup_40']
cox_columns += ['agegroup_50', 'agegroup_60', 'agegroup_70', 'female']
cox_columns += ['sofa_4-7', 'sofa_8-11', 'sofa_12-15', 'sofa_16-']
cox_df = descriptors.loc[(
    (descriptors['antibiotics_wo_sepsis'] == 1) &
    (descriptors['female'].isna() == 0)), cox_columns]
cph.fit(cox_df, 'lengthofstay', 'overleden_IC')

hr_v6 = cph.summary[
    ['exp(coef)', 'exp(coef) lower 95%', 'exp(coef) upper 95%', 'p']]
hr = pd.merge(hr, hr_v6, on='covariate', how='outer', suffixes=['_v5', '_v6'])

cph = CoxPHFitter()
cox_columns = ['lengthofstay', 'overleden_IC', 'agegroup_18', 'agegroup_40']
cox_columns += ['agegroup_50', 'agegroup_60', 'agegroup_70', 'female']
cox_columns += ['sofa_4-7', 'sofa_8-11', 'sofa_12-15', 'sofa_16-']
cox_df = descriptors.loc[(
    (descriptors['no_antibiotics'] == 1) &
    (descriptors['female'].isna() == 0)), cox_columns]
cph.fit(cox_df, 'lengthofstay', 'overleden_IC', step_size=0.1)

hr_v7 = cph.summary[
    ['exp(coef)', 'exp(coef) lower 95%', 'exp(coef) upper 95%', 'p']]
hr = pd.merge(hr, hr_v7, on='covariate', how='outer')

cph = CoxPHFitter()
cox_columns = ['lengthofstay', 'overleden_IC', 'agegroup_18', 'agegroup_40']
cox_columns += ['agegroup_50', 'agegroup_60', 'agegroup_70', 'female']
cox_columns += ['sofa_4-7', 'sofa_8-11', 'sofa_12-15', 'sofa_16-']
cox_df = descriptors.loc[(
    (descriptors['female'].isna() == 0)), cox_columns]
cph.fit(cox_df, 'lengthofstay', 'overleden_IC')

hr_v8 = cph.summary[
    ['exp(coef)', 'exp(coef) lower 95%', 'exp(coef) upper 95%', 'p']]
hr = pd.merge(hr, hr_v8, on='covariate', how='outer', suffixes=['_v7', '_v8'])

cph = CoxPHFitter()
cox_columns = ['lengthofstay', 'overleden_IC', 'agegroup_18', 'agegroup_40']
cox_columns += ['agegroup_50', 'agegroup_60', 'agegroup_70', 'female']
cox_columns += ['sepsis_wo_shock', 'septic_shock']
cox_columns += ['traumatic_brain_injury', 'heart_failure', 'shock']
cox_columns += ['aneurysm', 'respiratory_failure', 'trauma_other']
cox_df = descriptors.loc[(
    (descriptors['female'].isna() == 0)), cox_columns]
cph.fit(cox_df, 'lengthofstay', 'overleden_IC')

hr_v9 = cph.summary[
    ['exp(coef)', 'exp(coef) lower 95%', 'exp(coef) upper 95%', 'p']]
hr_v9.columns += '_v9'
hr = pd.merge(hr, hr_v9, on='covariate', how='outer')

hr_table = hr.round(3).astype(str)
hr_table['HR, including admission category (CI) (p-value)'] = (
    hr_table['exp(coef)_v1'] + ' (' + hr_table['exp(coef) lower 95%_v1'] +
    ', ' + hr_table['exp(coef) upper 95%_v1'] + ') (' + hr_table['p_v1'] + ')')
hr_table['HR, excluding admission category (CI) (p-value)'] = (
    hr_table['exp(coef)_v2'] + ' (' + hr_table['exp(coef) lower 95%_v2'] +
    ', ' + hr_table['exp(coef) upper 95%_v2'] + ') (' + hr_table['p_v2'] + ')')
hr_table['HR, emergency admissions only (CI) (p-value)'] = (
    hr_table['exp(coef)_v3'] + ' (' + hr_table['exp(coef) lower 95%_v3'] +
    ', ' + hr_table['exp(coef) upper 95%_v3'] + ') (' + hr_table['p_v3'] + ')')
hr_table['HR, sepsis without shock (CI) (p-value)'] = (
    hr_table['exp(coef)_v4'] + ' (' + hr_table['exp(coef) lower 95%_v4'] +
    ', ' + hr_table['exp(coef) upper 95%_v4'] + ') (' + hr_table['p_v4'] + ')')
hr_table['HR, septic shock (CI) (p-value)'] = (
    hr_table['exp(coef)_v5'] + ' (' + hr_table['exp(coef) lower 95%_v5'] +
    ', ' + hr_table['exp(coef) upper 95%_v5'] + ') (' + hr_table['p_v5'] + ')')
hr_table['HR, antibiotics without sepsis (CI) (p-value)'] = (
    hr_table['exp(coef)_v6'] + ' (' + hr_table['exp(coef) lower 95%_v6'] +
    ', ' + hr_table['exp(coef) upper 95%_v6'] + ') (' + hr_table['p_v6'] + ')')
hr_table['HR, no antibiotics (CI) (p-value)'] = (
    hr_table['exp(coef)_v7'] + ' (' + hr_table['exp(coef) lower 95%_v7'] +
    ', ' + hr_table['exp(coef) upper 95%_v7'] + ') (' + hr_table['p_v7'] + ')')
hr_table['HR, all admissions (CI) (p-value)'] = (
    hr_table['exp(coef)_v8'] + ' (' + hr_table['exp(coef) lower 95%_v8'] +
    ', ' + hr_table['exp(coef) upper 95%_v8'] + ') (' + hr_table['p_v8'] + ')')
hr_table['HR, other diagnoses (CI) (p-value)'] = (
    hr_table['exp(coef)_v9'] + ' (' + hr_table['exp(coef) lower 95%_v9'] +
    ', ' + hr_table['exp(coef) upper 95%_v9'] + ') (' + hr_table['p_v9'] + ')')

hr_table.iloc[:, -9:].to_csv(inputs.output_file_path + 'hr_table.csv')
