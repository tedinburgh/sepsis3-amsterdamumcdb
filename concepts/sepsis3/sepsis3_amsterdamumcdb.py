## Author: Tom Edinburgh
## v1: date 16/09/2021.

## This script is to calculate SOFA scores from the AmsterdamUMCdb database and
## write them to a csv file (that also contains other admission information).

## It's adapted from this AmsterdamUMCdb script, which calculates SOFA scores:
## https://github.com/AmsterdamUMC/AmsterdamUMCdb/blob/master/concepts/severityscores/sofa.ipynb)
## But this code using pandas, working directly from the AmsterdamUMCdb .csv
## data files, instead of via SQL queries (as in the GitHub above).
## These SQL queries are at the bottom of the script for reference (but aren't
## used in the script).

################################################################################

import pandas as pd
import numpy as np
import re

## Given AmsterdamUMCdb is largely in Dutch, this module may be useful. But it's
## not strictly necessary
try:
    from googletrans import Translator
except ImportError as e:
    print('Module googletrans does not exist')
    pass

## Make sure you have loaded this correctly, and are in the central directory of
## of AmsterdamUMCdb i.e. cd FILEPATH/AmsterdamUMCdb/
## Folders in this directory should include data, amsterdamumcdb, concepts etc.
import amsterdamumcdb as adb
dictionary = adb.get_dictionary()

################################################################################
## Load the data

## The idea of this script is that we end up with a binary variable indicating
## a sepsis episode for each day for each patient. These variables limit the
## window in which we care about, and can be set to None to include all data.
## However, doing this may end up killing your session due to memory costs.
MIN_TIME = -3 ## in days
MAX_TIME = 10 ## in days

## The directory 'data' should contain the base AmsterdamUMCdb .csv files.
## Access must be specifically requested from the database owners. The
## directory 'additional_files' will contain output .csv files from this script.
file_path = '../../data/'
additional_file_path = '../../data/additional_files/'

list_columns = ['admissionid', 'itemid', 'valueid', 'value']
list_columns += ['updatedat', 'measuredat', 'registeredby']
listitems = pd.read_csv(file_path + 'listitems.csv', usecols=list_columns)

drug_columns = ['admissionid', 'itemid', 'item', 'duration', 'rate', 'rateunit']
drug_columns += ['start', 'stop', 'dose', 'doserateperkg', 'doseunitid']
drug_columns += ['doserateunitid', 'ordercategoryid']
drugitems = pd.read_csv(file_path + 'drugitems.csv', usecols=drug_columns)

freetextitems = pd.read_csv(file_path + 'freetextitems.csv')

admissions_df = pd.read_csv(file_path + 'admissions.csv')

## From the script 'reason_for_admission.py', which must be run first
combined_diagnoses = pd.read_csv(
    additional_file_path + 'combined_diagnoses.csv')

## Out of the AmsterdamUMCdb data files, numerics is by far the biggest. So
## we have to be a bit more careful loading it in (by chunks).
numerics_cols = ['admissionid', 'itemid', 'item', 'value', 'unitid']
numerics_cols += ['measuredat', 'registeredby', 'islabresult', 'fluidout']
numerics_dtypes = ['int64', 'int64', 'str', 'float64', 'int64']
numerics_dtypes += ['int64', 'str', 'bool', 'float64']
numerics_dtypes = dict(zip(numerics_cols, numerics_dtypes))
numerics_csv = dict(
    encoding='latin-1', usecols=numerics_cols, dtype=numerics_dtypes,
    chunksize=10**6)

def numerics_read(itemid_list, start=None, end=None, admissions_df=None):
    ## Load numerics by chunk and retains rows with itemid in the given list.
    ## If either start/end are given as inputs to this function, then only
    ## rows with start < t < end are retained (where t refers to time since
    ## admission if admissions_df is given, else is just the given
    ## measurement time)
    numerics_list = []
    ii = 0
    if admissions_df is not None:
        admissions_df = admissions_df[['admissionid', 'admittedat']]
    with pd.read_csv(file_path + 'numericitems.csv', **numerics_csv) as reader:
        for chunk in reader:
            if ((ii % 100) == 0):
                print(ii)
            chunk = chunk.loc[chunk['itemid'].isin(itemid_list)]
            if admissions_df is not None:
                chunk = pd.merge(chunk, admissions_df, \
                    on='admissionid', how='left')
                chunk['diff_measuredat'] = \
                    chunk['measuredat'] - chunk['admittedat']
            else:
                chunk['diff_measuredat'] = chunk['measuredat']
            if start is not None:
                chunk = chunk.loc[chunk['diff_measuredat'] >= start]
            if end is not None:
                chunk = chunk.loc[chunk['diff_measuredat'] <= end]
            numerics_list.append(chunk)
            ii += 1
    numerics = pd.concat(numerics_list)
    ## This is the 'day' since admission, i.e. 24hr period centred at admission
    ## time.
    numerics['time'] = numerics['diff_measuredat'] // (1000*60*60*24)
    numerics.drop(columns=['diff_measuredat'], inplace=True)
    numerics.reset_index(drop=True, inplace=True)
    return numerics

## This list is taken from the AmsterdamUMCdb SOFA scores script.
numerics_sofa_itemid = [8845] # O2 l/min
numerics_sofa_itemid += [10387] # Zuurstof toediening (bloed)
numerics_sofa_itemid += [18587] # Zuurstof toediening
numerics_sofa_itemid += [6699] # FiO2 %: setting on Evita ventilator
## 12279 # O2 concentratie --measurement by Servo-i/Servo-U ventilator
numerics_sofa_itemid += [12279]
numerics_sofa_itemid += [12369] # SET %O2: used with BiPap Vision ventilator
numerics_sofa_itemid += [16246] # Zephyros FiO2: Non-invasive ventilation
numerics_sofa_itemid += [8794] # UrineCAD
numerics_sofa_itemid += [8796] # UrineSupraPubis
numerics_sofa_itemid += [8798] # UrineSpontaan
numerics_sofa_itemid += [8800] # UrineIncontinentie
numerics_sofa_itemid += [8803] # UrineUP
numerics_sofa_itemid += [10743] # Nefrodrain li Uit
numerics_sofa_itemid += [10745] # Nefrodrain re Uit
numerics_sofa_itemid += [19921] # UrineSplint Li
numerics_sofa_itemid += [19922] # UrineSplint Re]
numerics_sofa_itemid += [6846] # PCO2
numerics_sofa_itemid += [9990] # pCO2 (bloed)
numerics_sofa_itemid += [21213] # PCO2 (bloed) - kPa
numerics_sofa_itemid += [7433] # PO2
numerics_sofa_itemid += [9996] # PO2 (bloed)
numerics_sofa_itemid += [21214] # PO2 (bloed) - kPa
numerics_sofa_itemid += [9964] # Thrombo's (bloed)
numerics_sofa_itemid += [6797] # Thrombocyten
numerics_sofa_itemid += [10409] # Thrombo's citr. bloed (bloed)
numerics_sofa_itemid += [14252] # Thrombo CD61 (bloed)
numerics_sofa_itemid += [6813] # Bili Totaal
numerics_sofa_itemid += [9945] # Bilirubine (bloed)

numerics_abp_itemid = [6642] # ABP gemiddeld
numerics_abp_itemid += [6679] # Niet invasieve bloeddruk gemiddeld
numerics_abp_itemid += [8843] # ABP gemiddeld II

## 6836: Kreatinine µmol/l (erroneously documented as µmol)
numerics_creatinine_itemid = [6836]
numerics_creatinine_itemid += [9941] # Kreatinine (bloed) µmol/l
numerics_creatinine_itemid += [14216] # KREAT enzym. (bloed) µmol/l

end_time = 1000*60*60*24*MAX_TIME if MAX_TIME is not None else None
start_time = 1000*60*60*24*MIN_TIME if MIN_TIME is not None else None
numerics_sofa = numerics_read(
    numerics_sofa_itemid + numerics_abp_itemid, admissions_df=admissions_df,
    end=end_time, start=start_time)
## We need a baseline creatinine, so look back further.
numerics_creatinine = numerics_read(
    numerics_creatinine_itemid, admissions_df=admissions_df,
    end=end_time, start=-1000*60*60*24*365)

################################################################################
## Add admission time etc. to dataframes, and crop drug/list etc. tables to
## min/max times

admissions_add = admissions_df[['admissionid', 'admittedat', 'dischargedat']]
listitems = pd.merge(listitems, admissions_add, on='admissionid', how='left')
listitems['time'] = (listitems['measuredat'] - listitems['admittedat'])
listitems['time'] //= (1000*60*60*24) ## Convert to 'day'

drugitems = pd.merge(drugitems, admissions_add, on='admissionid', how='left')
drugitems['start_time'] = (drugitems['start'] - drugitems['admittedat'])
drugitems['stop_time'] = (drugitems['stop'] - drugitems['admittedat'])
drugitems[['start_time', 'stop_time']] //= (1000*60*60*24) ## Convert to 'day'

freetextitems = pd.merge(
    freetextitems, admissions_add, on='admissionid', how='left')
freetextitems['time'] = (
    freetextitems['measuredat'] - freetextitems['admittedat'])
freetextitems['time'] //= (1000*60*60*24)

## Cut only to the time window specified earlier
if MIN_TIME is not None:
    listitems = listitems.loc[(listitems['time'] >= MIN_TIME)]
    drugitems = drugitems.loc[(drugitems['start_time'] >= MIN_TIME)]
    freetextitems = freetextitems.loc[(freetextitems['time'] >= MIN_TIME)]
if MAX_TIME is not None:
    listitems = listitems.loc[(listitems['time'] <= MAX_TIME)]
    ## Keep this as start time rather than end time (i.e. started drug treatment
    ## on that particular 'day')
    drugitems = drugitems.loc[(drugitems['start_time'] <= MAX_TIME)]
    freetextitems = freetextitems.loc[(freetextitems['time'] <= MAX_TIME)]

def systeem_flag(df):
    return df['registeredby'].str.contains('systeem', flags=re.IGNORECASE)

################################################################################
## Begin the SOFA scores (again following AmsterdamUMCdb!)

########################
## Respiration score

## Get PaO2/FiO2 ratio
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
## Get PaO2 and FiO2 values
## Simultaneously retrieve PaCO2 and the 'nearest' FiO2 from the ventilator or
## estimated FiO2 based on applied oxygen device. Ideally documentation of
## measurements should be at the same time, but since this is not guaranteed
## allow a window.
## In more recent data PaCO2 and PaO2 were documented in kPa instead of mmHg.
fio2_itemids = [8845, 10387, 18587, 6699, 12279, 12369, 16246]
fio2_table = numerics_sofa.loc[
    (numerics_sofa['value'] > 0) &
    (numerics_sofa['itemid'].isin(fio2_itemids))]
fio2_table = pd.merge(
    fio2_table, oxy_flow_listitems.drop(columns=['time']),
    on=['admissionid', 'measuredat'], how='left')
fio2_table.rename(columns={'value_y': 'O2_device'}, inplace=True)
fio2_table.rename(columns={'value_x': 'value'}, inplace=True)
fio2_table['ventilatory_support'] = False
fio2_table.loc[
        fio2_table['itemid'].isin([6699, 12279, 12369, 16246]),
    'ventilatory_support'] = True

fio2_table['fio2'] = 0.21
fio2_ind = fio2_table['ventilatory_support'] & (~fio2_table['value'].isnull())
fio2_table.loc[fio2_ind, 'fio2'] = fio2_table.loc[fio2_ind, 'value']

valueids1 = [1, 2, 3, 4, 7, 8, 9, 18, 19]
## Diep Nasaal, Nasaal, Kapje, Kunstneus, O2-bril, Kinnebak, Nebulizer,
## Spreekcanule, Spreekklepje
temp_ind = (
    (~fio2_table['ventilatory_support']) &
    (fio2_table['valueid'].isin(valueids1)))
fio2_table.loc[
        temp_ind & (fio2_table['value'] >= 1) & (fio2_table['value'] < 2),
    'fio2'] = 0.22
fio2_table.loc[
        temp_ind & (fio2_table['value'] >= 2) & (fio2_table['value'] < 3),
    'fio2'] = 0.25
fio2_table.loc[
        temp_ind & (fio2_table['value'] >= 3) & (fio2_table['value'] < 4),
    'fio2'] = 0.27
fio2_table.loc[
        temp_ind & (fio2_table['value'] >= 4) & (fio2_table['value'] < 5),
    'fio2'] = 0.30
fio2_table.loc[temp_ind & (fio2_table['value'] >= 5), 'fio2'] = 0.35

valueids2 = [1, 3, 4, 8, 9, 18, 19]
## Diep Nasaal, Kapje, Kunstneus, Kinnebak, Nebulizer, Spreekcanule,
## Spreekklepje
temp_ind = (
    (~fio2_table['ventilatory_support']) &
    (fio2_table['valueid'].isin(valueids2)))
fio2_table.loc[
    temp_ind & (fio2_table['value'] >= 6) & (fio2_table['value'] < 7),
    'fio2'] = 0.40
fio2_table.loc[
        temp_ind & (fio2_table['value'] >= 7) & (fio2_table['value'] < 8),
    'fio2'] = 0.45
fio2_table.loc[temp_ind & (fio2_table['value'] >= 8), 'fio2'] = 0.50

valueids3 = [10, 11, 13, 14, 15, 16, 17]
## Waterset, Trach.stoma, Ambu, Guedel, DL-tube, CPAP, Non-Rebreathing masker
temp_ind = (
    (~fio2_table['ventilatory_support']) &
    (fio2_table['valueid'].isin(valueids3)))
fio2_table.loc[
        temp_ind & (fio2_table['value'] >= 6) & (fio2_table['value'] < 7),
    'fio2'] = 0.60
fio2_table.loc[
        temp_ind & (fio2_table['value'] >= 7) & (fio2_table['value'] < 8),
    'fio2'] = 0.70
fio2_table.loc[
        temp_ind & (fio2_table['value'] >= 8) & (fio2_table['value'] < 9),
    'fio2'] = 0.80
fio2_table.loc[
        temp_ind & (fio2_table['value'] >= 9) & (fio2_table['value'] < 10),
    'fio2'] = 0.85
fio2_table.loc[temp_ind & (fio2_table['value'] >= 10), 'fio2'] = 0.90
fio2_table.rename(columns={'measuredat': 'fio2_measuredat'}, inplace=True)
fio2_cols = ['admissionid', 'fio2_measuredat', 'fio2']
fio2_cols += ['ventilatory_support', 'time']

## This is initially pao2 and then merged with fio2 from above.
oxygenation = numerics_sofa.loc[
    numerics_sofa['itemid'].isin([7433, 9996, 21214])]
## itemids above are: PO2, PO2 (bloed), PO2 (bloed) - kPa
## Conversion from kPa to mmHg
oxygenation.loc[oxygenation['unitid'] == 152, 'value'] *= 7.50061683
oxygenation.rename(columns={'value': 'pao2'}, inplace=True)
oxygenation['manual_entry'] = True
oxygenation.loc[systeem_flag(oxygenation), 'manual_entry'] = False

f = freetextitems.loc[freetextitems['itemid'] == 11646]
f = f[['admissionid', 'measuredat', 'value']]
f.rename(columns={'value': 'specimen_source'}, inplace=True)
oxygenation = pd.merge(
    oxygenation, f, on=['admissionid','measuredat'], how='left')

oxygenation = oxygenation.loc[
    oxygenation['specimen_source'].isnull() |
    oxygenation['specimen_source'].str.contains('art', flags=re.IGNORECASE)]

oxygenation = pd.merge(
    oxygenation, fio2_table[fio2_cols], on=['admissionid','time'], how='left')
oxygenation['FiO2_time_difference'] = (
    oxygenation['fio2_measuredat'] - oxygenation['measuredat'])
## Keep fio2 only if no earlier than 60 minutes before pao2 measurement
oxygenation = oxygenation.loc[oxygenation['FiO2_time_difference'] > -1000*60*60]
## and no later than 15 minutes after pao2 measuremen
oxygenation = oxygenation.loc[oxygenation['FiO2_time_difference'] < 1000*60*15]
## Convert to days (not discrete)
oxygenation['FiO2_time_difference'] /= (1000*60*60*24)
oxygenation['abs_measuredat'] = oxygenation['FiO2_time_difference'].abs()

## Sort by the smallest fio2 time difference for each patient and timestamp
oxygenation = oxygenation.sort_values(
    by=['admissionid','measuredat','abs_measuredat'])
## Discard duplicates of the same patient and timestamp (keeping only the
## smallest fio2 time difference)
oxygenation = oxygenation.loc[
    (~(oxygenation[['admissionid', 'measuredat']].duplicated()))]
oxygenation['priority'] = 1

sofa_resp_cols = ['admissionid', 'pao2', 'specimen_source', 'manual_entry']
sofa_resp_cols += ['time', 'fio2', 'ventilatory_support']
sofa_resp_cols += ['FiO2_time_difference', 'priority']
sofa_respiration = oxygenation[sofa_resp_cols]
sofa_respiration.head()

## Remove extreme outliers (in the AmsterdamUMCdb script, histograms are
## plotted to identify these outliers by eye, here we just copy those values)
sofa_respiration.loc[(sofa_respiration['fio2'] > 100), 'fio2'] = np.nan
## Convert FiO2 in % to fraction
sofa_respiration.loc[
        (sofa_respiration['fio2'] <= 100) &
        (sofa_respiration['fio2'] >= 20),
    'fio2'] /= 100
## Remove lower outliers, most likely incorrectly labeled as 'arterial'
## instead of '(mixed/central) venous'
sofa_respiration.loc[sofa_respiration['pao2'] < 50, 'pao2'] = np.nan
sofa_respiration = sofa_respiration.dropna(subset=['pao2'])

## Get the PF ratio
sofa_respiration.loc[:,'pf_ratio'] = (
    sofa_respiration['pao2'] / sofa_respiration['fio2'])
sofa_respiration['ventilatory_support'] = (
    sofa_respiration['ventilatory_support'] == True)

## Calculate SOFA respiration score:
sofa_respiration['sofa_respiration_score'] = 0
sofa_respiration.loc[
        (sofa_respiration['pf_ratio'] < 400) &
        (sofa_respiration['pf_ratio'] >= 300),
    'sofa_respiration_score'] = 1
sofa_respiration.loc[
        (sofa_respiration['pf_ratio'] < 300),
    'sofa_respiration_score'] = 2
sofa_respiration.loc[
        (sofa_respiration['pf_ratio'] < 200) &
        (sofa_respiration['pf_ratio'] >= 100) &
        sofa_respiration['ventilatory_support'],
    'sofa_respiration_score'] = 3
sofa_respiration.loc[
        (sofa_respiration['pf_ratio'] < 100) &
        sofa_respiration['ventilatory_support'],
    'sofa_respiration_score'] = 4

sofa_respiration.head()

########################
## Coagulation score (platelets (thrombocytes))

sofa_platelets = numerics_sofa.loc[
    numerics_sofa['itemid'].isin([9964, 6797, 10409, 14252])]
## itemids above are: Thrombo's (bloed), Thrombocyten, Thrombo's citr. bloed
## (bloed), Thrombo CD61 (bloed)
sofa_platelets['manual_entry'] = True
sofa_platelets.loc[systeem_flag(sofa_platelets), 'manual_entry'] = False
sofa_platelets_columns = ['admissionid', 'itemid', 'item', 'value']
sofa_platelets_columns += ['registeredby', 'manual_entry', 'time']
sofa_platelets = sofa_platelets[sofa_platelets_columns]

## Calculate SOFA coagulation score:
sofa_platelets['sofa_coagulation_score'] = 0
sofa_platelets.loc[
        (sofa_platelets['value'] < 150) & (sofa_platelets['value'] >= 100),
    'sofa_coagulation_score'] = 1
sofa_platelets.loc[
        (sofa_platelets['value'] < 100) & (sofa_platelets['value'] >= 50),
    'sofa_coagulation_score'] = 3
sofa_platelets.loc[
        (sofa_platelets['value'] < 50) & (sofa_platelets['value'] >= 20),
    'sofa_coagulation_score'] = 3
sofa_platelets.loc[
        (sofa_platelets['value'] < 20),
    'sofa_coagulation_score'] = 4

sofa_platelets.head()

########################
## Liver score (bilirubin)

sofa_bilirubin = numerics_sofa.loc[numerics_sofa['itemid'].isin([6813, 9945])]
sofa_bilirubin['manual_entry'] = True
sofa_bilirubin.loc[systeem_flag(sofa_bilirubin), 'manual_entry'] = False
sofa_bilirubin_columns = ['admissionid', 'itemid', 'item', 'value']
sofa_bilirubin_columns += ['registeredby', 'manual_entry', 'time']
sofa_bilirubin = sofa_bilirubin[sofa_bilirubin_columns]

## Calculate SOFA liver score:
sofa_bilirubin['sofa_liver_score'] = 0
sofa_bilirubin.loc[
        (sofa_bilirubin['value'] >= 20) & (sofa_bilirubin['value'] < 33),
    'sofa_liver_score'] = 1
sofa_bilirubin.loc[
        (sofa_bilirubin['value'] >= 33) & (sofa_bilirubin['value'] < 102),
    'sofa_liver_score'] = 2
sofa_bilirubin.loc[
        (sofa_bilirubin['value'] >= 102) & (sofa_bilirubin['value'] < 204),
    'sofa_liver_score'] = 3
sofa_bilirubin.loc[
        (sofa_bilirubin['value'] >= 204),
    'sofa_liver_score'] = 4

sofa_bilirubin.head()

########################
## Cardiovascular score

cv_drugitems = [7179] # Dopamine (Inotropin)
cv_drugitems += [7178] # Dobutamine (Dobutrex)
cv_drugitems += [6818] # Adrenaline (Epinefrine)
cv_drugitems += [7229] # Noradrenaline (Norepinefrine)

sofa_cardiovascular = drugitems.loc[
    (drugitems['itemid'].isin(cv_drugitems)) &
    (drugitems['rate'] > 0.1) &
    (drugitems['ordercategoryid'] == 65)]

admissions_add = admissions_df[['admissionid', 'admittedat', 'weightgroup']]
sofa_cardiovascular = pd.merge(
    sofa_cardiovascular, admissions_add, on='admissionid', how='left')

weight_group_dict = {
    '59-': 55, '60-69': 65, '70-79': 75, '80-89': 85, '90-99': 95,
    '100-109': 105, '110+': 115, np.nan: 80}
sofa_cardiovascular['patientweight'] = (
    sofa_cardiovascular['weightgroup'].replace(weight_group_dict))

## Want to add extra rows to the dataframe, for when drug administration
## happened over consecutive 'days' (i.e. make an entry for each 'day' that
## the drug administration window overlaps with)
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

## Calculate gamma, as per AmsterdamUMCdb script
sofa_cardiovascular['gamma'] = (sofa_cardiovascular['dose'] / 80)
valid_weight_ind = sofa_cardiovascular['patientweight'] > 0
sofa_cardiovascular.loc[valid_weight_ind, 'gamma'] = (
    sofa_cardiovascular.loc[valid_weight_ind, 'dose'] /
    sofa_cardiovascular.loc[valid_weight_ind, 'patientweight'])
sofa_cardiovascular.loc[sofa_cardiovascular['doserateperkg'] == 1, 'gamma'] = (
    sofa_cardiovascular.loc[sofa_cardiovascular['doserateperkg'] == 1, 'dose'])
sofa_cardiovascular.loc[
        sofa_cardiovascular['doseunitid'] == 10,
    'gamma'] *= 1000
sofa_cardiovascular.loc[
        sofa_cardiovascular['doserateunitid'] == 5,
    'gamma'] /= 60

## Mean ABP
mean_abp = numerics_sofa.loc[numerics_sofa['itemid'].isin(numerics_abp_itemid)]
mean_abp['validated'] = True
mean_abp.loc[mean_abp['registeredby'].isnull(), 'validated'] = False
mean_abp.head()

## Remove extreme outliers, most likely data entry errors or measurement errors
mean_abp.loc[(mean_abp['value'] > 165), 'value'] = np.nan
mean_abp.loc[(mean_abp['value'] <= 30), 'value'] = np.nan
mean_abp_cols = ['admissionid', 'itemid', 'item', 'value', 'validated', 'time']
mean_abp = mean_abp[mean_abp_cols]
mean_abp = mean_abp.dropna(subset=['value'])
## Use mean_abp 'cleansed' dataframe
cv_groupby_cols = ['admissionid', 'itemid', 'item', 'time']
sofa_cardiovascular_map = mean_abp.groupby(cv_groupby_cols).agg(
        lowest_mean_abp=pd.NamedAgg(column='value', aggfunc='min')
    ).reset_index()

## Calculate SOFA cardiovascular score:
sofa_cardiovascular_map['sofa_cardiovascular_score'] = 0
## MAP < 70
sofa_cardiovascular_map.loc[
        (sofa_cardiovascular_map['lowest_mean_abp'] < 70),
    'sofa_cardiovascular_score'] = 1
sofa_cardiovascular_map.head()

sofa_cardiovascular_meds = sofa_cardiovascular.groupby(cv_groupby_cols).agg(
        total_duration=pd.NamedAgg(column='duration', aggfunc='sum'),
        max_gamma=pd.NamedAgg(column='gamma', aggfunc='max')
    ).reset_index()
sofa_cardiovascular_meds.head()

sofa_cardiovascular_meds['sofa_cardiovascular_score'] = 0
## Dopamine (itemid 7179) <= 5 or dobutamine (itemid 7178) any dose
sofa_cardiovascular_meds.loc[
        ((sofa_cardiovascular_meds['itemid'] == 7179) &
            (sofa_cardiovascular_meds['max_gamma'] <= 5)) |
        (sofa_cardiovascular_meds['itemid'] == 7178),
    'sofa_cardiovascular_score'] = 2
## Dopamine (itemid 7179) > 5, epinephrine (itemid 6818) <= 0.1,
## norepinephrine (itemid 7229) <= 0.1
sofa_cardiovascular_meds.loc[
        ((sofa_cardiovascular_meds['itemid'] == 7179) &
            (sofa_cardiovascular_meds['max_gamma'] > 5) &
            (sofa_cardiovascular_meds['max_gamma'] < 15)) |
        ((sofa_cardiovascular_meds['itemid'] == 6818) &
            (sofa_cardiovascular_meds['max_gamma'] <= 0.1)) |
        ((sofa_cardiovascular_meds['itemid'] == 7229) &
            (sofa_cardiovascular_meds['max_gamma'] <= 0.1)),
    'sofa_cardiovascular_score'] = 3
## Dopamine (itemid 7179) > 15, epinephrine (itemid 6818) > 0.1,
## norepinephrine (itemid 7229) > 0.1
sofa_cardiovascular_meds.loc[
        ((sofa_cardiovascular_meds['itemid'] == 7179) &
            (sofa_cardiovascular_meds['max_gamma'] > 15)) |
        ((sofa_cardiovascular_meds['itemid'] == 6818) &
            (sofa_cardiovascular_meds['max_gamma'] > 0.1)) |
        ((sofa_cardiovascular_meds['itemid'] == 7229) &
            (sofa_cardiovascular_meds['max_gamma'] > 0.1)),
    'sofa_cardiovascular_score'] = 4
sofa_cardiovascular_meds.head()

## Combine the scores from MAP and cardiovascular medication
sofa_cardiovascular = pd.concat(
        [sofa_cardiovascular_map, sofa_cardiovascular_meds], sort=False,
    ).sort_values(by=['admissionid','time']).reset_index(drop=True)

sofa_cardiovascular.head()

########################
## Glasgow Coma Scale score

eyes_itemids = [6732] # Actief openen van de ogen
eyes_itemids += [13077] # A_Eye
eyes_itemids += [14470] # RA_Eye
eyes_itemids += [16628] # MCA_Eye
eyes_itemids += [19635] # E_EMV_NICE_24uur
eyes_itemids += [19638] # E_EMV_NICE_Opname
motor_itemids = [6734] # Beste motore reactie van de armen
motor_itemids += [13072] # A_Motoriek
motor_itemids += [14476] # RA_Motoriek
motor_itemids += [16634] # MCA_Motoriek
motor_itemids += [19636] # M_EMV_NICE_24uur
motor_itemids += [19639] # M_EMV_NICE_Opname
verbal_itemids = [6735] # Beste verbale reactie
verbal_itemids += [13066] # A_Verbal
verbal_itemids += [14482] # RA_Verbal
verbal_itemids += [16640] # MCA_Verbal
verbal_itemids += [19637] # V_EMV_NICE_24uur
verbal_itemids += [19640] # V_EMV_NICE_Opname

## GCS eyes component
gcs_components = listitems.loc[listitems['itemid'].isin(eyes_itemids)]
gcs_components['eyes_score'] = 0
## Actief openen van de ogen
gcs_components.loc[gcs_components['itemid'] == 6732, 'eyes_score'] = (
    5 - gcs_components.loc[gcs_components['itemid'] == 6732, 'valueid'])
## A_Eye
gcs_components.loc[gcs_components['itemid'] == 13077, 'eyes_score'] = (
    gcs_components.loc[gcs_components['itemid'] == 13077, 'valueid'])
## RA_Eye, MCA_Eye, E_EMV_NICE_24uur
gcs_components.loc[
        gcs_components['itemid'].isin([14470, 16628, 19635]), 'eyes_score'] = (
    gcs_components.loc[
        gcs_components['itemid'].isin([14470, 16628, 19635]), 'valueid'] - 4)
## E_EMV_NICE_Opname
gcs_components.loc[gcs_components['itemid'] == 19638, 'eyes_score'] = (
    gcs_components.loc[gcs_components['itemid'] == 19638, 'valueid'] - 8)

## Preference, ranked by discipline
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
## Only keep the lowest score for the discipline of smallest rank
## (i.e. ICV_Medisch Staflid supercedes ICV_Medisch etc.)
gcs_components = gcs_components.loc[
    ~gcs_components[['admissionid', 'time']].duplicated()]
gcs_components.drop(
    columns=['itemid', 'valueid', 'value', 'admittedat', 'registeredby'],
    inplace=True)

gcs_cols = ['admissionid', 'measuredat', 'itemid', 'valueid', 'registeredby']
## Add GCS motor score
gcs_components = pd.merge(
    gcs_components,
    listitems.loc[listitems['itemid'].isin(motor_itemids), gcs_cols],
    on=['admissionid', 'measuredat'],
    how='left')
gcs_components['motor_score'] = 0
## Beste motore reactie van de armen
gcs_components.loc[gcs_components['itemid'] == 6734, 'motor_score'] = (
    7 - gcs_components.loc[gcs_components['itemid'] == 6734, 'valueid'])
## A_Motoriek
gcs_components.loc[gcs_components['itemid'] == 13072, 'motor_score'] = (
    gcs_components.loc[gcs_components['itemid'] == 13072, 'valueid'])
## RA_Motoriek, MCA_Motoriek, M_EMV_NICE_24uur
gcs_components.loc[
        gcs_components['itemid'].isin([14476, 16634, 19636]), 'motor_score'] = (
    gcs_components.loc[
        gcs_components['itemid'].isin([14476, 16634, 19636]), 'valueid'] - 6)
## M_EMV_NICE_Opname
gcs_components.loc[gcs_components['itemid'] == 19639, 'motor_score'] = (
    gcs_components.loc[gcs_components['itemid'] == 19639, 'valueid'] - 12)

## As above, add in preference by discipline
gcs_components['preference'] = 8
gcs_ind = gcs_components['registeredby'].isin(gcs_preferences.keys())
gcs_components.loc[gcs_ind, 'preference'] = (
    gcs_components.loc[gcs_ind, 'registeredby'].replace(gcs_preferences))
## Give higher preference by discipline (eye score should be the same for
## each admission id+time here)
gcs_components.sort_values(
    by=['admissionid', 'time', 'preference', 'motor_score'], inplace=True)
## Only keep the highest score for the discipline of smallest rank
gcs_components = gcs_components.loc[
    ~gcs_components[['admissionid', 'time']].duplicated()]
gcs_components.drop(columns=['itemid', 'valueid', 'registeredby'], inplace=True)
## Motor score is a float (due to pandas merge, so convert to int)
gcs_components['motor_score'] = gcs_components['motor_score'].astype(int)

## Add GCS verbal score
gcs_components = pd.merge(
    gcs_components,
    listitems.loc[listitems['itemid'].isin(verbal_itemids), gcs_cols],
    on=['admissionid', 'measuredat'],
    how='left')
gcs_components['verbal_score'] = 0
## Beste verbale reactie
gcs_components.loc[gcs_components['itemid'] == 6735, 'verbal_score'] = (
    6 - gcs_components.loc[gcs_components['itemid'] == 6735, 'valueid'])
## A_Verbal
gcs_components.loc[gcs_components['itemid'] == 13066, 'verbal_score'] = (
    gcs_components.loc[gcs_components['itemid'] == 13066, 'valueid'])
## RA_Verbal, MCA_Verbal
gcs_components.loc[
        gcs_components['itemid'].isin([14482, 16640]), 'verbal_score'] = (
    gcs_components.loc[
        gcs_components['itemid'].isin([14482, 16640]), 'valueid'] - 5)
## V_EMV_NICE_24uur
gcs_components.loc[gcs_components['itemid'] == 19637, 'verbal_score'] = (
    gcs_components.loc[gcs_components['itemid'] == 19637, 'valueid'] - 9)
## V_EMV_NICE_Opname
gcs_components.loc[gcs_components['itemid'] == 19640, 'verbal_score'] = (
    gcs_components.loc[gcs_components['itemid'] == 19640, 'valueid'] - 15)

## As above, add in preference by discipline
gcs_components['preference'] = 8
gcs_ind = gcs_components['registeredby'].isin(gcs_preferences.keys())
gcs_components.loc[gcs_ind, 'preference'] = (
    gcs_components.loc[gcs_ind, 'registeredby'].replace(gcs_preferences))
## Give higher preference by discipline
gcs_components.sort_values(
    by=['admissionid', 'time', 'preference', 'verbal_score'], inplace=True)
## Only keep the highest score for the discipline of smallest rank
gcs_components = gcs_components.loc[
    ~gcs_components[['admissionid', 'time']].duplicated()]
gcs_components.drop(columns=['itemid', 'valueid', 'registeredby'], inplace=True)
## Verbal score minimum is 1
gcs_components.loc[gcs_components['verbal_score'] < 1, 'verbal_score'] = 1
gcs_components['verbal_score'] = gcs_components['verbal_score'].astype(int)

## Combine the scores for total GCS score
gcs_components['min_gcs'] = (gcs_components['eyes_score'] +
    gcs_components['motor_score'] + gcs_components['verbal_score'])
gcs_components.head()
sofa_cns = gcs_components[['admissionid', 'time', 'min_gcs']]

## Calculate SOFA cardiovascular score:
sofa_cns['sofa_cns_score'] = 0
## MAP < 70
sofa_cns.loc[
        (sofa_cns['min_gcs'] >= 13) & (sofa_cns['min_gcs'] < 15),
    'sofa_cns_score'] = 1
sofa_cns.loc[
        (sofa_cns['min_gcs'] >= 10) & (sofa_cns['min_gcs'] < 13),
    'sofa_cns_score'] = 2
sofa_cns.loc[
        (sofa_cns['min_gcs'] >= 6) & (sofa_cns['min_gcs'] < 10),
    'sofa_cns_score'] = 3
sofa_cns.loc[(sofa_cns['min_gcs'] < 6), 'sofa_cns_score'] = 4

sofa_cns.head()

########################
## Renal score

## Get urineoutput
sofa_ruo_itemid = [8794, 8796, 8798, 8800, 8803, 10743, 10745, 19921, 19922]
## Dataframe is called sofa_renal_urine_output in amsterdamUMCdb sofa script
sofa_urine_out = numerics_sofa.loc[
    numerics_sofa['itemid'].isin(sofa_ruo_itemid)]
sofa_urine_out.drop(
    columns=['unitid', 'measuredat', 'islabresult', 'fluidout', 'admittedat'],
    inplace=True)
sofa_urine_out.head()

## Probably decimal error when entering volumes > 2500
sofa_urine_out.loc[(sofa_urine_out['value'] > 2500), 'value'] /= 10
## Remove extreme outliers, most likely data entry error)
sofa_urine_out.loc[(sofa_urine_out['value'] > 4500), 'value'] = np.nan
sofa_urine_out = sofa_urine_out.dropna(subset=['value'])
## Get urine output per 24 hours
sofa_daily_urine_out = sofa_urine_out.groupby(['admissionid','time']).agg(
        daily_urine_output=pd.NamedAgg(column='value', aggfunc='sum')
    ).reset_index()
sofa_daily_urine_out.head()

## Calculate SOFA renal score for urine output:
sofa_daily_urine_out['sofa_renal_score'] = 0
## Urine output < 500 ml/day
sofa_daily_urine_out.loc[
        (sofa_daily_urine_out['daily_urine_output'] < 500) &
        (sofa_daily_urine_out['daily_urine_output'] > 200),
     'sofa_renal_score'] = 3
## Urine output < 200 ml/day
sofa_daily_urine_out.loc[
        (sofa_daily_urine_out['daily_urine_output'] < 200),
    'sofa_renal_score'] = 4
sofa_daily_urine_out.head()

## Get serum creatinine (baseline from -365 days from admission)
baseline_creatinine = numerics_creatinine.groupby(['admissionid']).agg(
        baseline_creatinine=pd.NamedAgg(column='value', aggfunc='min')
    ).reset_index()
## Max creatinine on each day (but only from MIN_TIME rather than -365 days)
max_creatinine = numerics_creatinine.loc[
    numerics_creatinine['time'] >= MIN_TIME]
max_creatinine = max_creatinine.groupby(['admissionid','time']).agg(
        max_creatinine=pd.NamedAgg(column='value', aggfunc='max')
    ).reset_index()
## Merge baseline on admissionid only and max on both admissionid and time
creatinine = pd.merge(
    numerics_creatinine, baseline_creatinine, on='admissionid', how='left')
creatinine = pd.merge(
    creatinine, max_creatinine, on=['admissionid','time'], how='right')

creatinine['manual_entry'] = True
creatinine.loc[systeem_flag(creatinine), 'manual_entry'] = False

creatinine['acute_renal_failure'] = False
## AKI definition: 3 fold increase
creatinine.loc[
        (creatinine['baseline_creatinine'] > 0) &
        (creatinine['max_creatinine'] / creatinine['baseline_creatinine'] > 3),
    'acute_renal_failure'] = True
## AKI definition: increase to >= 354 umol/l AND at least 44 umol/l increase
creatinine.loc[
        (creatinine['max_creatinine'] >= 354) &
        ((creatinine['max_creatinine'] -
            creatinine['baseline_creatinine']) >= 44),
    'acute_renal_failure'] = True

creatinine.drop(
    columns=['unitid', 'measuredat', 'islabresult', 'fluidout', 'admittedat'],
    inplace=True)

creatinine.head()

## Looking at the data it's relevatively easy to spot most lab collection errors
## (i.e. single outliers between relatively normal values)
## TO DO: algorithm to remove these errors, but not 'real' outliers
## (untouched from AmsterdamUMCdb)
## Remove extreme outliers, most likely data entry errors (manual_entry = True)
creatinine.loc[
        (creatinine['value'] < 30) & (creatinine['manual_entry']),
    'value'] = np.nan
creatinine = creatinine.dropna(subset=['value'])

## Get highest creatinine per 24 hours
## Use creatinine 'cleansed' dataframe from APACHE score
sofa_renal_creatinine = creatinine.groupby(['admissionid','time']).agg(
        max_creatinine=pd.NamedAgg(column='value', aggfunc='max')
    ).reset_index()
sofa_renal_creatinine.head()
## Calculate SOFA renal score for creatinine:
sofa_renal_creatinine['sofa_renal_score'] = 0
## Creatinine 110-170 umol/l
sofa_renal_creatinine.loc[
        (sofa_renal_creatinine['max_creatinine'] >= 110) &
        (sofa_renal_creatinine['max_creatinine'] < 171),
     'sofa_renal_score'] = 1
## Creatinine 171-299 umol/l
sofa_renal_creatinine.loc[
        (sofa_renal_creatinine['max_creatinine'] >= 171) &
        (sofa_renal_creatinine['max_creatinine'] < 300),
    'sofa_renal_score'] = 2
## Creatinine 300-440 umol/l
sofa_renal_creatinine.loc[
        (sofa_renal_creatinine['max_creatinine'] >= 300) &
        (sofa_renal_creatinine['max_creatinine'] <= 440),
     'sofa_renal_score'] = 3
## Creatinine >440 umol/l
sofa_renal_creatinine.loc[
        (sofa_renal_creatinine['max_creatinine'] > 440),
    'sofa_renal_score'] = 4
sofa_renal_creatinine.head()

## Combine the scores from creatinine and urine output
sofa_renal = pd.concat(
        [sofa_renal_creatinine, sofa_daily_urine_out], sort=False,
    ).sort_values(by=['admissionid','time'])

sofa_renal.head()

########################
## Final SOFA scores

## Merge the scores
sofa = admissions_df['admissionid']
## Max respiration score
scores = sofa_respiration.groupby(['admissionid', 'time']).agg(
    sofa_respiration_score=pd.NamedAgg(
        column='sofa_respiration_score', aggfunc='max'))
scores = scores.sort_values(by=['admissionid', 'time']).reset_index()
sofa = pd.merge(sofa, scores, on='admissionid', how='left')

## Max coagulation score
scores = sofa_platelets.groupby(['admissionid', 'time']).agg(
    sofa_coagulation_score=pd.NamedAgg(
        column='sofa_coagulation_score', aggfunc='max'))

## We only want to keep columns where the 'day' (24hr periods post-admission)
## match, but some 'days' in the new SOFA component score may not be in the
## previous SOFA component scores (and we want to keep this)
## So for these rows, set all previous SOFA component scores to nan and then
## only keep matching 'days'
## Here time_x is the combined previous SOFA component scores 'day' and time_y
## is the new SOFA component score 'day'
def join_sofa(sofa, scores):
    scores = scores.sort_values(by=['admissionid', 'time']).reset_index()
    sofa = pd.concat([sofa, scores[['admissionid', 'time']]])
    sofa = sofa.loc[~sofa[['admissionid', 'time']].duplicated()]
    sofa = sofa.sort_values(by=['admissionid', 'time']).reset_index(drop=True)
    sofa = pd.merge(sofa, scores, on=['admissionid','time'], how='left')
    return sofa

sofa = join_sofa(sofa, scores)

## Max liver score
scores = sofa_bilirubin.groupby(['admissionid', 'time']).agg(
    sofa_liver_score=pd.NamedAgg(column='sofa_liver_score', aggfunc='max'))
sofa = join_sofa(sofa, scores)

## Max cardiovascular score
scores = sofa_cardiovascular.groupby(['admissionid', 'time']).agg(
    sofa_cardiovascular_score=pd.NamedAgg(
        column='sofa_cardiovascular_score', aggfunc='max'))
sofa = join_sofa(sofa, scores)

## Max central nervous system score
scores = sofa_cns.groupby(['admissionid', 'time']).agg(
    sofa_cns_score=pd.NamedAgg(column='sofa_cns_score', aggfunc='max'))
sofa = join_sofa(sofa, scores)

## Max renal score
scores = sofa_renal.groupby(['admissionid', 'time']).agg(
    sofa_renal_score=pd.NamedAgg(column='sofa_renal_score', aggfunc='max'))
sofa = join_sofa(sofa, scores)

## Calculate total score (add al values in columns)
total_scores = sofa.set_index(['admissionid','time']).sum(
    axis=1, skipna=True).to_frame('sofa_total_score')
sofa = pd.merge(sofa, total_scores, on=['admissionid','time'], how='left')
sofa.head()

sofa.dropna(subset=['time'], inplace=True)

## save as .csv file
sofa.to_csv(additional_file_path + 'sofa.csv', index=False)

## SOFA scores (as per AmsterdamUMCdb) completed

################################################################################
## Next for Sepsis-3 definition is antibiotics escalation

########################
## Antibiotics

antibiotics = pd.DataFrame(columns=['itemid', 'rank'])
antibiotics.loc[0] = [7185, 1] # Doxycycline (Vibramycine)
antibiotics.loc[1] = [9142, 2] # Tetracycline
antibiotics.loc[2] = [19764, 4] # Tigecycline (Tygacil)
## not included: Chlooramfenicol (Globenicol) (cream/ointment) 6932
antibiotics.loc[3] = [9047, 2] # Chlooramfenicol
antibiotics.loc[4] = [6847, 1] # Amoxicilline (Clamoxyl/Flemoxin)
antibiotics.loc[5] = [9128, 3] # Piperacilline (Pipcil)
antibiotics.loc[6] = [6871, 1] # Benzylpenicilline (Penicilline)
antibiotics.loc[7] = [9037, 1] # Feneticilline (Broxil)
## not included: Flucloxacilline (Stafoxil/Floxapen) 9070 (cream/ointment)
antibiotics.loc[8] = [9029, 2] # Amoxicilline/Clavulaanzuur (Augmentin)
antibiotics.loc[9] = [9152, 1] # Cefazoline (Kefzol)
antibiotics.loc[10] = [9151, 1] # Cefuroxim (Zinacef)
antibiotics.loc[11] = [6917, 2] # Ceftazidim (Fortum)
antibiotics.loc[12] = [9133, 2] # Ceftriaxon (Rocephin)
antibiotics.loc[13] = [9030, 4] # Aztreonam (Azactam)
antibiotics.loc[14] = [8127, 4] # Meropenem (Meronem)
antibiotics.loc[15] = [9070, 2] # Zilversulfadiazine (Flammazine)
antibiotics.loc[16] = [7208, 2] # Erythromycine (Erythrocine)
antibiotics.loc[17] = [8546, 2] # Claritromycine (Klacid)
antibiotics.loc[18] = [13057, 2] # Azitromycine (Zithromax)
antibiotics.loc[19] = [6958, 2] # Clindamycine (Dalacin)
antibiotics.loc[20] = [7044, 2] # Tobramycine (Obracin)
antibiotics.loc[21] = [13094, 2] # Tobramycine oogzalf (Tobrex)
antibiotics.loc[22] = [7235, 2] # Gentamicine (Garamycin)
## not included: Dexamethason/gentamicine oogzalf (Dexamytrex) 13102 (eye)
antibiotics.loc[23] = [9109, 2] # Neomycine sulfaat
antibiotics.loc[24] = [6834, 2] # Amikacine (Amukin)
## not included: Ofloxacine (Trafloxal) oogdruppels 12997 (eye drops)
antibiotics.loc[25] = [6948, 2] # Ciprofloxacine (Ciproxin)
antibiotics.loc[26] = [9117, 2] # Norfloxacine (Noroxin)
antibiotics.loc[27] = [12398, 2] # Levofloxacine (Tavanic)
antibiotics.loc[28] = [12956, 2] # Moxifloxacin (Avelox)
antibiotics.loc[29] = [7064, 2] # Vancomycine
antibiotics.loc[30] = [8549, 4] # Belcomycine (Colistinesulfaat) 4 x dgs
antibiotics.loc[31] = [10584, 4] # Belcomycine (Colistinesulfaat) 6 x dgs
antibiotics.loc[32] = [20175, 4] # Colistine
antibiotics.loc[33] = [20176, 4] # Colistine Inhalatie
## not included: Fusidinezuur (Fucidin) 9075 - prophylactic after
##  cardiothoracic surgery
## not included: Fusidinezuur oogdruppels (Fusithalmic) 13045 (eye drops)
antibiotics.loc[34] = [8942, 1] # Metronidazol-Flagyl
antibiotics.loc[35] = [7187, 1] # Metronidazol (Flagyl)
antibiotics.loc[36] = [14236, 2] # Nitrofurantoïne ( Furadantine)
antibiotics.loc[37] = [19137, 4] # Linezolid (Zyvoxid)
antibiotics.loc[38] = [19773, 4] # Daptomycine (Cubicin)
antibiotics.loc[39] = [8394, 1] # Co-Trimoxazol (Bactrimel)
antibiotics.loc[40] = [9052, 1] # Co-trimoxazol forte (Bactrimel)
antibiotics.loc[41] = [6919, 2] # Cefotaxim (Claforan)

## As part of the selective digestive decontamination, patients expected to stay
## at least 24-48hrs in ICU receive 16 (at least between 10 and 20) doses of
## cefotaxime across 4 days. If the clinician suspects an infection, this should
## be switched to  ceftriaxone (which has a similar spectrum). If cefotaxime is
## continued after these initial doses, assume the clinician has suspected an
## infection and kept cefotaxime.
## Similarly, for cardiothoracic surgery patients, antibiotic usage will be
## prophylactic (vancomycin/fusidinezuur annd cefazolin).
## Fusidinezuur is always prophylactic, vancomycin prophylactic in day 1.
## Prophylactic antibiotic usage is picked up again later in the script.

drugitems_abx = drugitems.loc[drugitems['itemid'].isin(antibiotics['itemid'])]

drugitems_abx.loc[
        (drugitems_abx['stop_time'] <= 4) &
        (drugitems_abx['itemid'] == 6919),
    'itemid'] = np.nan
drugitems_abx.loc[
        (drugitems_abx['stop_time'] <= 1) &
        (drugitems_abx['itemid'] == 7064),
    'itemid'] = np.nan
drugitems_abx.dropna(subset=['itemid'], inplace=True)

drugitems_abx = pd.merge(
    drugitems_abx, antibiotics[['itemid', 'rank']], on='itemid', how='left')
drugitems_abx['intravenous'] = drugitems_abx['ordercategoryid'].isin([15,65,55])
## Sepsis-3 (Shah et al.) say antibiotics considered on during 24hr period if
## administration occurred within 24hr period or within 12hrs before 24hr period
drugitems_abx['start_time'] = (
    drugitems_abx['start'] - drugitems_abx['admittedat'] - 1000*60*60*12)
drugitems_abx['start_time'] //= (1000*60*60*24)
## The same as before
drugitems_abx['stop_time'] = drugitems_abx['stop'] - drugitems_abx['admittedat']
drugitems_abx['stop_time'] //= (1000*60*60*24)
## As with the SOFA cardiovascular score, we want to add extra rows to the
## dataframe, for when drug administration happened over consecutive 'days'
## (i.e. make an entry for each 'day' tha the drug administration window
## overlaps with)
n_days = drugitems_abx['stop_time'] - drugitems_abx['start_time'] + 1
drugitems_abx = drugitems_abx.loc[
    drugitems_abx.index.repeat(n_days)].reset_index(drop=True)
drugitems_abx['time'] = np.hstack([np.arange(x) for x in n_days])
drugitems_abx['time'] += drugitems_abx['start_time']

abx_cols = ['admissionid', 'time', 'intravenous', 'rank', 'item', 'itemid']
drugitems_abx = drugitems_abx[abx_cols]
drugitems_abx = drugitems_abx.loc[~drugitems_abx.duplicated()]

max_rank = drugitems_abx.groupby(['admissionid', 'time']).agg(
        max_rank=pd.NamedAgg(column='rank', aggfunc='max')
    ).reset_index()
drugitems_abx = pd.merge(
    drugitems_abx, max_rank, on=['admissionid', 'time'], how='left')
drugitems_abx['rank_diff'] = drugitems_abx['max_rank'] - drugitems_abx['rank']

## Now compute the antibiotic escalation from day to day. Antibiotic escalation
## occurs if a new drug of higher ranking is administered (i.e. the max_rank
## increases) or if the number of drugs at the highest ranking is increased
## (i.e. n_max_rank increases), when this is accompanied by at least one
## antibiotics given intravenously
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
## If two consecutive rows correspond to different patients, then we don't want
## antibiotics escalation defined from the first patient to the second!
abx_escalation.loc[
        (~abx_escalation['admissionid'].duplicated()),
    ['max_rank_increase', 'n_max_rank_increase', 'time_diff']] = np.nan
abx_escalation['antibiotic_escalation'] = False
## Any new entry is assumed to be the first antibiotics given to that patient
abx_escalation.loc[
        (~abx_escalation['admissionid'].duplicated()),
    'antibiotic_escalation'] = True
abx_escalation.loc[
        (abx_escalation['max_rank_increase'] > 0) |
        (abx_escalation['n_max_rank_increase'] > 0) |
        (abx_escalation['time_diff'] != 1),
    'antibiotic_escalation'] = True

abx_escalation.drop(
    columns=['max_rank_increase', 'n_max_rank_increase', 'time_diff'],
    inplace=True)

################################################################################
## Now we are in a position to compute sepsis episodes according to the Sepsis-3
## definition of increase in SOFA score accompanied by antibiotics escalation

## Similar to SOFA scores, we want to include antibiotic escalation where no
## SOFA scores exist yet
sepsis = pd.concat([sofa, abx_escalation[['admissionid', 'time']]])
sepsis = sepsis.loc[~sepsis[['admissionid', 'time']].duplicated()]
sepsis = sepsis.sort_values(by=['admissionid', 'time']).reset_index(drop=True)
sepsis.loc[
        sepsis['sofa_total_score'].isna() & (sepsis['time'] < 0),
    'sofa_total_score'] = 0
sepsis = sepsis.dropna(subset=['time']).reset_index(drop=True)

## Shah et al. (DOI: 10.1097/CCM.0000000000005169) in Sepsis-3 review
## remove any patients with less than 3 SOFA scores in the first 24hrs,
## mark these patients in case we want to do the same
sofa_components = ['sofa_respiration_score', 'sofa_coagulation_score']
sofa_components += ['sofa_liver_score', 'sofa_cardiovascular_score']
sofa_components += ['sofa_cns_score', 'sofa_renal_score']
sepsis['n_sofa_scores'] = sepsis[sofa_components].notna().sum(axis=1)

sepsis['discard'] = False
sepsis.loc[
        (sepsis['n_sofa_scores'] < 3) &
        (sepsis['time'] == 0),
    'discard'] = True
sepsis.loc[
        sepsis['admissionid'].isin(sepsis.loc[sepsis['discard'],'admissionid']),
    'discard'] = True

## Merge sofa scores with antibiotic escalation ressults
sepsis = pd.merge(
    sepsis, abx_escalation, on=['admissionid', 'time'], how='left')

sepsis['time_diff'] = sepsis['time'].diff()
sepsis.loc[~sepsis['admissionid'].duplicated(), 'time_diff'] = np.nan

## SOFA increase corresponds to three time periods: previous and current day,
## current and subsequent day, previous and subsequent day (see Shah et al. for
## more details)
## Change in SOFA between previous and current 24hr periods
sepsis['sofa_diff0'] = sepsis['sofa_total_score'].diff()
sepsis[['sofa_diff1', 'sofa_diff2']] = np.nan
## Change in SOFA between current and subsequent 24hr periods
sepsis['sofa_diff1'].iloc[:-1] = sepsis['sofa_diff0'].iloc[1:]
## Change in SOFA between previous and subsequent 24hr periods
sepsis['sofa_diff2'].iloc[:-1] = (
    sepsis['sofa_total_score'].diff(periods=2).iloc[1:])

## If this is the first entry for the patient, then there is no 'previous day'
## So the difference in these cases is simply the (total score - 0)
## (we assume a SOFA score of 0 on admission)
sepsis.loc[~sepsis['admissionid'].duplicated(), ['sofa_diff0', 'sofa_diff2']] =(
    sepsis.loc[~sepsis['admissionid'].duplicated(), 'sofa_total_score'])
## On the last day for that patient, sofa_diff1 and sofa_diff2 should not exist
sepsis.loc[
        (~sepsis['admissionid'].duplicated(keep='last')),
    ['sofa_diff1', 'sofa_diff2']] = np.nan

## A sepsis episode is defined as antibiotics escalation ('infection')
## accompanied by SOFA increase of 2 or more.
sepsis['sepsis_episode'] = (
    (sepsis['antibiotic_escalation']) &
    ((sepsis['sofa_diff0'] >= 2) |
        (sepsis['sofa_diff1'] >= 2) |
        (sepsis['sofa_diff2'] >= 2)))

## Next: find patients who had elective surgery, these may have been given
## antibiotics for prophylactic use in the first 24hr after admission
## and are not classified as sepsis even if having a high SOFA score
## (subsequent antibiotic escalation and SOFA increase is marked as sepsis)
surgical_cols = ['admissionid', 'surgical', 'urgency']
sepsis = pd.merge(
    sepsis, combined_diagnoses[surgical_cols], on='admissionid', how='left')
sepsis['elective_surgery'] = (
    (sepsis['surgical'] == 1) &
    (sepsis['urgency'] == 0))
sepsis['emergency_surgery'] = (
    (sepsis['surgical'] == 1) &
    (sepsis['urgency'] == 1))

## Assume prophylactic treatment rather than sepsis if elective surgery
## accompanied by high SOFA/antibiotics (including a day after)
sepsis['prophylaxis'] = (sepsis['sepsis_episode'] & sepsis['elective_surgery'])
sepsis['prophylaxis'].iloc[:-1] |= (
    (sepsis['sepsis_episode'].iloc[:-1] &
    (sepsis.loc[1:,'elective_surgery'].values)))
## If prophylaxis is assumed, then not a sepsis episode
sepsis['sepsis_episode'] &= ~sepsis['prophylaxis']
sepsis['infection'] = (sepsis['antibiotic_escalation'] & ~sepsis['prophylaxis'])

## Lactate (for septic shock: max lactate of 2mmol/L + cardiovascular SOFA
## score of >=3 which corresponds to use of vasopressors)
lactate_itemids = [10053, 6837, 9580]
numerics_lactate = numerics_read(
    lactate_itemids, admissions_df=admissions_df,
    end=end_time, start=start_time)
lactate_max = numerics_lactate.groupby(['admissionid', 'time']).agg(
    max_lactate=pd.NamedAgg(column='value', aggfunc='max')).reset_index()

sepsis = pd.merge(
    sepsis, lactate_max[['admissionid', 'time', 'max_lactate']],
    on=['admissionid', 'time'], how='left')

sepsis['septic_shock'] = (
    (sepsis['sepsis_episode'] &
    (sepsis['max_lactate'] > 2) &
    (sepsis['sofa_cardiovascular_score'] >= 3)))

## This gives a binary variable for each patient if they had at least one
## sepsis episode.
sepsis_patients = sepsis.groupby(['admissionid']).agg(
        sepsis=pd.NamedAgg(column='sepsis_episode', aggfunc='any'),
        septic_shock=pd.NamedAgg(column='septic_shock', aggfunc='any')
    ).reset_index()

sepsis_cols = ['admissionid', 'time', 'sofa_total_score']
sepsis_cols += ['antibiotic_escalation', 'prophylaxis']
sepsis_cols += ['infection', 'sepsis_episode', 'septic_shock']
sepsis_output = sepsis[sepsis_cols]

## Save to csv
sepsis_output.to_csv(additional_file_path + 'sepsis3.csv', index=False)

################################################################################
## Here are the SQL queries from AmsterdamUMCdb SOFA script, from which this
## code is developed.

## FiO2/PaO2 SQL
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
                12279, --O2 concentratie --measurement by Servo-i/Servo-U ventilator
                12369, --SET %O2: used with BiPap Vision ventilator
                16246 --Zephyros FiO2: Non-invasive ventilation
            ) THEN TRUE
            ELSE FALSE
        END AS ventilatory_support,
        CASE
            WHEN n.itemid IN (
                --FiO2 settings on respiratory support
                6699, --FiO2 %: setting on Evita ventilator
                12279, --O2 concentratie --measurement by Servo-i/Servo-U ventilator
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
                            WHEN n.value >= 1 AND n.value < 2 THEN 0.22 -- not defined by NICE
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
            12279, --O2 concentratie --measurement by Servo-i/Servo-U ventilator
            12369, --SET %O2: used with BiPap Vision ventilator
            16246 --Zephyros FiO2: Non-invasive ventilation
        )
    --measurements within 24 hours of ICU stay:
    AND (n.measuredat - a.admittedat) <= 1000*60*60*24 AND (n.measuredat - a.admittedat) >= 0
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
        (fio2_table.measuredat - pao2.measuredat)/(60*1000) AS FiO2_time_difference,
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
        fio2_table.measuredat > pao2.measuredat - 60*60*1000 AND --no earlier than 60 minutes before pao2 measurement
        fio2_table.measuredat < pao2.measuredat + 15*60*1000 --no later than 15 minutes after pao2 measurement
    WHERE
        pao2.itemid IN (
            7433, --PO2
            9996, --PO2 (bloed)
            21214 --PO2 (bloed) - kPa
        )
    --measurements within 24 hours of ICU stay (use 30 minutes before admission to allow for time differences):
    AND (pao2.measuredat - a.admittedat) <= 1000*60*60*24 AND (pao2.measuredat - a.admittedat) >= -(1000*60*30) AND
    (f.value LIKE '%art.%' OR f.value IS NULL)  -- source is arterial or undefined (assume arterial)
)
SELECT * FROM oxygenation
WHERE priority = 1
"""


## Platelets SQL
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
--measurements within 24 hours of ICU stay (use 30 minutes before admission to allow for time differences):
AND (n.measuredat - a.admittedat) <= 1000*60*60*24 AND (n.measuredat - a.admittedat) >= -(1000*60*30)
"""
##
##
## Bilirubin SQL
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
--measurements within 24 hours of ICU stay (use 30 minutes before admission to allow for time differences):
AND (n.measuredat - a.admittedat) <= 1000*60*60*24 AND (n.measuredat - a.admittedat) >= -(1000*60*30)
"""


## Cardiovascular SQL
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
    WHEN doserateperkg == B'0' AND doseunitid = 11 AND doserateunitid = 4 --unit: µg/min -> µg/kg/min
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
    WHEN doserateperkg = B'0' AND doseunitid = 10 AND doserateunitid = 5 --unit: mg/uur  -> µg/kg/min
        THEN CASE
            WHEN patientweight > 0
            THEN dose*1000/patientweight/60
            ELSE dose*1000/80/60 --mean weight
        END
    WHEN doserateperkg = B'1' AND doseunitid = 11 AND doserateunitid = 4 --unit: µg/kg/min (no conversion needed)
        THEN dose
    WHEN doserateperkg = B'1' AND doseunitid = 11 AND doserateunitid = 5 --unit: µg/kg/uur -> µg/kg/min
        THEN dose/60
    END AS gamma
    FROM dosing
    WHERE
    -- medication given within 24 hours of ICU stay:
    start_time <= 24*60 AND stop_time >= 0
    ORDER BY admissionid, start_time
"""


## Mean ABP SQL
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
##
##
## Glasgow Coma Scale SQL
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
            WHEN 6734 THEN 7 - motor.valueid    --Beste motore reactie van de armen
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
    AND (eyes.measuredat - a.admittedat) <= 1000*60*60*24 AND (eyes.measuredat - a.admittedat) >= 0
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


## Urine output SQL
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
-- measurements within 24 hours of ICU stay (use 30 minutes before admission to allow for time differences):
AND (n.measuredat - a.admittedat) <= 1000*60*60*24 AND (n.measuredat - a.admittedat) >= 0
"""


## Creatinine SQL
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
    (n.measuredat - a.admittedat)/(60*60*1000) > -(365*24) AND (n.measuredat - a.admittedat) < (24*60*60*1000)
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
    (n.measuredat - a.admittedat) > 0 AND (n.measuredat - a.admittedat) < (7*24*60*60*1000)
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
        WHEN baseline_creatinine > 0 AND m.max_creatinine_7days/baseline_creatinine > 3 THEN TRUE
        -- AKI definition: increase to >= 354 umol/l AND at least 44 umol/l increase:
        WHEN max_creatinine_7days >= 354 AND max_creatinine_7days - baseline_creatinine >= 44 THEN TRUE
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
-- measurements within 24 hours of ICU stay (use 30 minutes before admission to allow for time differences):
AND (n.measuredat - a.admittedat) <= 1000*60*60*24 AND (n.measuredat - a.admittedat) >= -(1000*60*30)
"""
