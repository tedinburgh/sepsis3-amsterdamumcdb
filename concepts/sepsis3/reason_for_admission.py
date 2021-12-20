## Author: Tom Edinburgh
## v1: date 04/10/2021.

## This script computes the sepsis definition as per the AmsterdamUMCdb script,
## and write the results to .csv file. The main SQL query from this script
## (which has been adapted into pandas via the raw files) is at the bottom of
## the script. See:
## https://github.com//AmsterdamUMC//AmsterdamUMCdb/blob/master/concepts/diagnosis/reason_for_admission.ipynb

################################################################################

import pandas as pd
import numpy as np
import re

import amsterdamumcdb as adb
dictionary = adb.get_dictionary()

file_path = '../../data/'
additional_file_path = '../../data/additional_files/'

list_columns = ['admissionid', 'itemid', 'valueid', 'value']
list_columns += ['measuredat', 'updatedat', 'registeredby']
listitems = pd.read_csv(file_path + 'listitems.csv', usecols=list_columns)

drug_columns = ['admissionid', 'itemid', 'item', 'duration', 'rate', 'rateunit']
drug_columns += ['start', 'stop', 'dose', 'doserateperkg', 'doseunitid']
drug_columns += ['doserateunitid', 'ordercategoryid']
drugitems = pd.read_csv(file_path + 'drugitems.csv', usecols=drug_columns)

procedureorder_columns = ['admissionid', 'itemid', 'item', 'registeredat']
procedureorderitems = pd.read_csv(
    file_path + 'procedureorderitems.csv', usecols=procedureorder_columns)

freetextitems = pd.read_csv(file_path + 'freetextitems.csv')

admissions_df = pd.read_csv(file_path + 'admissions.csv')

################################################################################

## Get itemids for diagnosis related variables
## Main group (level 0)
level0_itemid = [13110] # D_Hoofdgroep
level0_itemid += [16651] # DMC_Hoofdgroep, Medium Care
level0_itemid += [18588] # Apache II Hoofdgroep
level0_itemid += [16997] # APACHE IV Groepen
level0_itemid += [18669] # NICE APACHEII diagnosen
level0_itemid += [18671] # NICE APACHEIV diagnosen

## Subgroup (level 1)
level1_itemid = [13111] # D_Subgroep_Thoraxchirurgie
level1_itemid += [16669] # DMC_Subgroep_Thoraxchirurgie
level1_itemid += [13112] # D_Subgroep_Algemene chirurgie
level1_itemid += [16665] # DMC_Subgroep_Algemene chirurgie
level1_itemid += [13113] # D_Subgroep_Neurochirurgie
level1_itemid += [16667] # DMC_Subgroep_Neurochirurgie
level1_itemid += [13114] # D_Subgroep_Neurologie
level1_itemid += [16668] # DMC_Subgroep_Neurologie
level1_itemid += [13115] # D_Subgroep_Interne geneeskunde
level1_itemid += [16666] # DMC_Subgroep_Interne geneeskunde

## Surgical
surgical_itemid = [13116] # D_Thoraxchirurgie_CABG en Klepchirurgie
surgical_itemid += [16671] # DMC_Thoraxchirurgie_CABG en Klepchirurgie
surgical_itemid += [13117] # D_Thoraxchirurgie_Cardio anders
surgical_itemid += [16672] # DMC_Thoraxchirurgie_Cardio anders
surgical_itemid += [13118] # D_Thoraxchirurgie_Aorta chirurgie
surgical_itemid += [16670] # DMC_Thoraxchirurgie_Aorta chirurgie
surgical_itemid += [13119] # D_Thoraxchirurgie_Pulmonale chirurgie
surgical_itemid += [16673] # DMC_Thoraxchirurgie_Pulmonale chirurgie
## Not surgical: 13141 # D_Algemene chirurgie_Algemeen
## Not surgical: 16642 # DMC_Algemene chirurgie_Algemeen
surgical_itemid += [13121] # D_Algemene chirurgie_Buikchirurgie
surgical_itemid += [16643] # DMC_Algemene chirurgie_Buikchirurgie
surgical_itemid += [13123] # D_Algemene chirurgie_Endocrinologische chirurgie
surgical_itemid += [16644] # DMC_Algemene chirurgie_Endocrinologische chirurgie
surgical_itemid += [13145] # D_Algemene chirurgie_KNO/Overige
surgical_itemid += [16645] # DMC_Algemene chirurgie_KNO/Overige
surgical_itemid += [13125] # D_Algemene chirurgie_Orthopedische chirurgie
surgical_itemid += [16646] # DMC_Algemene chirurgie_Orthopedische chirurgie
surgical_itemid += [13122] # D_Algemene chirurgie_Transplantatie chirurgie
surgical_itemid += [16647] # DMC_Algemene chirurgie_Transplantatie chirurgie
surgical_itemid += [13124] # D_Algemene chirurgie_Trauma
surgical_itemid += [16648] # DMC_Algemene chirurgie_Trauma
surgical_itemid += [13126] # D_Algemene chirurgie_Urogenitaal
surgical_itemid += [16649] # DMC_Algemene chirurgie_Urogenitaal
surgical_itemid += [13120] # D_Algemene chirurgie_Vaatchirurgie
surgical_itemid += [16650] # DMC_Algemene chirurgie_Vaatchirurgie
surgical_itemid += [13128] # D_Neurochirurgie _Vasculair chirurgisch
surgical_itemid += [16661] # DMC_Neurochirurgie _Vasculair chirurgisch
surgical_itemid += [13129] # D_Neurochirurgie _Tumor chirurgie
surgical_itemid += [16660] # DMC_Neurochirurgie _Tumor chirurgie
surgical_itemid += [13130] # D_Neurochirurgie_Overige
surgical_itemid += [16662] # DMC_Neurochirurgie_Overige
surgical_itemid += [18596] # Apache II Operatief  Gastr-intenstinaal
surgical_itemid += [18597] # Apache II Operatief Cardiovasculair
surgical_itemid += [18598] # Apache II Operatief Hematologisch
surgical_itemid += [18599] # Apache II Operatief Metabolisme
surgical_itemid += [18600] # Apache II Operatief Neurologisch
surgical_itemid += [18601] # Apache II Operatief Renaal
surgical_itemid += [18602] # Apache II Operatief Respiratoir
surgical_itemid += [17008] # APACHEIV Post-operative cardiovascular
surgical_itemid += [17009] # APACHEIV Post-operative gastro-intestinal
surgical_itemid += [17010] # APACHEIV Post-operative genitourinary
surgical_itemid += [17011] # APACHEIV Post-operative hematology
surgical_itemid += [17012] # APACHEIV Post-operative metabolic
surgical_itemid += [17013] # APACHEIV Post-operative musculoskeletal /skin
surgical_itemid += [17014] # APACHEIV Post-operative neurologic
surgical_itemid += [17015] # APACHEIV Post-operative respiratory
surgical_itemid += [17016] # APACHEIV Post-operative transplant
surgical_itemid += [17017] # APACHEIV Post-operative trauma

## Additional variables (level 2)
level2_itemid = surgical_itemid.copy()
level2_itemid += [13141] # D_Algemene chirurgie_Algemeen
level2_itemid += [16642] # DMC_Algemene chirurgie_Algemeen
level2_itemid += [13133] # D_Interne Geneeskunde_Cardiovasculair
level2_itemid += [16653] # DMC_Interne Geneeskunde_Cardiovasculair
level2_itemid += [13134] # D_Interne Geneeskunde_Pulmonaal
level2_itemid += [16658] # DMC_Interne Geneeskunde_Pulmonaal
level2_itemid += [13135] # D_Interne Geneeskunde_Abdominaal
level2_itemid += [16652] # DMC_Interne Geneeskunde_Abdominaal
level2_itemid += [13136] # D_Interne Geneeskunde_Infectieziekten
level2_itemid += [16655] # DMC_Interne Geneeskunde_Infectieziekten
level2_itemid += [13137] # D_Interne Geneeskunde_Metabool
level2_itemid += [16656] # DMC_Interne Geneeskunde_Metabool
level2_itemid += [13138] # D_Interne Geneeskunde_Renaal
level2_itemid += [16659] # DMC_Interne Geneeskunde_Renaal
level2_itemid += [13139] # D_Interne Geneeskunde_Hematologisch
level2_itemid += [16654] # DMC_Interne Geneeskunde_Hematologisch
level2_itemid += [13140] # D_Interne Geneeskunde_Overige
level2_itemid += [16657] # DMC_Interne Geneeskunde_Overige
level2_itemid += [13131] # D_Neurologie_Vasculair neurologisch
level2_itemid += [16664] # DMC_Neurologie_Vasculair neurologisch
level2_itemid += [13132] # D_Neurologie_Overige
level2_itemid += [16663] # DMC_Neurologie_Overige
level2_itemid += [13127] # D_KNO/Overige
level2_itemid += [18589] # Apache II Non-Operatief Cardiovasculair
level2_itemid += [18590] # Apache II Non-Operatief Gastro-intestinaal
level2_itemid += [18591] # Apache II Non-Operatief Hematologisch
level2_itemid += [18592] # Apache II Non-Operatief Metabolisme
level2_itemid += [18593] # Apache II Non-Operatief Neurologisch
level2_itemid += [18594] # Apache II Non-Operatief Renaal
level2_itemid += [18595] # Apache II Non-Operatief Respiratoir
level2_itemid += [16998] # APACHE IV Non-operative cardiovascular
level2_itemid += [16999] # APACHE IV Non-operative Gastro-intestinal
level2_itemid += [17000] # APACHE IV Non-operative genitourinary
level2_itemid += [17001] # APACHEIV  Non-operative haematological
level2_itemid += [17002] # APACHEIV  Non-operative metabolic
level2_itemid += [17003] # APACHEIV Non-operative musculo-skeletal
level2_itemid += [17004] # APACHEIV Non-operative neurologic
level2_itemid += [17005] # APACHEIV Non-operative respiratory
level2_itemid += [17006] # APACHEIV Non-operative transplant
level2_itemid += [17007] # APACHEIV Non-operative trauma
## Both NICE APACHEII/IV also count towards surgical if valueid in correct range
level2_itemid += [18669] # NICE APACHEII diagnosen
level2_itemid += [18671] # NICE APACHEIV diagnosen

string_aggfunc = lambda x: '; '.join(v for v in x.unique())

diagnosis_groups = listitems.loc[listitems['itemid'].isin(level0_itemid)]
nice_ind = diagnosis_groups['itemid'].isin([18669, 18671])
diagnosis_groups.loc[nice_ind, 'value'] = \
    diagnosis_groups.loc[nice_ind, 'value'].apply(lambda x: x.split(' - ')[0])
diagnosis_groups.rename(columns={
        'value': 'diagnosis_group',
        'valueid': 'diagnosis_group_id'
    }, inplace=True)
diagnosis_groups['typeid'] = np.nan
diagnosis_groups.loc[diagnosis_groups['itemid'] == 18671, 'typeid'] = 6
diagnosis_groups.loc[diagnosis_groups['itemid'] == 18669, 'typeid'] = 5
diagnosis_groups.loc[
        (diagnosis_groups['itemid'] >= 16998) &
        (diagnosis_groups['itemid'] <= 17017),
    'typeid'] = 4
diagnosis_groups.loc[
        (diagnosis_groups['itemid'] >= 18589) &
        (diagnosis_groups['itemid'] <= 18602),
    'typeid'] = 3
diagnosis_groups.loc[
        (diagnosis_groups['itemid'] >= 13116) &
        (diagnosis_groups['itemid'] <= 13145),
    'typeid'] = 2
diagnosis_groups.loc[
        (diagnosis_groups['itemid'] >= 16642) &
        (diagnosis_groups['itemid'] <= 16673),
    'typeid'] = 1
## Add row number (ordered by typeid/measuredat)
diagnosis_groups = diagnosis_groups[~diagnosis_groups.duplicated()]
## See ticket #49 in AmsterdamUMCdb GitHub for discussion of sorting
diagnosis_groups = diagnosis_groups.sort_values(
    by=['admissionid', 'typeid', 'updatedat'], ascending=False)
n_entries = diagnosis_groups.groupby(
    ['admissionid']).size().sort_index(ascending=False)
diagnosis_groups['rownum'] = np.hstack([np.arange(x) for x in n_entries])
diagnosis_groups.sort_index(inplace=True)


diagnosis_subgroups = listitems.loc[listitems['itemid'].isin(level1_itemid)]
diagnosis_subgroups.rename(columns={
        'value': 'diagnosis_subgroup',
        'valueid': 'diagnosis_subgroup_id'
    }, inplace=True)
diagnosis_subgroups = diagnosis_subgroups.sort_values(
    by=['admissionid', 'updatedat'], ascending=False)
n_entries = diagnosis_subgroups.groupby(
    ['admissionid']).size().sort_index(ascending=False)
diagnosis_subgroups['rownum'] = np.hstack([np.arange(x) for x in n_entries])
diagnosis_subgroups.sort_index(inplace=True)


diagnoses = listitems.loc[listitems['itemid'].isin(level2_itemid)]
## There is a bug without the next line, as valueid 29 is
## 'Non-operatief Cardiovasculair -Coronair lijden', which is missing the
## extra space after the main diagnosis type
nice_ind = diagnoses['itemid'].isin([18669, 18671])
diagnoses.loc[nice_ind, 'value'] = diagnoses.loc[nice_ind, 'value'].apply(
    lambda x: x.replace(' -', ' - ').split(' - ')[1])
diagnoses.rename(columns={
        'value': 'diagnosis',
        'valueid': 'diagnosis_id'
    }, inplace=True)
diagnoses['surgical'] = diagnoses['itemid'].isin(surgical_itemid) * 1
diagnoses.loc[
        (diagnoses['itemid'] == 18669) &
        (diagnoses['diagnosis_id'] >= 1) &
        (diagnoses['diagnosis_id'] <= 26),
    'surgical'] = 1
diagnoses.loc[
        (diagnoses['itemid'] == 18671) &
        (diagnoses['diagnosis_id'] >= 222) &
        (diagnoses['diagnosis_id'] <= 452),
    'surgical'] = 1

diagnoses['typeid'] = np.nan
diagnoses.loc[diagnoses['itemid'] == 18671, 'typeid'] = 6
diagnoses.loc[diagnoses['itemid'] == 18669, 'typeid'] = 5
diagnoses.loc[
        (diagnoses['itemid'] >= 16998) &
        (diagnoses['itemid'] <= 17017),
    'typeid'] = 4
diagnoses.loc[
        (diagnoses['itemid'] >= 18589) &
        (diagnoses['itemid'] <= 18602),
    'typeid'] = 3
diagnoses.loc[
        (diagnoses['itemid'] >= 13116) &
        (diagnoses['itemid'] <= 13145),
    'typeid'] = 2
diagnoses.loc[
        (diagnoses['itemid'] >= 16642) &
        (diagnoses['itemid'] <= 16673),
    'typeid'] = 1
diagnoses['diagnosis_id'] = diagnoses['diagnosis_id'].astype(str)
## NOTE: this is different from AmsterdamUMCdb script. In their script,
## 'surgical' is not incldued in list, and if there are more than one entries
## that tie for the same admissionid/typeid/measuredat, then it's not clear
## which comes first (randomly assigned). This can lead to discrepancies between
## implementations of this code.
diagnoses = diagnoses.groupby(['admissionid', 'typeid', 'updatedat']).agg(
        surgical=pd.NamedAgg(column='surgical', aggfunc='any'),
        diagnosis=pd.NamedAgg(column='diagnosis', aggfunc=string_aggfunc),
        diagnosis_id=pd.NamedAgg(column='diagnosis_id', aggfunc=string_aggfunc)
    ).reset_index()
diagnoses['diagnosis_type'] = diagnoses['typeid'].replace({
    6: 'NICE APACHE IV',
    5: 'NICE APACHE II',
    4: 'APACHE IV',
    3: 'APACHE II',
    2: 'Legacy ICU',
    1: 'Legacy MCU'})
diagnoses = diagnoses.sort_values(
    by=['admissionid', 'typeid', 'updatedat'], ascending=False)
n_entries = diagnoses.groupby(
    ['admissionid']).size().sort_index(ascending=False)
diagnoses['rownum'] = np.hstack([np.arange(x) for x in n_entries])
diagnoses.sort_index(inplace=True)

################################################################################
## Sepsis definition
sepsis = listitems.loc[listitems['itemid'] == 15808]
sepsis['sepsis_at_admission'] = sepsis['valueid'].replace({2: 0})
sepsis['typeid'] = np.nan
## As above! Discrepancies occur if sort_values used instead of groupby
sepsis = sepsis.groupby(['admissionid', 'updatedat']).agg(
        sepsis_at_admission=pd.NamedAgg(
            column='sepsis_at_admission', aggfunc='any')
    ).reset_index()
sepsis = sepsis.sort_values(by=['admissionid', 'updatedat'], ascending=False)
n_entries = sepsis.groupby(['admissionid']).size().sort_index(ascending=False)
sepsis['rownum'] = np.hstack([np.arange(x) for x in n_entries])
sepsis.sort_index(inplace=True)

########################
## Sepsis antibiotics
sepsis_abx_itemid = [6834] # Amikacine (Amukin)
sepsis_abx_itemid += [6847] # Amoxicilline (Clamoxyl/Flemoxin)
sepsis_abx_itemid += [6871] # Benzylpenicilline (Penicilline)
sepsis_abx_itemid += [6917] # Ceftazidim (Fortum)
## Not included: 6919 # Cefotaxim (Claforan): prophylaxis
sepsis_abx_itemid += [6948] # Ciprofloxacine (Ciproxin)
sepsis_abx_itemid += [6953] # Rifampicine (Rifadin)
sepsis_abx_itemid += [6958] # Clindamycine (Dalacin)
sepsis_abx_itemid += [7044] # Tobramycine (Obracin)
## Not included: 7064 # Vancomycine: prophylaxis for valve surgery
sepsis_abx_itemid += [7123] # Imipenem (Tienam)
sepsis_abx_itemid += [7185] # Doxycycline (Vibramycine)
## Not included: 7187 # Metronidazol (Flagyl): often used for GI surgical
##  prophylaxis
## Not included: 7208 # Erythromycine (Erythrocine): often used for
##  gastroparesis
sepsis_abx_itemid += [7227] # Flucloxacilline (Stafoxil/Floxapen)
sepsis_abx_itemid += [7231] # Fluconazol (Diflucan)
sepsis_abx_itemid += [7232] # Ganciclovir (Cymevene)
sepsis_abx_itemid += [7233] # Flucytosine (Ancotil)
sepsis_abx_itemid += [7235] # Gentamicine (Garamycin)
sepsis_abx_itemid += [7243] # Foscarnet trinatrium (Foscavir)
sepsis_abx_itemid += [7450] # Amfotericine B (Fungizone)
## Not included: 7504 # X nader te bepalen --non-stock medication
sepsis_abx_itemid += [8127] # Meropenem (Meronem)
sepsis_abx_itemid += [8229] # Myambutol (ethambutol)
sepsis_abx_itemid += [8374] # Kinine dihydrocloride
## Not included: 8375 # Immunoglobuline (Nanogam): not anbiotic
## Not included: 8394 # Co-Trimoxazol (Bactrimel): often
##  prophylactic (unless high dose)
sepsis_abx_itemid += [8547] # Voriconazol(VFEND)
## Not included: 9029 # Amoxicilline/Clavulaanzuur (Augmentin): often used for
##  ENT surgical prophylaxis
sepsis_abx_itemid += [9030] # Aztreonam (Azactam)
sepsis_abx_itemid += [9047] # Chlooramfenicol
## Not included: 9075 # Fusidinezuur (Fucidin): prophylaxis
sepsis_abx_itemid += [9128] # Piperacilline (Pipcil)
sepsis_abx_itemid += [9133] # Ceftriaxon (Rocephin)
## Not included: 9151 # Cefuroxim (Zinacef): often used for GI/transplant
##  surgical prophylaxis
## Not included: 9152 # Cefazoline (Kefzol): prophylaxis for cardiac surgery
sepsis_abx_itemid += [9458] # Caspofungine
sepsis_abx_itemid += [9542] # Itraconazol (Trisporal)
## Not included: 9602 # Tetanusimmunoglobuline: prophylaxis/not antibiotic
sepsis_abx_itemid += [12398] # Levofloxacine (Tavanic)
sepsis_abx_itemid += [12772] # Amfotericine B lipidencomplex  (Abelcet)
sepsis_abx_itemid += [15739] # Ecalta (Anidulafungine)
sepsis_abx_itemid += [16367] # Research Anidulafungin/placebo
sepsis_abx_itemid += [16368] # Research Caspofungin/placebo
sepsis_abx_itemid += [18675] # Amfotericine B in liposomen (Ambisome )
sepsis_abx_itemid += [19137] # Linezolid (Zyvoxid)
sepsis_abx_itemid += [19764] # Tigecycline (Tygacil)
sepsis_abx_itemid += [19773] # Daptomycine (Cubicin)
sepsis_abx_itemid += [20175] # Colistine

sepsis_abx = drugitems.loc[
        (drugitems['itemid'].isin(sepsis_abx_itemid)) &
        (drugitems['start'] < 24*60*60*1000)]
sepsis_abx = sepsis_abx[['admissionid', 'item']]
sepsis_abx['sepsis_abx_bool'] = 1
sepsis_abx = sepsis_abx.groupby(['admissionid', 'sepsis_abx_bool']).agg(
        sepsis_abx_given=pd.NamedAgg(column='item', aggfunc=string_aggfunc)
    ).reset_index()

## Other antibiotics
other_abx_itemid = [7064] # Vancomycine: prophylaxis for valve surgery
## Not included: 7187 # Metronidazol (Flagyl): often used for GI surgical
##  prophylaxis
other_abx_itemid += [7187]
## Not included: 8394 # Co-Trimoxazol (Bactrimel): often prophylactic
##  (unless high dose)
other_abx_itemid += [8394]
## Not included: 9029 # Amoxicilline/Clavulaanzuur (Augmentin): often used for
##  ENT surgical prophylaxis
other_abx_itemid += [9029]
## Not included: 9151 # Cefuroxim (Zinacef): often used for GI surgical
##  prophylaxis
other_abx_itemid += [9151]
other_abx_itemid += [9152] # Cefazoline (Kefzol): prophylaxis

other_abx = drugitems.loc[
        (drugitems['itemid'].isin(other_abx_itemid)) &
        (drugitems['start'] < 24*60*60*1000)]
other_abx = other_abx[['admissionid', 'item']]
other_abx['other_abx_bool'] = 1
other_abx = other_abx.groupby(['admissionid', 'other_abx_bool']).agg(
        other_abx_given=pd.NamedAgg(column='item', aggfunc=string_aggfunc)
    ).reset_index()

########################
## Cultures drawn
cultures_itemid = [9189] # Bloedkweken afnemen
## Not included: 8097 # Sputumkweek afnemen: often used routinely
## Not included: 8418 # Urinekweek afnemen
## Not included: 8588 # MRSA kweken afnemen
cultures_itemid += [9190] # Cathetertipkweek afnemen
## Not included: 9191 # Drainvochtkweek afnemen
## Not included: 9192 # Faeceskweek afnemen: Clostridium
## Not included: 9193 # X-Kweek nader te bepalen
## Not included: 9194 # Liquorkweek afnemen
## Not included: 9195 # Neuskweek afnemen
## Not included: 9197 # Perineumkweek afnemen: often used routinely
## Not included: 9198 # Rectumkweek afnemen: often used routinely
cultures_itemid += [9200] # Wondkweek afnemen
cultures_itemid += [9202] # Ascitesvochtkweek afnemen
## Not included: 9203 # Keelkweek afnemen: often used routinely
## Not included: 9204 # SDD-kweken afnemen: often used routinely
cultures_itemid += [9205] # Legionella sneltest (urine)
## Not included: 1302 # SDD Inventarisatiekweken afnemen: often used routinely
## Not included: 19663 # Research Neuskweek COUrSe
## Not included: 19664 # Research Sputumkweek COUrSe

cultures = procedureorderitems.loc[
    procedureorderitems['itemid'].isin(cultures_itemid) &
    (procedureorderitems['registeredat'] < 6*60*60*1000)]
cultures = cultures[['admissionid', 'item']]
cultures['sepsis_cultures_bool'] = 1
cultures = cultures.groupby(['admissionid', 'sepsis_cultures_bool']).agg(
        sepsis_cultures_drawn=pd.NamedAgg(column='item', aggfunc=string_aggfunc)
    ).reset_index()

########################
## Merge for a single admissions dataframe
diagnosis_groups_cols = ['admissionid', 'diagnosis_group', 'diagnosis_group_id']
## This happens at the end of the SQL query
diagnosis_groups_add = diagnosis_groups.loc[(diagnosis_groups['rownum'] == 0)]
diagnosis_groups_add = diagnosis_groups_add[diagnosis_groups_cols]

diagnosis_subgroups_cols = ['admissionid', 'diagnosis_subgroup']
diagnosis_subgroups_cols += ['diagnosis_subgroup_id']
## This happens at the end of the SQL query
diagnosis_subgroups_add = diagnosis_subgroups.loc[
    (diagnosis_subgroups['rownum'] == 0)]
diagnosis_subgroups_add = diagnosis_subgroups_add[diagnosis_subgroups_cols]

diagnoses_cols = ['admissionid', 'diagnosis_type', 'diagnosis']
diagnoses_cols += ['diagnosis_id', 'surgical']
## This happens at the end of the SQL query
diagnoses_add = diagnoses.loc[(diagnoses['rownum'] == 0), diagnoses_cols]

sepsis_cols = ['admissionid', 'sepsis_at_admission']
sepsis_add = sepsis.loc[(sepsis['rownum'] == 0), sepsis_cols]

combined_diagnoses = pd.merge(
    admissions_df, diagnoses_add, on='admissionid', how='left')
combined_diagnoses = pd.merge(
    combined_diagnoses, diagnosis_groups_add, on='admissionid', how='left')
combined_diagnoses = pd.merge(
    combined_diagnoses, diagnosis_subgroups_add, on='admissionid', how='left')
combined_diagnoses = pd.merge(
    combined_diagnoses, sepsis_add, on='admissionid', how='left')
combined_diagnoses = pd.merge(
    combined_diagnoses, sepsis_abx, on='admissionid', how='left')
combined_diagnoses = pd.merge(
    combined_diagnoses, other_abx, on='admissionid', how='left')
combined_diagnoses = pd.merge(
    combined_diagnoses, cultures, on='admissionid', how='left')

## Up to this point covers the SQL query below

################################################################################
## Fixes for 'non-emergency medical' admissions (i.e. urgency=0, surgical=0)

surgical_dg = ['Algemene chirurgie', 'Thoraxchirurgie', 'Neurochirurgie']
surgical_dsg = ['Buikchirurgie', 'Vaatchirurgie', 'Orthopedische chirurgie']
surgical_dsg += ['CABG en Klepchirurgie', 'Pulmonale chirurgie']
surgical_dsg += ['Tumor chirurgie', 'Vasculair chirurgisch', 'Aorta chirurgie']
surgical_dsg += ['Transplantatie chirurgie', 'Endocrinologische chirurgie']
combined_diagnoses.loc[(
        (combined_diagnoses['urgency'] == 0) &
        (combined_diagnoses['surgical'] == 0) &
        (combined_diagnoses['diagnosis_group'].isin(surgical_dg) |
            combined_diagnoses['diagnosis_subgroup'].isin(surgical_dsg))),
    'surgical'] = True

combined_diagnoses.loc[(
        (combined_diagnoses['urgency'] == 0) &
        (combined_diagnoses['surgical'] == 0)),
    'urgency'] = 1

################################################################################

re_sepsis_surg = r'sepsis|pneumoni|GI perforation|perforation/rupture|'
re_sepsis_surg += r'infection|abscess|GI Vascular ischemia|diverticular|'
re_sepsis_surg += 'appendectomy|peritonitis'
re_sepsis_med = r'sepsis|septic|infect|pneumoni|cholangitis|pancr|endocarditis|'
re_sepsis_med += r'meningitis|GI perforation|abces|abscess|darm ischaemie|GI '
re_sepsis_med += 'vascular|fasciitis|inflammatory|peritonitis'

combined_diagnoses['sepsis_surgical'] = (combined_diagnoses['surgical'] == 1)
combined_diagnoses['sepsis_surgical'] &= (
    combined_diagnoses['diagnosis'].str.contains(
        re_sepsis_surg, na=False, flags=re.IGNORECASE))
combined_diagnoses['sepsis_surgical'] &= (
    ~(combined_diagnoses['sepsis_at_admission'] == 0))

combined_diagnoses['sepsis_med'] = (combined_diagnoses['surgical'] == 0)
combined_diagnoses['sepsis_med'] &= (
    combined_diagnoses['diagnosis'].str.contains(
        re_sepsis_med, na=False, flags=re.IGNORECASE))

## Medical admissions with sepsis
combined_diagnoses['sepsis'] = combined_diagnoses['sepsis_med']
combined_diagnoses['sepsis'] |= (
    (combined_diagnoses['sepsis_at_admission'] == 1) |
    (combined_diagnoses['sepsis_abx_bool'] == 1))
combined_diagnoses['sepsis'] |= (
    (combined_diagnoses['other_abx_bool'] == 1) &
    (combined_diagnoses['sepsis_cultures_bool'] == 1))
combined_diagnoses['sepsis'] &= (
    ~(combined_diagnoses['sepsis_at_admission'] == 0))

## Add surgical admissions with sepsis
combined_diagnoses['sepsis'] |= combined_diagnoses['sepsis_surgical']

################################################################################
## Also compute all the reasons for admission for comparison to AmsterdamUMCdb
## script.

re_cardiosurg = r'(CABG|AVR|hartchirurgie|heart surgery|'
re_cardiosurg += r'Chron. cardiovasculaire ziekte|hartkleppen|cardiovascula|'
re_cardiosurg += r'MVP|MVR|mitral|tricuspid|pericard|aortic.*valve|lobectom|'
re_cardiosurg += r'segment|thorax|Bentall|aorta-ascendens|aorta-boog|'
re_cardiosurg += r'aorta-wortel|aorta-descendens|lung|pneumectomie|bullectom|'
re_cardiosurg += r'respiratoir neoplasm|thoracoscop|thoracotom(y|ie)|respirato|'
re_cardiosurg += r'vrije wand ruptuur|VSR|ASD|pleurectom|intracardiac|'
re_cardiosurg += r'aneurysmectom|congenital defect repair)(?! for esophag)'

re_neurosurg = r'neuro|tentorieel|cranial|subdur|epidur|subarachn|intracerbr|'
re_neurosurg += r'hoofdtrauma|SAB|S.A.H.|CNS|Hoofd|seizures|craniotom|'
re_neurosurg += r'cranioplast|spinal|dwarslaesie|ventriculstom|transphenoidal|'
re_neurosurg += r'brain|A.V.M.|Arteriovenous malformation|brughoek|'
re_neurosurg += r'shunts and revisions|stereotactic|Burr hole|cerebrospinal'

re_vascsurg = r'vaatchirurgie|vasc.*surg|EVAR|aorta vervanging|perifeer vasc|'
re_vascsurg += r'embolectom|aneurysm|carotid|endovasc|dissectie|endarterectom|'
re_vascsurg += r'thrombectomy|dilatation|PTCA|all other bypass|transplantectom|'
re_vascsurg += r'femoral-popliteal|aorto-femoral|femoral-femoral'

re_gisurg = r'oesophagus|esophageal|maag|stomach|darm|intestin|gastro-intestin|'
re_gisurg += r'pancreatitis|laparotom|gastro-intestinale perforatie|galblaas|'
re_gisurg += r'Bleeding-.*GI|other GI|colon|rectal|GI.*surgery|GI obstruction|'
re_gisurg += r'Whipple|diverticular|appendectomy|peritonitis|cholecystectomy|'
re_gisurg += r'exenteration'

re_uro = r'(?<!ne)(urolog|cystectomy|genitourinary surgery|prostatectom|'
re_uro += r'ileal\-conduit|orchiectomy|bladder repair|nefrectom|nephrectom'
re_uro += r'|renaal neopsplama)'

re_gensurg = r'mond/keel|cancer oral|cancer laryngeal|spondylodes|'
re_gensurg += r'Fusion-spinal|devices for spine|orthopedic|renaal|metabol|'
re_gensurg += r'endocrin|thyroid|hip replacement|knee replacement|adrenalectom|'
re_gensurg += r'tracheostomy|bewaking|amputation|apnea-sleep|mastectomy|'
re_gensurg += r'|lymph node dissection|cosmetic|fracture-pathological|bewaking'

re_trauma_surg = r'(?<!non-)(?<!see )(trauma|hypotherm|'
re_trauma_surg += r'smoke inhalation)(?!, see trauma)(?! see)(?!: see)'

re_tx_surg = r'niertransplantatie|kidney transplant|renaal transplantatie|'
re_tx_surg += r'pancreastransplantatie'

re_respfailure_surg = r'resp.*insuff|na respiratoir arrest|arrest, respiratory|'
re_respfailure_surg += r'atelectas|pneumoni|na ok'

re_obgyn = r'hysterectom|Cesarean|ectopic pregnancy'
re_cardiacarrest_surg = r'Cardiac arrest.*respiratory arrest|na reanimatie'
re_hepa = r'lever'
re_surg_other = r'diagnose anders|respiratoir|cardiovasculair|niet operatief'

re_surgical_medical = r'Bloeding tractus digestivus|Haemorragische shock|'
re_surgical_medical += r'Gastro-intestinale bloeding|Bleeding, upper GI|Renaal|'
re_surgical_medical += r'hematologisch|Hematologische maligniteit|'
re_surgical_medical += r'Haematologisch'

re_surgical = '(' + re_cardiosurg + '|' + re_neurosurg + '|' + re_vascsurg + '|'
re_surgical += re_gisurg + '|' + re_uro + '|' + re_obgyn + '|' + re_gensurg
re_surgical += '|' + re_trauma_surg + '|' + re_tx_surg + '|' + re_hepa + '|'
re_surgical += re_surg_other + '|' + re_surgical_medical + '|'
re_surgical += re_respfailure_surg + '|' + re_sepsis_surg + '|'
re_surgical += re_cardiacarrest_surg + ')' # Deleted extra '|'

## Medications
re_respfailure_med = r'(?<! without )(resp.*insuff|pneumoni|respirato|'
re_respfailure_med += r'luchtweg obstructie|obstruction-airway|'
re_respfailure_med += r'chronisch obstructieve longziekte|emphysema|asthma|'
re_respfailure_med += r'aspiratie|aspiration|longembolie|pulmonary|pulmonaire'
re_respfailure_med += r'long|lung|atelectas|ALI|ARDS|pleural|cancer, lung|'
re_respfailure_med += r'pneumothorax|verdrinking|near drowning|weaning|'
re_respfailure_med += r'hemothorax|PCP)'

re_cardiacarrest_med = r'Cardiac arrest.*respiratory arrest|na reanimatie'

re_cardio = r'cardiogene shock|shock, cardiogenic|angina|ritme|rhythm|'
re_cardio += r'cardiovascular|cardiovasculair|myocardial|endocarditis|'
re_cardio += r'coronair|coronary|cardiomyopath|tamponade|pericardial|CHF|'
re_cardio += r'papillary muscle|^MI|hartkleppen|hart falen|decompensatio cordis'

re_neuro = r'(?<!see )(insult|seizure|CVA|observatie neurologische status|'
re_neuro += r'intracerebraal haematoom|intracranial|intracerebr|subdur|'
re_neuro += r'subarachno|epidur|coma|neurologisch|neurologic|CZS|S\.A\.B\.|'
re_neuro += r'neurologie|hoofdtrauma|head|neuro|muscula|spinal|meningitis|'
re_neuro += r'encephalitis|myasthenia|vaatspasme protocol|Guillian-Barre|'
re_neuro += r'encephalopath|musculoskeletal)(?!-see Neurological System)'

re_bleeding = r'bloeding tractus digestivus|gastro-intestinale bloeding|'
re_bleeding += r'gastro-intestinaal|bleeding, upper GI|bleeding, lower GI|'
re_bleeding += r'bleeding, GI|ulcer'

re_trauma_med = r'(?<!non-)(?<!see )(trauma|hypotherm|'
re_trauma_med += r'smoke inhalation)(?!, see trauma)(?! see)(?!: see)'

re_hemonc = r'malign|hematolog|cancer|bloeding|TTP|HUS|leukemi|pancytopen|'
re_hemonc += r'coagulopath|anemia|neutropen|lymph|sickle'

re_endo_med = r'metabolisme|keto-acidose|diabetic|metabolic|endocrine|'
re_endo_med += r'hypertens|acid-base|hypoglyc|thyroid'

re_vasc_med = r'aneurysm|vascular|cardiovascular medical|Thrombus,arterial|'
re_vasc_med += r'ascular medical|dissect|EVAR|embolectom'

re_gi_med = r'ileus|GI medical|GI obstruction'
re_tox = r'intox|overdosis|overdose|toxicity|withdrawal|drug'
re_shock_med = r'hypovolemische shock|shock|haemorr|hemorrha|anaphylaxis'
re_nefro_med = r'renaal|renal|tubulus|genitourinary|urolog|nefr'
re_hepa_med = r'lever|hepatic'
re_obgyn_med = r'obstetrie|postpartum|eclampsia'
re_mon_med = r'apnea, sleep|monitoring|bewaking|observatie'
re_tx_med = r'transplant'

re_medical = '(' + re_respfailure_med + '|' + re_cardiacarrest_med + '|'
re_medical += re_sepsis_med + '|' + re_cardio + '|' + re_neuro + '|'
re_medical += re_bleeding + '|' + re_gi_med + '|' + re_tox + '|' + re_trauma_med
re_medical += '|' + re_hemonc + '|' + re_endo_med + '|' + re_shock_med + '|'
re_medical += re_nefro_med + '|' + re_hepa_med + '|' + re_obgyn_med + '|'
## Deleted extra '|'
re_medical += re_vasc_med + '|' + re_mon_med + '|' + re_tx_med + ')'

## Compute reasons for admission
combined_diagnoses['cardiosurg'] = (
    (combined_diagnoses['surgical'] == 1) &
    (combined_diagnoses['diagnosis'].str.contains(
        re_cardiosurg, na=False, flags=re.IGNORECASE)))

combined_diagnoses['respfailure'] = (
    (combined_diagnoses['surgical'] == 1) &
    (combined_diagnoses['diagnosis'].str.contains(
        re_respfailure_surg, na=False, flags=re.IGNORECASE)))

combined_diagnoses['respfailure'] |= (
    (combined_diagnoses['surgical'] == 0) &
    (combined_diagnoses['diagnosis'].str.contains(
        re_respfailure_med, na=False, flags=re.IGNORECASE)))

combined_diagnoses['neurosurg'] = (
    (combined_diagnoses['surgical'] == 1) &
    (combined_diagnoses['diagnosis'].str.contains(
        re_neurosurg, na=False, flags=re.IGNORECASE)))

combined_diagnoses['gisurg'] =  (
    (combined_diagnoses['surgical'] == 1) &
    (combined_diagnoses['diagnosis'].str.contains(
        re_gisurg, na=False, flags=re.IGNORECASE)))

combined_diagnoses['cardiacarrest'] = (
    (combined_diagnoses['surgical'] == 1) &
    (combined_diagnoses['diagnosis'].str.contains(
        re_cardiacarrest_surg, na=False, flags=re.IGNORECASE)))
combined_diagnoses['cardiacarrest'] |= (
    (combined_diagnoses['surgical'] == 0) &
    (combined_diagnoses['diagnosis'].str.contains(
        re_cardiacarrest_med, na=False, flags=re.IGNORECASE)))

combined_diagnoses['vascsurg'] = (
    (combined_diagnoses['surgical'] == 1) &
    (combined_diagnoses['diagnosis'].str.contains(
        re_vascsurg, na=False, flags=re.IGNORECASE)))
## exclude cerebral aneurysms
combined_diagnoses['vascsurg'] &= (
    combined_diagnoses['diagnosis_group'].str.contains(
        'neuro', na=False, flags=re.IGNORECASE) == 0)

combined_diagnoses['trauma'] = (
    (combined_diagnoses['surgical'] == 1) &
    (combined_diagnoses['diagnosis'].str.contains(
        re_trauma_surg, na=False, flags=re.IGNORECASE)))
combined_diagnoses['trauma'] |= (
    (combined_diagnoses['surgical'] == 0) &
    (combined_diagnoses['diagnosis'].str.contains(
        re_trauma_med, na=False, flags=re.IGNORECASE)))

combined_diagnoses['neuro'] = (
    (combined_diagnoses['surgical'] == 0) &
    (combined_diagnoses['diagnosis'].str.contains(
        re_neuro, na=False, flags=re.IGNORECASE)))
## exclude trauma cases
combined_diagnoses['neuro'] &= (
    ~(combined_diagnoses['diagnosis'].str.contains(
        re_trauma_surg, na=False, flags=re.IGNORECASE)))

combined_diagnoses['cardio'] = (
    (combined_diagnoses['surgical'] == 0) &
    (combined_diagnoses['diagnosis'].str.contains(
        re_cardio, na=False, flags=re.IGNORECASE)))

################################################################################
## Save file

combined_diagnoses.to_csv(
    additional_file_path + 'combined_diagnoses.csv', index=False)

################################################################################

def counts_fun(df):
    cols = ['diagnosis_group', 'diagnosis']
    counts = df.groupby(cols).size().to_frame('n_admissions')
    return counts.sort_values('n_admissions', ascending=False).reset_index()

cardiosurg_counts = counts_fun(
    combined_diagnoses.loc[combined_diagnoses['cardiosurg']])
respfailure_counts = counts_fun(
    combined_diagnoses.loc[combined_diagnoses['respfailure']])
neurosurg_counts = counts_fun(
    combined_diagnoses.loc[combined_diagnoses['neurosurg']])
gisurg_counts = counts_fun(
    combined_diagnoses.loc[combined_diagnoses['gisurg']])
cardiacarrest_counts = counts_fun(
    combined_diagnoses.loc[combined_diagnoses['cardiacarrest']])
vascsurg_counts = counts_fun(
    combined_diagnoses.loc[combined_diagnoses['vascsurg']])
trauma_counts = counts_fun(combined_diagnoses.loc[combined_diagnoses['trauma']])
neuro_counts = counts_fun(combined_diagnoses.loc[combined_diagnoses['neuro']])
cardio_counts = counts_fun(combined_diagnoses.loc[combined_diagnoses['cardio']])

################################################################################
## Build reason for admission table

icu_ind = combined_diagnoses['location'].str.contains('IC')
mcu_ind = combined_diagnoses['location'].str.contains('^MC$')

## Number of patients
no_patients_total = combined_diagnoses['patientid'].nunique()
no_patients_ICU = combined_diagnoses.loc[icu_ind, 'patientid'].nunique()
no_patients_MCU = combined_diagnoses.loc[mcu_ind, 'patientid'].nunique()

## Number of admissions
no_admissions_total = combined_diagnoses['admissionid'].nunique()
no_admissions_ICU = combined_diagnoses.loc[icu_ind, 'admissionid'].nunique()
no_admissions_MCU = combined_diagnoses.loc[mcu_ind, 'admissionid'].nunique()

## Surgical admissions
surgical_total = (combined_diagnoses['surgical'] == 1).sum()
surgical_pct_total = 100 * surgical_total / no_admissions_total
surgical_ICU = ((combined_diagnoses['surgical'] == 1) & icu_ind).sum()
surgical_pct_ICU = 100 * surgical_ICU / no_admissions_ICU
surgical_MCU = ((combined_diagnoses['surgical'] == 1) & mcu_ind).sum()
surgical_pct_MCU = 100 * surgical_MCU / no_admissions_MCU

## Medical admissions
medical_total = (combined_diagnoses['surgical'] == 0).sum()
medical_pct_total = 100 * medical_total / no_admissions_total
medical_ICU = ((combined_diagnoses['surgical'] == 0) & icu_ind).sum()
medical_pct_ICU = 100 * medical_ICU / no_admissions_ICU
medical_MCU = ((combined_diagnoses['surgical'] == 0) & mcu_ind).sum()
medical_pct_MCU = 100 * medical_MCU / no_admissions_MCU

## Surgical/medical unknown
surgical_na_total = combined_diagnoses['surgical'].isna().sum()
surgical_na_pct_total = (100 * surgical_na_total / no_admissions_total)
surgical_na_ICU = (combined_diagnoses['surgical'].isna() & icu_ind).sum()
surgical_na_pct_ICU = 100 * surgical_na_ICU / no_admissions_ICU
surgical_na_MCU = (combined_diagnoses['surgical'].isna() & mcu_ind).sum()
surgical_na_pct_MCU = 100 * surgical_na_MCU / no_admissions_MCU

## Urgency
urgency_total = (combined_diagnoses['urgency'] == 1).sum()
urgency_pct_total = 100 * urgency_total / no_admissions_total
urgency_ICU = ((combined_diagnoses['urgency'] == 1) & icu_ind).sum()
urgency_pct_ICU = 100 * urgency_ICU / no_admissions_ICU
urgency_MCU = ((combined_diagnoses['urgency'] == 1) & mcu_ind).sum()
urgency_pct_MCU = 100 * urgency_MCU / no_admissions_MCU

## Cardiac surgery
cardiosurg_total = (combined_diagnoses['cardiosurg']).sum()
cardiosurg_pct_total = 100 * cardiosurg_total / no_admissions_total
cardiosurg_ICU = (combined_diagnoses['cardiosurg'] & icu_ind).sum()
cardiosurg_pct_ICU = 100 * cardiosurg_ICU / no_admissions_ICU
cardiosurg_MCU = (combined_diagnoses['cardiosurg'] & mcu_ind).sum()
cardiosurg_pct_MCU = 100 * cardiosurg_MCU / no_admissions_MCU

## Neurosurgery
neurosurg_total = (combined_diagnoses['neurosurg']).sum()
neurosurg_pct_total = 100 * neurosurg_total / no_admissions_total
neurosurg_ICU = (combined_diagnoses['neurosurg'] & icu_ind).sum()
neurosurg_pct_ICU = 100 * neurosurg_ICU / no_admissions_ICU
neurosurg_MCU = (combined_diagnoses['neurosurg'] & mcu_ind).sum()
neurosurg_pct_MCU = 100 * neurosurg_MCU / no_admissions_MCU

## Gastro-intestinal surgery
gisurg_total = (combined_diagnoses['gisurg']).sum()
gisurg_pct_total = 100 * gisurg_total / no_admissions_total
gisurg_ICU = (combined_diagnoses['gisurg'] & icu_ind).sum()
gisurg_pct_ICU = 100 * gisurg_ICU / no_admissions_ICU
gisurg_MCU = (combined_diagnoses['gisurg'] & mcu_ind).sum()
gisurg_pct_MCU = 100 * gisurg_MCU / no_admissions_MCU

## Cardiac arrest
cardiacarrest_total = (combined_diagnoses['cardiacarrest']).sum()
cardiacarrest_pct_total = 100 * cardiacarrest_total / no_admissions_total
cardiacarrest_ICU = (combined_diagnoses['cardiacarrest'] & icu_ind).sum()
cardiacarrest_pct_ICU = 100 * cardiacarrest_ICU / no_admissions_ICU
cardiacarrest_MCU = (combined_diagnoses['cardiacarrest'] & mcu_ind).sum()
cardiacarrest_pct_MCU = 100 * gisurg_MCU / no_admissions_MCU

## Vascular surgery
vascsurg_total = (combined_diagnoses['vascsurg']).sum()
vascsurg_pct_total = 100 * vascsurg_total / no_admissions_total
vascsurg_ICU = (combined_diagnoses['vascsurg'] & icu_ind).sum()
vascsurg_pct_ICU = 100 * vascsurg_ICU / no_admissions_ICU
vascsurg_MCU = (combined_diagnoses['vascsurg'] & mcu_ind).sum()
vascsurg_pct_MCU = 100 * vascsurg_MCU / no_admissions_MCU

## Trauma
trauma_total = (combined_diagnoses['trauma']).sum()
trauma_pct_total = 100 * trauma_total / no_admissions_total
trauma_ICU = (combined_diagnoses['trauma'] & icu_ind).sum()
trauma_pct_ICU = 100 * trauma_ICU / no_admissions_ICU
trauma_MCU = (combined_diagnoses['trauma'] & mcu_ind).sum()
trauma_pct_MCU = 100 * trauma_MCU / no_admissions_MCU

## Sepsis
sepsis_total = (combined_diagnoses['sepsis']).sum()
sepsis_pct_total = 100 * sepsis_total / no_admissions_total
sepsis_ICU = (combined_diagnoses['sepsis'] & icu_ind).sum()
sepsis_pct_ICU = 100 * sepsis_ICU / no_admissions_ICU
sepsis_MCU = (combined_diagnoses['sepsis'] & mcu_ind).sum()
sepsis_pct_MCU = 100 * sepsis_MCU / no_admissions_MCU

## Respiratory failure
respfailure_total = (combined_diagnoses['respfailure']).sum()
respfailure_pct_total = 100 * respfailure_total / no_admissions_total
respfailure_ICU = (combined_diagnoses['respfailure'] & icu_ind).sum()
respfailure_pct_ICU = 100 * respfailure_ICU / no_admissions_ICU
respfailure_MCU = (combined_diagnoses['respfailure'] & mcu_ind).sum()
respfailure_pct_MCU = 100 * respfailure_MCU / no_admissions_MCU

## Cardio
cardio_total = (combined_diagnoses['cardio']).sum()
cardio_pct_total = 100 * cardio_total / no_admissions_total
cardio_ICU = (combined_diagnoses['cardio'] & icu_ind).sum()
cardio_pct_ICU = 100 * cardio_ICU / no_admissions_ICU
cardio_MCU = (combined_diagnoses['cardio'] & mcu_ind).sum()
cardio_pct_MCU = 100 * cardio_MCU / no_admissions_MCU

## Neuro
neuro_total = (combined_diagnoses['neuro']).sum()
neuro_pct_total = 100 * neuro_total / no_admissions_total
neuro_ICU = (combined_diagnoses['neuro'] & icu_ind).sum()
neuro_pct_ICU = 100 * neuro_ICU / no_admissions_ICU
neuro_MCU = (combined_diagnoses['neuro'] & mcu_ind).sum()
neuro_pct_MCU = 100 * neuro_MCU / no_admissions_MCU

reason_cols = ['Distinct patients', 'ICU admissions', 'Surgical admissions']
reason_cols += ['Medical admissions', 'Surgical/medical unknown']
reason_cols += ['Urgent admissions', 'Cardiothoracic', 'Sepsis']
reason_cols += ['Respiratory failure', 'Neurosurgery', 'Trauma']
reason_cols += ['Gastro-intestinal', 'Vascular surgery', 'Cardiac arrest']
reason_cols += ['Neurological', 'Cardiac disorder']
reason_for_admission = pd.DataFrame(
    index=reason_cols, columns=['Total', 'ICU', 'MCU'])

reason_for_admission.loc['Distinct patients', 'Total'] = no_patients_total
reason_for_admission.loc['Distinct patients', 'ICU'] = no_patients_ICU
reason_for_admission.loc['Distinct patients', 'MCU'] = no_patients_MCU

reason_for_admission.loc['ICU admissions', 'Total'] = no_admissions_total
reason_for_admission.loc['ICU admissions', 'ICU'] = no_admissions_ICU
reason_for_admission.loc['ICU admissions', 'MCU'] = no_admissions_MCU

reason_for_admission.loc['Surgical admissions', 'Total'] = surgical_total
reason_for_admission.loc['Surgical admissions', 'ICU'] = surgical_ICU
reason_for_admission.loc['Surgical admissions', 'MCU'] = surgical_MCU

reason_for_admission.loc['Medical admissions', 'Total'] = medical_total
reason_for_admission.loc['Medical admissions', 'ICU'] = medical_ICU
reason_for_admission.loc['Medical admissions', 'MCU'] = medical_MCU

reason_for_admission.loc['Surgical/medical unknown', 'Total'] = surgical_na_total
reason_for_admission.loc['Surgical/medical unknown', 'ICU'] = surgical_na_ICU
reason_for_admission.loc['Surgical/medical unknown', 'MCU'] = surgical_na_MCU

reason_for_admission.loc['Urgent admissions', 'Total'] = urgency_total
reason_for_admission.loc['Urgent admissions', 'ICU'] = urgency_ICU
reason_for_admission.loc['Urgent admissions', 'MCU'] = urgency_MCU

reason_for_admission.loc['Cardiothoracic', 'Total'] = cardiosurg_total
reason_for_admission.loc['Cardiothoracic', 'ICU'] = cardiosurg_ICU
reason_for_admission.loc['Cardiothoracic', 'MCU'] = cardiosurg_MCU

reason_for_admission.loc['Sepsis', 'Total'] = sepsis_total
reason_for_admission.loc['Sepsis', 'ICU'] = sepsis_ICU
reason_for_admission.loc['Sepsis', 'MCU'] = sepsis_MCU

reason_for_admission.loc['Respiratory failure', 'Total'] = respfailure_total
reason_for_admission.loc['Respiratory failure', 'ICU'] = respfailure_ICU
reason_for_admission.loc['Respiratory failure', 'MCU'] = respfailure_MCU

reason_for_admission.loc['Neurosurgery', 'Total'] = neurosurg_total
reason_for_admission.loc['Neurosurgery', 'ICU'] = neurosurg_ICU
reason_for_admission.loc['Neurosurgery', 'MCU'] = neurosurg_MCU

reason_for_admission.loc['Trauma', 'Total'] = trauma_total
reason_for_admission.loc['Trauma', 'ICU'] = trauma_ICU
reason_for_admission.loc['Trauma', 'MCU'] = trauma_MCU

reason_for_admission.loc['Gastro-intestinal', 'Total'] = gisurg_total
reason_for_admission.loc['Gastro-intestinal', 'ICU'] = gisurg_ICU
reason_for_admission.loc['Gastro-intestinal', 'MCU'] = gisurg_MCU

reason_for_admission.loc['Vascular surgery', 'Total'] = vascsurg_total
reason_for_admission.loc['Vascular surgery', 'ICU'] = vascsurg_ICU
reason_for_admission.loc['Vascular surgery', 'MCU'] = vascsurg_MCU

reason_for_admission.loc['Cardiac arrest', 'Total'] = cardiacarrest_total
reason_for_admission.loc['Cardiac arrest', 'ICU'] = cardiacarrest_ICU
reason_for_admission.loc['Cardiac arrest', 'MCU'] = cardiacarrest_MCU

reason_for_admission.loc['Neurological', 'Total'] = neuro_total
reason_for_admission.loc['Neurological', 'ICU'] = neuro_ICU
reason_for_admission.loc['Neurological', 'MCU'] = neuro_MCU

reason_for_admission.loc['Cardiac disorder', 'Total'] = cardio_total
reason_for_admission.loc['Cardiac disorder', 'ICU'] = cardio_ICU
reason_for_admission.loc['Cardiac disorder', 'MCU'] = cardio_MCU

########################
diagnoses_sql = """
WITH diagnosis_groups AS (
    SELECT admissionid,
        item,
        CASE
            WHEN itemid IN (
                18669, --NICE APACHEII diagnosen
                18671 --NICE APACHEIV diagnosen
            )
            THEN split_part(value, ' - ', 1)
            -- 'e.g. 'Non-operative cardiovascular - Anaphylaxis' -> Non-operative cardiovascular
            ELSE value
        END as diagnosis_group,
        valueid as diagnosis_group_id,
        ROW_NUMBER() OVER(PARTITION BY admissionid
        ORDER BY
            CASE --prefer NICE > APACHE IV > II > D
                WHEN itemid = 18671 THEN 6 --NICE APACHEIV diagnosen
                WHEN itemid = 18669 THEN 5 --NICE APACHEII diagnosen
                WHEN itemid BETWEEN 16998 AND 17017 THEN 4 --APACHE IV diagnosis
                WHEN itemid BETWEEN 18589 AND 18602 THEN 3 --APACHE II diagnosis
                WHEN itemid BETWEEN 13116 AND 13145 THEN 2 --D diagnosis ICU
                WHEN itemid BETWEEN 16642 AND 16673 THEN 1 --DMC diagnosis Medium Care
            END DESC,
        measuredat DESC) AS rownum
    FROM listitems
    WHERE itemid IN (
        --MAIN GROUP - LEVEL 0
        13110, --D_Hoofdgroep
        16651, --DMC_Hoofdgroep, Medium Care

        18588, --Apache II Hoofdgroep
        16997, --APACHE IV Groepen

        18669, --NICE APACHEII diagnosen
        18671 --NICE APACHEIV diagnosen
    )
),diagnosis_subgroups AS (
    SELECT admissionid,
        item,
        value as diagnosis_subgroup,
        valueid as diagnosis_subgroup_id,
        ROW_NUMBER() OVER(PARTITION BY admissionid
        ORDER BY measuredat DESC) AS rownum
    FROM listitems
    WHERE itemid IN (
        --SUB GROUP - LEVEL 1
        13111, --D_Subgroep_Thoraxchirurgie
        16669, --DMC_Subgroep_Thoraxchirurgie
        13112, --D_Subgroep_Algemene chirurgie
        16665, --DMC_Subgroep_Algemene chirurgie
        13113, --D_Subgroep_Neurochirurgie
        16667, --DMC_Subgroep_Neurochirurgie
        13114, --D_Subgroep_Neurologie
        16668, --DMC_Subgroep_Neurologie
        13115, --D_Subgroep_Interne geneeskunde
        16666 --DMC_Subgroep_Interne geneeskunde
    )
), diagnoses AS (
    SELECT admissionid,
        item,
        CASE
            WHEN itemid IN (
                18669, --NICE APACHEII diagnosen
                18671 --NICE APACHEIV diagnosen
            )
            THEN split_part(value, ' - ', 2)
            -- 'e.g. 'Non-operative cardiovascular - Anaphylaxis' -> Anaphylaxis
            ELSE value
        END as diagnosis,
        CASE
            WHEN itemid IN (
                --SURGICAL
                13116, --D_Thoraxchirurgie_CABG en Klepchirurgie
                16671, --DMC_Thoraxchirurgie_CABG en Klepchirurgie
                13117, --D_Thoraxchirurgie_Cardio anders
                16672, --DMC_Thoraxchirurgie_Cardio anders
                13118, --D_Thoraxchirurgie_Aorta chirurgie
                16670, --DMC_Thoraxchirurgie_Aorta chirurgie
                13119, --D_Thoraxchirurgie_Pulmonale chirurgie
                16673, --DMC_Thoraxchirurgie_Pulmonale chirurgie

                --Not surgical: 13141, --D_Algemene chirurgie_Algemeen
                --Not surgical: 16642, --DMC_Algemene chirurgie_Algemeen
                13121, --D_Algemene chirurgie_Buikchirurgie
                16643, --DMC_Algemene chirurgie_Buikchirurgie
                13123, --D_Algemene chirurgie_Endocrinologische chirurgie
                16644, --DMC_Algemene chirurgie_Endocrinologische chirurgie
                13145, --D_Algemene chirurgie_KNO/Overige
                16645, --DMC_Algemene chirurgie_KNO/Overige
                13125, --D_Algemene chirurgie_Orthopedische chirurgie
                16646, --DMC_Algemene chirurgie_Orthopedische chirurgie
                13122, --D_Algemene chirurgie_Transplantatie chirurgie
                16647, --DMC_Algemene chirurgie_Transplantatie chirurgie
                13124, --D_Algemene chirurgie_Trauma
                16648, --DMC_Algemene chirurgie_Trauma
                13126, --D_Algemene chirurgie_Urogenitaal
                16649, --DMC_Algemene chirurgie_Urogenitaal
                13120, --D_Algemene chirurgie_Vaatchirurgie
                16650, --DMC_Algemene chirurgie_Vaatchirurgie

                13128, --D_Neurochirurgie _Vasculair chirurgisch
                16661, --DMC_Neurochirurgie _Vasculair chirurgisch
                13129, --D_Neurochirurgie _Tumor chirurgie
                16660, --DMC_Neurochirurgie _Tumor chirurgie
                13130, --D_Neurochirurgie_Overige
                16662, --DMC_Neurochirurgie_Overige

                18596, --Apache II Operatief  Gastr-intenstinaal
                18597, --Apache II Operatief Cardiovasculair
                18598, --Apache II Operatief Hematologisch
                18599, --Apache II Operatief Metabolisme
                18600, --Apache II Operatief Neurologisch
                18601, --Apache II Operatief Renaal
                18602, --Apache II Operatief Respiratoir

                17008, --APACHEIV Post-operative cardiovascular
                17009, --APACHEIV Post-operative gastro-intestinal
                17010, --APACHEIV Post-operative genitourinary
                17011, --APACHEIV Post-operative hematology
                17012, --APACHEIV Post-operative metabolic
                17013, --APACHEIV Post-operative musculoskeletal /skin
                17014, --APACHEIV Post-operative neurologic
                17015, --APACHEIV Post-operative respiratory
                17016, --APACHEIV Post-operative transplant
                17017 --APACHEIV Post-operative trauma

            ) THEN 1
            WHEN itemid = 18669 AND valueid BETWEEN 1 AND 26 THEN 1 --NICE APACHEII diagnosen
            WHEN itemid = 18671 AND valueid BETWEEN 222 AND 452 THEN 1 --NICE APACHEIV diagnosen
            ELSE 0
        END AS surgical,
        valueid as diagnosis_id,
        CASE
                WHEN itemid = 18671 THEN 'NICE APACHE IV'
                WHEN itemid = 18669 THEN 'NICE APACHE II'
                WHEN itemid BETWEEN 16998 AND 17017 THEN 'APACHE IV'
                WHEN itemid BETWEEN 18589 AND 18602 THEN 'APACHE II'
                WHEN itemid BETWEEN 13116 AND 13145 THEN 'Legacy ICU'
                WHEN itemid BETWEEN 16642 AND 16673 THEN 'Legacy MCU'
        END AS diagnosis_type,
        ROW_NUMBER() OVER(PARTITION BY admissionid
        ORDER BY
            CASE --prefer NICE > APACHE IV > II > D
                WHEN itemid = 18671 THEN 6 --NICE APACHEIV diagnosen
                WHEN itemid = 18669 THEN 5 --NICE APACHEII diagnosen
                WHEN itemid BETWEEN 16998 AND 17017 THEN 4 --APACHE IV diagnosis
                WHEN itemid BETWEEN 18589 AND 18602 THEN 3 --APACHE II diagnosis
                WHEN itemid BETWEEN 13116 AND 13145 THEN 2 --D diagnosis ICU
                WHEN itemid BETWEEN 16642 AND 16673 THEN 1 --DMC diagnosis Medium Care
            END DESC,
            measuredat DESC) AS rownum
    FROM listitems
    WHERE itemid IN (
        -- Diagnosis - LEVEL 2
        --SURGICAL
        13116, --D_Thoraxchirurgie_CABG en Klepchirurgie
        16671, --DMC_Thoraxchirurgie_CABG en Klepchirurgie
        13117, --D_Thoraxchirurgie_Cardio anders
        16672, --DMC_Thoraxchirurgie_Cardio anders
        13118, --D_Thoraxchirurgie_Aorta chirurgie
        16670, --DMC_Thoraxchirurgie_Aorta chirurgie
        13119, --D_Thoraxchirurgie_Pulmonale chirurgie
        16673, --DMC_Thoraxchirurgie_Pulmonale chirurgie

        13141, --D_Algemene chirurgie_Algemeen
        16642, --DMC_Algemene chirurgie_Algemeen
        13121, --D_Algemene chirurgie_Buikchirurgie
        16643, --DMC_Algemene chirurgie_Buikchirurgie
        13123, --D_Algemene chirurgie_Endocrinologische chirurgie
        16644, --DMC_Algemene chirurgie_Endocrinologische chirurgie
        13145, --D_Algemene chirurgie_KNO/Overige
        16645, --DMC_Algemene chirurgie_KNO/Overige
        13125, --D_Algemene chirurgie_Orthopedische chirurgie
        16646, --DMC_Algemene chirurgie_Orthopedische chirurgie
        13122, --D_Algemene chirurgie_Transplantatie chirurgie
        16647, --DMC_Algemene chirurgie_Transplantatie chirurgie
        13124, --D_Algemene chirurgie_Trauma
        16648, --DMC_Algemene chirurgie_Trauma
        13126, --D_Algemene chirurgie_Urogenitaal
        16649, --DMC_Algemene chirurgie_Urogenitaal
        13120, --D_Algemene chirurgie_Vaatchirurgie
        16650, --DMC_Algemene chirurgie_Vaatchirurgie

        13128, --D_Neurochirurgie _Vasculair chirurgisch
        16661, --DMC_Neurochirurgie _Vasculair chirurgisch
        13129, --D_Neurochirurgie _Tumor chirurgie
        16660, --DMC_Neurochirurgie _Tumor chirurgie
        13130, --D_Neurochirurgie_Overige
        16662, --DMC_Neurochirurgie_Overige

        18596, --Apache II Operatief  Gastr-intenstinaal
        18597, --Apache II Operatief Cardiovasculair
        18598, --Apache II Operatief Hematologisch
        18599, --Apache II Operatief Metabolisme
        18600, --Apache II Operatief Neurologisch
        18601, --Apache II Operatief Renaal
        18602, --Apache II Operatief Respiratoir

        17008, --APACHEIV Post-operative cardiovascular
        17009, --APACHEIV Post-operative gastro-intestinal
        17010, --APACHEIV Post-operative genitourinary
        17011, --APACHEIV Post-operative hematology
        17012, --APACHEIV Post-operative metabolic
        17013, --APACHEIV Post-operative musculoskeletal /skin
        17014, --APACHEIV Post-operative neurologic
        17015, --APACHEIV Post-operative respiratory
        17016, --APACHEIV Post-operative transplant
        17017, --APACHEIV Post-operative trauma

        --MEDICAL
        13133, --D_Interne Geneeskunde_Cardiovasculair
        16653, --DMC_Interne Geneeskunde_Cardiovasculair
        13134, --D_Interne Geneeskunde_Pulmonaal
        16658, --DMC_Interne Geneeskunde_Pulmonaal
        13135, --D_Interne Geneeskunde_Abdominaal
        16652, --DMC_Interne Geneeskunde_Abdominaal
        13136, --D_Interne Geneeskunde_Infectieziekten
        16655, --DMC_Interne Geneeskunde_Infectieziekten
        13137, --D_Interne Geneeskunde_Metabool
        16656, --DMC_Interne Geneeskunde_Metabool
        13138, --D_Interne Geneeskunde_Renaal
        16659, --DMC_Interne Geneeskunde_Renaal
        13139, --D_Interne Geneeskunde_Hematologisch
        16654, --DMC_Interne Geneeskunde_Hematologisch
        13140, --D_Interne Geneeskunde_Overige
        16657, --DMC_Interne Geneeskunde_Overige

        13131, --D_Neurologie_Vasculair neurologisch
        16664, --DMC_Neurologie_Vasculair neurologisch
        13132, --D_Neurologie_Overige
        16663, --DMC_Neurologie_Overige
        13127, --D_KNO/Overige

        18589, --Apache II Non-Operatief Cardiovasculair
        18590, --Apache II Non-Operatief Gastro-intestinaal
        18591, --Apache II Non-Operatief Hematologisch
        18592, --Apache II Non-Operatief Metabolisme
        18593, --Apache II Non-Operatief Neurologisch
        18594, --Apache II Non-Operatief Renaal
        18595, --Apache II Non-Operatief Respiratoir

        16998, --APACHE IV Non-operative cardiovascular
        16999, --APACHE IV Non-operative Gastro-intestinal
        17000, --APACHE IV Non-operative genitourinary
        17001, --APACHEIV  Non-operative haematological
        17002, --APACHEIV  Non-operative metabolic
        17003, --APACHEIV Non-operative musculo-skeletal
        17004, --APACHEIV Non-operative neurologic
        17005, --APACHEIV Non-operative respiratory
        17006, --APACHEIV Non-operative transplant
        17007, --APACHEIV Non-operative trauma

        --NICE: surgical/medical combined in same parameter
        18669, --NICE APACHEII diagnosen
        18671 --NICE APACHEIV diagnosen
    )
), sepsis AS (
    SELECT
        admissionid,
        CASE valueid
            WHEN 1 THEN 1 --'Ja'
            WHEN 2 THEN 0 --'Nee'
        END as sepsis_at_admission,
        ROW_NUMBER() OVER(
            PARTITION BY
                admissionid
            ORDER BY
                measuredat DESC) AS rownum
    FROM listitems
    WHERE
        itemid = 15808
), sepsis_antibiotics AS ( --non prophylactic antibiotics
    SELECT
        admissionid,
        CASE
            WHEN COUNT(*) > 0 THEN 1
            ELSE 0
        END AS sepsis_antibiotics_bool,
        STRING_AGG(DISTINCT item, '; ') AS sepsis_antibiotics_given
    FROM drugitems
    WHERE
        itemid IN (
            6834, --Amikacine (Amukin)
            6847, --Amoxicilline (Clamoxyl/Flemoxin)
            6871, --Benzylpenicilline (Penicilline)
            6917, --Ceftazidim (Fortum)
            --6919, --Cefotaxim (Claforan) -> prophylaxis
            6948, --Ciprofloxacine (Ciproxin)
            6953, --Rifampicine (Rifadin)
            6958, --Clindamycine (Dalacin)
            7044, --Tobramycine (Obracin)
            --7064, --Vancomycine -> prophylaxis for valve surgery
            7123, --Imipenem (Tienam)
            7185, --Doxycycline (Vibramycine)
            --7187, --Metronidazol (Flagyl) -> often used for GI surgical prophylaxis
            --7208, --Erythromycine (Erythrocine) -> often used for gastroparesis
            7227, --Flucloxacilline (Stafoxil/Floxapen)
            7231, --Fluconazol (Diflucan)
            7232, --Ganciclovir (Cymevene)
            7233, --Flucytosine (Ancotil)
            7235, --Gentamicine (Garamycin)
            7243, --Foscarnet trinatrium (Foscavir)
            7450, --Amfotericine B (Fungizone)
            --7504, --X nader te bepalen --non-stock medication
            8127, --Meropenem (Meronem)
            8229, --Myambutol (ethambutol)
            8374, --Kinine dihydrocloride
            --8375, --Immunoglobuline (Nanogam) -> not anbiotic
            --8394, --Co-Trimoxazol (Bactrimel) -> often prophylactic (unless high dose)
            8547, --Voriconazol(VFEND)
            --9029, --Amoxicilline/Clavulaanzuur (Augmentin) -> often used for ENT surgical prophylaxis
            9030, --Aztreonam (Azactam)
            9047, --Chlooramfenicol
            --9075, --Fusidinezuur (Fucidin) -> prophylaxis
            9128, --Piperacilline (Pipcil)
            9133, --Ceftriaxon (Rocephin)
            --9151, --Cefuroxim (Zinacef) -> often used for GI/transplant surgical prophylaxis
            --9152, --Cefazoline (Kefzol) -> prophylaxis for cardiac surgery
            9458, --Caspofungine
            9542, --Itraconazol (Trisporal)
            --9602, --Tetanusimmunoglobuline -> prophylaxis/not antibiotic
            12398, --Levofloxacine (Tavanic)
            12772, --Amfotericine B lipidencomplex  (Abelcet)
            15739, --Ecalta (Anidulafungine)
            16367, --Research Anidulafungin/placebo
            16368, --Research Caspofungin/placebo
            18675, --Amfotericine B in liposomen (Ambisome )
            19137, --Linezolid (Zyvoxid)
            19764, --Tigecycline (Tygacil)
            19773, --Daptomycine (Cubicin)
            20175 --Colistine
        )
        AND start < 24*60*60*1000 --within 24 hours (to correct for antibiotics administered before ICU)
    GROUP BY admissionid
), other_antibiotics AS ( --'prophylactic' antibiotics that may be used for sepsis
    SELECT
        admissionid,
        CASE
            WHEN COUNT(*) > 0 THEN 1
            ELSE 0
        END AS other_antibiotics_bool,
        STRING_AGG(DISTINCT item, '; ') AS other_antibiotics_given
    FROM drugitems
    WHERE
        itemid IN (
            7064, --Vancomycine -> prophylaxis for valve surgery
            7187, --Metronidazol (Flagyl) -> often used for GI surgical prophylaxis
            8394, --Co-Trimoxazol (Bactrimel) -> often prophylactic (unless high dose)
            9029, --Amoxicilline/Clavulaanzuur (Augmentin) -> often used for ENT surgical prophylaxis
            9151, --Cefuroxim (Zinacef) -> often used for GI surgical prophylaxis
            9152 --Cefazoline (Kefzol) -> prophylaxis
        )
        AND start < 24*60*60*1000 --within 24 hours (to correct for antibiotics administered before ICU)
    GROUP BY admissionid
), cultures AS (
    SELECT
        admissionid,
        CASE
            WHEN COUNT(*) > 0 THEN 1
            ELSE 0
        END AS sepsis_cultures_bool,
        STRING_AGG(DISTINCT item, '; ') AS sepsis_cultures_drawn
    FROM procedureorderitems
    WHERE
        itemid IN (
        --8097, --Sputumkweek afnemen -> often used routinely
        --8418, --Urinekweek afnemen
        --8588, --MRSA kweken afnemen
        9189, --Bloedkweken afnemen
        9190, --Cathetertipkweek afnemen
        --9191, --Drainvochtkweek afnemen
        --9192, --Faeceskweek afnemen -> Clostridium
        --9193, --X-Kweek nader te bepalen
        --9194, --Liquorkweek afnemen
        --9195, --Neuskweek afnemen
        --9197, --Perineumkweek afnemen -> often used routinely
        -9198, --Rectumkweek afnemen -> often used routinely
        9200, --Wondkweek afnemen
        9202, --Ascitesvochtkweek afnemen
        --9203, --Keelkweek afnemen -> often used routinely
        --9204, --SDD-kweken afnemen -> often used routinely
        9205 --Legionella sneltest (urine)
        --1302, --SDD Inventarisatiekweken afnemen -> often used routinely
        --19663, --Research Neuskweek COUrSe
        --19664, --Research Sputumkweek COUrSe
        )
        AND registeredat < 6*60*60*1000 --within 6 hours
    GROUP BY admissionid
)
SELECT
    admissions.*
    , diagnosis_type
    , diagnosis, diagnosis_id
    , diagnosis_subgroup
    , diagnosis_subgroup_id
    , diagnosis_group
    , diagnosis_group_id
    , surgical
    , sepsis_at_admission
    , sepsis_antibiotics_bool
    , sepsis_antibiotics_given
    , other_antibiotics_bool
    , other_antibiotics_given
    , sepsis_cultures_bool
    , sepsis_cultures_drawn
FROM admissions
LEFT JOIN diagnoses on admissions.admissionid = diagnoses.admissionid
LEFT JOIN diagnosis_subgroups on admissions.admissionid = diagnosis_subgroups.admissionid
LEFT JOIN diagnosis_groups on admissions.admissionid = diagnosis_groups.admissionid
LEFT JOIN sepsis on admissions.admissionid = sepsis.admissionid
LEFT JOIN sepsis_antibiotics on admissions.admissionid = sepsis_antibiotics.admissionid
LEFT JOIN other_antibiotics on admissions.admissionid = other_antibiotics.admissionid
LEFT JOIN cultures on admissions.admissionid = cultures.admissionid
WHERE --only last updated record
    (diagnoses.rownum = 1 OR diagnoses.rownum IS NULL) AND
    (diagnosis_subgroups.rownum = 1 OR diagnosis_subgroups.rownum IS NULL) AND
    (diagnosis_groups.rownum = 1 OR diagnosis_groups.rownum IS NULL) AND
    (sepsis.rownum = 1 OR sepsis.rownum IS NULL)
;
"""
