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
import argparse
import amsterdamumcdb as adb

import matplotlib.pyplot as plt
import lifelines
from matplotlib.backends.backend_pdf import PdfPages
from palettable.colorbrewer import qualitative as ql
from palettable.colorbrewer import diverging as dv
from palettable.colorbrewer import sequential as sq
from matplotlib.patches import Rectangle, Polygon
from matplotlib.collections import PatchCollection

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
        default='../../data/additional_files/figures/',
        help='''File path to the directory that will contain the outputs of
            this script, including all figures (default: %(default)s)''',
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

# inputs.additional_file_path = '../../data/additional_files_t0_v4_161022/'
# inputs.output_file_path = '../../data/additional_files_t0_v4_161022/figures/'

###############################################################################
# Load intermediate files

dictionary = adb.get_dictionary()

sofa = pd.read_csv(inputs.additional_file_path + 'sofa.csv')
# Including subscores
sepsis = pd.read_csv(inputs.additional_file_path + 'sepsis_all.csv')
drugitems_abx_all = pd.read_csv(
    inputs.additional_file_path + 'drugitems_abx_all.csv')

descriptors_all = pd.read_csv(inputs.additional_file_path + 'descriptors.csv')
descriptors_all = descriptors_all.loc[(descriptors_all['discard'] == 0)]
descriptors_icu = descriptors_all.loc[(descriptors_all['location_IC'] == 1)]
descriptors_mcu = descriptors_all.loc[(descriptors_all['location_IC'] == 0)]

abx_usage = pd.read_csv(inputs.additional_file_path + 'abx_usage.csv')
abx_usage.rename(columns={'start_time': 'time'}, inplace=True)
abx_usage_icu = abx_usage.loc[(abx_usage['location_IC'] == 1)]
abx_usage_mcu = abx_usage.loc[(abx_usage['location_IC'] == 0)]

abx_np_usage = pd.read_csv(
    inputs.additional_file_path + 'abx_non_prophylactic_usage.csv')
abx_np_usage.rename(columns={'start_time': 'time'}, inplace=True)
abx_np_usage['location_IC'] = abx_np_usage['location_IC'].fillna(False)
abx_np_usage_icu = abx_np_usage.loc[(abx_np_usage['location_IC'] == 1)]
abx_np_usage_mcu = abx_np_usage.loc[(abx_np_usage['location_IC'] == 0)]

###############################################################################
# Prep for figures: turn colors to HEX codes

Dark2_6_hex = [None] * len(ql.Dark2_6.mpl_colors)
RdYlBu_4_hex = [None] * len(dv.RdYlBu_4.mpl_colors)
for ii in range(len(ql.Dark2_6.mpl_colors)):
    c = [int(256 * x) for x in ql.Dark2_6.mpl_colors[ii]]
    Dark2_6_hex[ii] = '#{0:02x}{1:02x}{2:02x}'.format(c[0], c[1], c[2])
for ii in range(len(dv.RdYlBu_4.mpl_colors)):
    c = [int(256 * x) for x in dv.RdYlBu_4.mpl_colors[ii]]
    RdYlBu_4_hex[ii] = '#{0:02x}{1:02x}{2:02x}'.format(c[0], c[1], c[2])

###############################################################################
# Prep for figures: intermediary functions


def counts_fun(df):
    area_plot = df.groupby('time').agg(
            no_abx=pd.NamedAgg(column='no_abx', aggfunc='sum'),
            abx=pd.NamedAgg(column='abx', aggfunc='sum'),
            sepsis_episode=pd.NamedAgg(column='sepsis_episode', aggfunc='sum'),
            septic_shock=pd.NamedAgg(column='septic_shock', aggfunc='sum'),
            discharged=pd.NamedAgg(column='discharged', aggfunc='sum'),
            died=pd.NamedAgg(column='died', aggfunc='sum'),
        ).reset_index()
    area_plot.rename(columns={
        'time': 'Time since admission (d)',
        'no_abx': 'Not on antibiotics',
        'abx': 'On antibiotics (no sepsis)',
        'sepsis_episode': 'Sepsis without shock',
        'septic_shock': 'Septic shock',
        'discharged': 'Discharged',
        'died': 'Died'}, inplace=True)
    return area_plot


def kmf_plot(df, kmf, ax, color, ls='-', label=None, arc=False, min_shape=50):
    if df.shape[0] > min_shape:
        kmf.fit(
            df['lengthofstay'] / 24, df['overleden_IC'],
            label=label)
        # ax = kmf.plot_cumulative_density(
        ax = kmf.plot_survival_function(
            ax=ax, at_risk_counts=arc, color=color, ci_alpha=0.15, ls=ls)
    return


###############################################################################
# Area plot for sepsis status over time (up to day 14)

for unit_type in ['icu', 'mcu']:
    if unit_type == 'icu':
        descriptors = descriptors_icu
    elif unit_type == 'mcu':
        descriptors = descriptors_mcu

    n_admissions = descriptors.loc[(
        descriptors['discard'] == 0), 'admissionid'].nunique()
    n_time = (MAX_TIME - MIN_TIME)
    time_range = list(range(n_time))

    # Expand out the descriptors dataframe to account for status on each
    # day (including discharge/death)
    sepsis_status = pd.DataFrame(
        index=range(n_admissions * n_time), columns=['admissionid', 'time'])
    sepsis_status['admissionid'] = np.repeat(
        descriptors.loc[(descriptors['discard'] == 0), 'admissionid'].unique(),
        n_time)
    sepsis_status['time'] = np.tile(time_range, n_admissions)
    antibiotics_add = drugitems_abx_all.loc[
           (drugitems_abx_all['time'] >= 0) &
           (drugitems_abx_all['time'] < n_time) &
           (drugitems_abx_all['intravenous'])]
    antibiotics_add = antibiotics_add[['admissionid', 'time']].assign(abx=True)
    antibiotics_add = antibiotics_add[~antibiotics_add.duplicated()]
    sepsis_status = pd.merge(
        sepsis_status, antibiotics_add,
        on=['admissionid', 'time'], how='left')
    sepsis_status['no_abx'] = ~(sepsis_status['abx'] == 1)
    sepsis_status['abx'] = (sepsis_status['abx'] == 1)
    sepsis_add = sepsis.drop(columns=['sepsis_episode', 'septic_shock'])
    sepsis_add.rename(columns={
        'continued_sepsis_episode': 'sepsis_episode',
        'continued_septic_shock': 'septic_shock'}, inplace=True)
    sepsis_add = sepsis_add[
        ['admissionid', 'time', 'sepsis_episode', 'septic_shock']]
    # A sepsis episode on day -1 should also cover day 0 here (due to the
    # consecutive day SOFA score change)!
    sepsis_ind = sepsis_add.loc[
            (sepsis_add['admissionid'].diff(-1) == 0) &
            (sepsis_add['time'].diff(-1) == -1) &
            (sepsis_add['sepsis_episode'])
        ].index
    sepsis_add.loc[sepsis_ind + 1, 'sepsis_episode'] = True
    sepsis_ind = sepsis_add.loc[
            (sepsis_add['admissionid'].diff(-1) == 0) &
            (sepsis_add['time'].diff(-1) == -1) &
            (sepsis_add['septic_shock'])
        ].index
    sepsis_add.loc[sepsis_ind + 1, 'septic_shock'] = True
    sepsis_add = sepsis_add.loc[
           (sepsis_add['time'] >= 0) &
           (sepsis_add['time'] < n_time)]
    sepsis_status = pd.merge(
        sepsis_status, sepsis_add,
        on=['admissionid', 'time'], how='left')
    sepsis_columns = ['sepsis_episode', 'septic_shock']
    sepsis_status[sepsis_columns] = sepsis_status[sepsis_columns].fillna(False)

    descriptors_add = descriptors[
        ['admissionid', 'overleden_IC', 'admittedat', 'dischargedat']]
    sepsis_status = pd.merge(
        sepsis_status, descriptors_add, on='admissionid', how='left')

    sepsis_status['discharge_time'] = (
        sepsis_status['dischargedat'] - sepsis_status['admittedat'])
    sepsis_status['discharge_time'] /= (1000*60*60*24)
    sepsis_status['discharged'] = (
        (~sepsis_status['overleden_IC']) &
        (sepsis_status['time'] > sepsis_status['discharge_time']))
    sepsis_status['died'] = (
        (sepsis_status['overleden_IC']) &
        (sepsis_status['time'] > sepsis_status['discharge_time']))

    for column in ['abx', 'no_abx', 'sepsis_episode', 'septic_shock']:
        sepsis_status.loc[sepsis_status['died'], column] = False
        sepsis_status.loc[sepsis_status['discharged'], column] = False
    for column in ['abx', 'no_abx', 'sepsis_episode']:
        sepsis_status.loc[sepsis_status['septic_shock'], column] = False
    for column in ['abx', 'no_abx']:
        sepsis_status.loc[sepsis_status['sepsis_episode'], column] = False

    sepsis_admission_columns = ['admissionid', 'sepsis', 'septic_shock']
    sepsis_admission_columns += ['antibiotics_admission']
    descriptors_add = descriptors[sepsis_admission_columns]
    descriptors_add.rename(columns={
        'antibiotics_admission': 'abx_admission',
        'sepsis': 'sepsis_admission',
        'septic_shock': 'septic_shock_admission'}, inplace=True)
    sepsis_status = pd.merge(
        sepsis_status, descriptors_add, on='admissionid', how='left')
    sepsis_status = sepsis_status.fillna(False)
    sepsis_status['no_abx_admission'] = (sepsis_status['abx_admission'] == 0)
    sepsis_status.loc[(
            sepsis_status['septic_shock_admission'] == 1),
        ['sepsis_admission', 'abx_admission', 'no_abx_admission']] = False
    sepsis_status.loc[(
            sepsis_status['sepsis_admission'] == 1),
        ['abx_admission', 'no_abx_admission']] = False

    all_admission_counts = counts_fun(sepsis_status)
    sepsis_counts = counts_fun(
        sepsis_status.loc[sepsis_status['sepsis_admission']])
    septic_shock_counts = counts_fun(
        sepsis_status.loc[sepsis_status['septic_shock_admission']])
    abx_counts = counts_fun(
        sepsis_status.loc[sepsis_status['abx_admission']])
    no_abx_counts = counts_fun(
        sepsis_status.loc[sepsis_status['no_abx_admission']])

    ###########################################################################
    # Trajectories plot

    figure_name = unit_type + '/' + 'trajectories_' + unit_type + '.pdf'
    with PdfPages(inputs.output_file_path + figure_name) as pdf:
        fig, ax = plt.subplots(nrows=3, ncols=2, sharex=True, figsize=[10, 10])
        all_admission_counts.plot.bar(
            x='Time since admission (d)', ax=ax[0, 0],
            legend=False, width=1, stacked=True, align='edge',
            color=ql.Dark2_6.mpl_colors, ylabel='Number of admissions')
        fig.legend(bbox_to_anchor=[0.9, 0.9])
        ax[0, 0].set_ylim([0, n_admissions])
        ax[0, 0].set_title('(a) All ICU admissions', fontsize=10)
        ax[0, 1].axis('off')
        # Septic shock
        septic_shock_counts.plot.bar(
            x='Time since admission (d)', ax=ax[1, 0], rot=False,
            legend=False, width=1, stacked=True, align='edge',
            color=ql.Dark2_6.mpl_colors, ylabel='Number of admissions')
        ax[1, 0].set_ylim([0, septic_shock_counts.loc[0].sum()])
        ax[1, 0].set_title('(b) Septic shock at admission', fontsize=10)
        # Sepsis without shock
        sepsis_counts.plot.bar(
            x='Time since admission (d)', ax=ax[1, 1],
            legend=False, width=1, stacked=True, align='edge',
            color=ql.Dark2_6.mpl_colors, ylabel='Number of admissions')
        ax[1, 1].set_ylim([0, sepsis_counts.loc[0].sum()])
        ax[1, 1].set_title(
            '(c) Sepsis without shock at admission', fontsize=10)
        # Antibiotics
        abx_counts.plot.bar(
            x='Time since admission (d)', ax=ax[2, 0], rot=False,
            legend=False, width=1, stacked=True, align='edge',
            color=ql.Dark2_6.mpl_colors, ylabel='Number of admissions')
        ax[2, 0].set_ylim([0, abx_counts.loc[0].sum()])
        ax[2, 0].set_title(
            '(d) Antibiotics without sepsis at admission', fontsize=10)
        # No antibiotics
        no_abx_counts.plot.bar(
            x='Time since admission (d)', ax=ax[2, 1], rot=False,
            legend=False, width=1, stacked=True, align='edge',
            color=ql.Dark2_6.mpl_colors, ylabel='Number of admissions')
        ax[2, 1].set_ylim([0, no_abx_counts.loc[0].sum()])
        ax[2, 1].set_title('(e) Not on antibiotics at admission', fontsize=10)
        ymin, ymax = ax[0, 0].get_ylim()
        for x in range(n_time):
            ax[0, 0].plot(
                [x, x], [ymin, ymax],
                ls='--', color='black', lw=0.5)
        ax[0, 0].set_xlim([0, n_time])
        for row in range(1, 3):
            for col in range(2):
                ymin, ymax = ax[row, col].get_ylim()
                for x in range(n_time):
                    ax[row, col].plot(
                        [x, x], [ymin, ymax],
                        ls='--', color='black', lw=0.5)
                ax[row, col].set_xlim([0, n_time])
        plt.tight_layout()
        pdf.savefig(fig)
        plt.close()

###############################################################################
# Kaplan-Meier curves

for unit_type in ['icu', 'mcu']:
    if unit_type == 'icu':
        descriptors = descriptors_icu
    elif unit_type == 'mcu':
        descriptors = descriptors_mcu

    kmf = lifelines.KaplanMeierFitter()
    figure_name = unit_type + '/' + 'sepsis_km_' + unit_type + '.pdf'
    with PdfPages(inputs.output_file_path + figure_name) as pdf:
        fig, ax = plt.subplots(nrows=3, ncols=2, sharex=True, figsize=[10, 10])
        kmf_plot(
            descriptors.loc[descriptors['septic_shock']],
            ax=ax[0, 0], kmf=kmf, ls='solid',
            label='Septic shock', color=RdYlBu_4_hex[0])
        kmf_plot(
            descriptors.loc[descriptors['sepsis_wo_shock']],
            ax=ax[0, 0], kmf=kmf, ls='dashed',
            label='Sepsis without shock', color=RdYlBu_4_hex[1])
        kmf_plot(
            descriptors.loc[descriptors['antibiotics_wo_sepsis']],
            ax=ax[0, 0], kmf=kmf, ls=(0, (3, 1, 1, 1)),
            label='On antibiotics (no sepsis)', color=RdYlBu_4_hex[2])
        kmf_plot(
            descriptors.loc[descriptors['no_antibiotics']],
            ax=ax[0, 0], kmf=kmf, ls=(0, (3, 2, 1, 2, 1, 2)),
            label='Not on antibiotics', color=RdYlBu_4_hex[3])
        title_str = '(a) Survival curves (all ICU admissions) \n'
        title_str += 'by admission sepsis status'
        ax[0, 0].set_title(title_str, fontsize=10)
        ax[0, 0].set_ylabel('Survival probability')
        ax[0, 0].set_xlabel('Time since admission (d)')
        ax[0, 0].set_ylim([0, 1])
        ax[0, 0].set_xlim([0, n_time])
        ax[0, 0].legend(loc='lower left')
        # By SOFA
        kmf_plot(
            descriptors.loc[
                (descriptors['sofa_total_score'] >= 0) &
                (descriptors['sofa_total_score'] < 4)],
            ax=ax[0, 1], ls='solid', kmf=kmf,
            label='Max SOFA 0-3', color=Dark2_6_hex[0])
        kmf_plot(
            descriptors.loc[
                (descriptors['sofa_total_score'] >= 4) &
                (descriptors['sofa_total_score'] < 8)],
            ax=ax[0, 1], ls='dashed', kmf=kmf,
            label='Max SOFA 4-7', color=Dark2_6_hex[1])
        kmf_plot(
            descriptors.loc[
                (descriptors['sofa_total_score'] >= 8) &
                (descriptors['sofa_total_score'] < 12)],
            ax=ax[0, 1], ls=(0, (1, 1)), kmf=kmf,
            label='Max SOFA 8-11', color=Dark2_6_hex[2])
        kmf_plot(
            descriptors.loc[
                (descriptors['sofa_total_score'] >= 12) &
                (descriptors['sofa_total_score'] < 16)],
            ax=ax[0, 1], ls=(0, (3, 1, 1, 2)), kmf=kmf,
            label='Max SOFA 12-15', color=Dark2_6_hex[3])
        kmf_plot(
            descriptors.loc[
                (descriptors['sofa_total_score'] >= 16)],
            ax=ax[0, 1], ls=(0, (7, 1)), kmf=kmf,
            label='Max SOFA >15', color=Dark2_6_hex[4])
        ax[0, 1].set_title(
            '(b) Survival curves (all ICU admissions) ' +
            '\n by admission SOFA score',
            fontsize=10)
        ax[0, 1].set_xlabel('Time since admission (d)')
        ax[0, 1].set_ylim([0, 1])
        ax[0, 1].set_xlim([0, n_time])
        ax[0, 1].legend(loc='lower left')
        # Septic shock at admission
        kmf_plot(
            descriptors.loc[
                descriptors['septic_shock'] &
                (descriptors['sofa_total_score'] >= 0) &
                (descriptors['sofa_total_score'] < 4)],
            ax=ax[1, 0], ls='solid', kmf=kmf,
            label='Max SOFA 0-3', color=Dark2_6_hex[0])
        kmf_plot(
            descriptors.loc[
                descriptors['septic_shock'] &
                (descriptors['sofa_total_score'] >= 4) &
                (descriptors['sofa_total_score'] < 8)],
            ax=ax[1, 0], ls='dashed', kmf=kmf,
            label='Max SOFA 4-7', color=Dark2_6_hex[1])
        kmf_plot(
            descriptors.loc[
                descriptors['septic_shock'] &
                (descriptors['sofa_total_score'] >= 8) &
                (descriptors['sofa_total_score'] < 12)],
            ax=ax[1, 0], ls=(0, (1, 1)), kmf=kmf,
            label='Max SOFA 8-11', color=Dark2_6_hex[2])
        kmf_plot(
            descriptors.loc[
                descriptors['septic_shock'] &
                (descriptors['sofa_total_score'] >= 12) &
                (descriptors['sofa_total_score'] < 16)],
            ax=ax[1, 0], ls=(0, (3, 1, 1, 2)), kmf=kmf,
            label='Max SOFA 12-15', color=Dark2_6_hex[3])
        kmf_plot(
            descriptors.loc[
                descriptors['septic_shock'] &
                (descriptors['sofa_total_score'] >= 16)],
            ax=ax[1, 0], ls=(0, (7, 1)), kmf=kmf,
            label='Max SOFA >15', color=Dark2_6_hex[4])
        ax[1, 0].set_title(
            '(c) Survival curves for patients\nwith septic shock\n' +
            'by admission SOFA score',
            fontsize=10)
        ax[1, 0].set_xlabel('Time since admission (d)')
        ax[1, 0].set_ylim([0, 1])
        ax[1, 0].set_xlim([0, n_time])
        ax[1, 0].legend(loc='lower left')
        # if ax[1, 0].get_legend() is not None:
        #     ax[1, 0].get_legend().remove()
        # Sepsis without shock at admission
        kmf_plot(
            descriptors.loc[
                descriptors['sepsis_wo_shock'] &
                (descriptors['sofa_total_score'] >= 0) &
                (descriptors['sofa_total_score'] < 4)],
            ax=ax[1, 1], ls='solid', kmf=kmf,
            label='Max SOFA 0-3', color=Dark2_6_hex[0])
        kmf_plot(
            descriptors.loc[
                descriptors['sepsis_wo_shock'] &
                (descriptors['sofa_total_score'] >= 4) &
                (descriptors['sofa_total_score'] < 8)],
            ax=ax[1, 1], ls='dashed', kmf=kmf,
            label='Max SOFA 4-7', color=Dark2_6_hex[1])
        kmf_plot(
            descriptors.loc[
                descriptors['sepsis_wo_shock'] &
                (descriptors['sofa_total_score'] >= 8) &
                (descriptors['sofa_total_score'] < 12)],
            ax=ax[1, 1], ls=(0, (1, 1)), kmf=kmf,
            label='Max SOFA 8-11', color=Dark2_6_hex[2])
        kmf_plot(
            descriptors.loc[
                descriptors['sepsis_wo_shock'] &
                (descriptors['sofa_total_score'] >= 12) &
                (descriptors['sofa_total_score'] < 16)],
            ax=ax[1, 1], ls=(0, (3, 1, 1, 2)), kmf=kmf,
            label='Max SOFA 12-15', color=Dark2_6_hex[3])
        kmf_plot(
            descriptors.loc[
                descriptors['sepsis_wo_shock'] &
                (descriptors['sofa_total_score'] >= 16)],
            ax=ax[1, 1], ls=(0, (7, 1)), kmf=kmf,
            label='Max SOFA >15', color=Dark2_6_hex[4])
        ax[1, 1].set_title(
            '(d) Survival curves for patients\nwith sepsis without shock\n' +
            'by admission SOFA score',
            fontsize=10)
        ax[1, 1].set_ylabel('Survival probability')
        ax[1, 1].set_xlabel('Time since admission (d)')
        ax[1, 1].set_ylim([0, 1])
        ax[1, 1].set_xlim([0, n_time])
        ax[1, 1].legend(loc='lower left')
        # if ax[1, 1].get_legend() is not None:
        #     ax[1, 1].get_legend().remove()
        # Antibiotics without sepsis
        kmf_plot(
            descriptors.loc[
                descriptors['antibiotics_wo_sepsis'] &
                (descriptors['sofa_total_score'] >= 0) &
                (descriptors['sofa_total_score'] < 4)],
            ax=ax[2, 0], ls='solid', kmf=kmf,
            label='Max SOFA 0-3', color=Dark2_6_hex[0])
        kmf_plot(
            descriptors.loc[
                descriptors['antibiotics_wo_sepsis'] &
                (descriptors['sofa_total_score'] >= 4) &
                (descriptors['sofa_total_score'] < 8)],
            ax=ax[2, 0], ls='dashed', kmf=kmf,
            label='Max SOFA 4-7', color=Dark2_6_hex[1])
        kmf_plot(
            descriptors.loc[
                descriptors['antibiotics_wo_sepsis'] &
                (descriptors['sofa_total_score'] >= 8) &
                (descriptors['sofa_total_score'] < 12)],
            ax=ax[2, 0], ls=(0, (1, 1)), kmf=kmf,
            label='Max SOFA 8-11', color=Dark2_6_hex[2])
        kmf_plot(
            descriptors.loc[
                descriptors['antibiotics_wo_sepsis'] &
                (descriptors['sofa_total_score'] >= 12) &
                (descriptors['sofa_total_score'] < 16)],
            ax=ax[2, 0], ls=(0, (3, 1, 1, 2)), kmf=kmf,
            label='Max SOFA 12-15', color=Dark2_6_hex[3])
        kmf_plot(
            descriptors.loc[
                descriptors['antibiotics_wo_sepsis'] &
                (descriptors['sofa_total_score'] >= 16)],
            ax=ax[2, 0], ls=(0, (7, 1)), kmf=kmf,
            label='Max SOFA >15', color=Dark2_6_hex[4])
        ax[2, 0].set_title(
            '(e) Survival curves for patients\nreceiving antibiotics without' +
            'sepsis\nby admission SOFA score',
            fontsize=10)
        ax[2, 0].set_ylabel('Survival probability')
        ax[2, 0].set_xlabel('Time since admission (d)')
        ax[2, 0].set_ylim([0, 1])
        ax[2, 0].set_xlim([0, n_time])
        ax[2, 0].legend(loc='lower left')
        # if ax[2, 0].get_legend() is not None:
        #     ax[2, 0].get_legend().remove()
        # No antibiotics
        kmf_plot(
            descriptors.loc[
                descriptors['no_antibiotics'] &
                (descriptors['sofa_total_score'] >= 0) &
                (descriptors['sofa_total_score'] < 4)],
            ax=ax[2, 1], ls='solid', kmf=kmf,
            label='Max SOFA 0-3', color=Dark2_6_hex[0])
        kmf_plot(
            descriptors.loc[
                descriptors['no_antibiotics'] &
                (descriptors['sofa_total_score'] >= 4) &
                (descriptors['sofa_total_score'] < 8)],
            ax=ax[2, 1], ls='dashed', kmf=kmf,
            label='Max SOFA 4-7', color=Dark2_6_hex[1])
        kmf_plot(
            descriptors.loc[
                descriptors['no_antibiotics'] &
                (descriptors['sofa_total_score'] >= 8) &
                (descriptors['sofa_total_score'] < 12)],
            ax=ax[2, 1], ls=(0, (1, 1)), kmf=kmf,
            label='Max SOFA 8-11', color=Dark2_6_hex[2])
        kmf_plot(
            descriptors.loc[
                descriptors['no_antibiotics'] &
                (descriptors['sofa_total_score'] >= 12) &
                (descriptors['sofa_total_score'] < 16)],
            ax=ax[2, 1], ls=(0, (3, 1, 1, 2)), kmf=kmf,
            label='Max SOFA 12-15', color=Dark2_6_hex[3])
        kmf_plot(
            descriptors.loc[
                descriptors['no_antibiotics'] &
                (descriptors['sofa_total_score'] >= 16)],
            ax=ax[2, 1], ls=(0, (7, 1)), kmf=kmf,
            label='Max SOFA >15', color=Dark2_6_hex[4])
        ax[2, 1].set_title(
            '(f) Survival curves for patients\nnot receiving antibiotics\n' +
            'by admission SOFA score',
            fontsize=10)
        ax[2, 1].set_xlabel('Time since admission (d)')
        ax[2, 1].set_ylim([0, 1])
        ax[2, 1].set_xlim([0, n_time])
        ax[2, 1].legend(loc='lower left')
        # if ax[2, 1].get_legend() is not None:
        #     ax[2, 1].get_legend().remove()
        for row in range(3):
            for col in range(2):
                ymin, ymax = ax[row, col].get_ylim()
                for x in range(n_time):
                    ax[row, col].plot(
                        [x, x], [ymin, ymax],
                        ls='--', color='black', lw=0.5)
                ax[row, col].set_xlim([0, n_time])
        plt.tight_layout()
        pdf.savefig(fig)
        plt.close()

###############################################################################
# SOFA histograms

sofa_hist = pd.merge(
    sofa, sepsis[['admissionid', 'time', 'sepsis_episode']],
    on=['admissionid', 'time'], how='left')
descriptors_add = ['admissionid', 'sepsis', 'discard', 'location_IC']
descriptors_add += ['overleden_IC']
sofa_hist = pd.merge(
    sofa_hist, descriptors_all[descriptors_add], on='admissionid', how='left')
sofa_hist = sofa_hist.loc[(sofa_hist['discard'] == 0)]
sofa_hist.loc[(
        (sofa_hist['time'] == 0) & (sofa_hist['sepsis'] == 1)),
    'sepsis_episode'] = True
sofa_hist.drop(columns=['discard', 'sepsis'], inplace=True)
sofa_hist.loc[sofa_hist['sepsis_episode'], 'day_post_sepsis'] = 0
sofa_hist.loc[(
        (sofa_hist['overleden_IC'] == 1) &
        sofa_hist['admissionid'].isin(
            sofa_hist.loc[sofa_hist['sepsis_episode'], 'admissionid']) &
        (sofa_hist['admissionid'].duplicated(keep='last') == 0)),
    'day_pre_death'] = 0
sofa_hist.loc[(
        (sofa_hist['overleden_IC'] == 0) &
        sofa_hist['admissionid'].isin(
            sofa_hist.loc[sofa_hist['sepsis_episode'], 'admissionid']) &
        (sofa_hist['admissionid'].duplicated(keep='last') == 0)),
    'day_pre_discharge'] = 0

sofa_hist = sofa_hist.reset_index(drop=True)

jj = sofa_hist.loc[(sofa_hist['day_post_sepsis'] == 0)].index
n = sofa_hist.shape[0]
for ii in range(1, 7):
    kk = (
        sofa_hist.loc[(jj + ii)[(jj + ii) < n], 'admissionid'].values ==
        sofa_hist.loc[jj[(jj + ii) < n], 'admissionid'].values)
    sofa_hist.loc[jj[(jj + ii) < n][kk] + ii, 'day_post_sepsis'] = (
        sofa_hist.loc[jj[(jj + ii) < n][kk] + ii, 'time'].values -
        sofa_hist.loc[jj[(jj + ii) < n][kk], 'time'].values)
sofa_hist.loc[sofa_hist['day_post_sepsis'] > 6, 'day_post_sepsis'] = np.nan

jj = sofa_hist.loc[(sofa_hist['day_pre_death'] == 0)].index
for ii in range(1, 7):
    kk = (
        sofa_hist.loc[(jj - ii)[(jj - ii) >= 0], 'admissionid'].values ==
        sofa_hist.loc[jj[(jj - ii) >= 0], 'admissionid'].values)
    sofa_hist.loc[jj[(jj - ii) >= 0][kk] - ii, 'day_pre_death'] = (
        sofa_hist.loc[jj[(jj - ii) >= 0][kk] - ii, 'time'].values -
        sofa_hist.loc[jj[(jj - ii) >= 0][kk], 'time'].values)
sofa_hist.loc[sofa_hist['day_pre_death'] < -6, 'day_pre_death'] = np.nan

jj = sofa_hist.loc[(sofa_hist['day_pre_discharge'] == 0)].index
for ii in range(1, 7):
    kk = (
        sofa_hist.loc[(jj - ii)[(jj - ii) >= 0], 'admissionid'].values ==
        sofa_hist.loc[jj[(jj - ii) >= 0], 'admissionid'].values)
    sofa_hist.loc[jj[(jj - ii) >= 0][kk] - ii, 'day_pre_discharge'] = (
        sofa_hist.loc[jj[(jj - ii) >= 0][kk] - ii, 'time'].values -
        sofa_hist.loc[jj[(jj - ii) >= 0][kk], 'time'].values)
sofa_hist.loc[
    sofa_hist['day_pre_discharge'] < -6, 'day_pre_discharge'] = np.nan

subscore_labels = ['Cardiovascular', 'Respiration', 'Coagulation', 'Renal']
subscore_labels += ['Liver', 'CNS']
sofa_subscores = ['sofa_' + x.lower() + '_score' for x in subscore_labels]
color_list = [sq.Blues_7.mpl_colors[1:], sq.Greens_7.mpl_colors[1:]]
color_list += [sq.Oranges_7.mpl_colors[1:], sq.Purples_7.mpl_colors[1:]]
color_list += [sq.Reds_7.mpl_colors[1:], sq.PuRd_7.mpl_colors[1:]]

sofa_hist_icu = sofa_hist.loc[(sofa_hist['location_IC'] == 1)]
sofa_hist_mcu = sofa_hist.loc[(sofa_hist['location_IC'] == 0)]

sofa_hist_df = pd.DataFrame(index=range(6))
sofa_hist_df['x'] = [str(x) for x in range(5)] + ['NA']

for unit_type in ['icu', 'mcu']:
    if unit_type == 'icu':
        sofa_hist = sofa_hist_icu
    elif unit_type == 'mcu':
        sofa_hist = sofa_hist_mcu

    figure_name = unit_type + '/' + 'sofa_pre_death_' + unit_type + '.pdf'
    with PdfPages(inputs.output_file_path + figure_name) as pdf:
        fig, ax = plt.subplots(
            nrows=6, ncols=6, sharex=True, sharey=True, figsize=[10, 10])
        for ii in range(6):
            for jj in range(6):
                temp_ind = (sofa_hist['day_pre_death'] == -ii)
                temp_vals = sofa_hist.loc[
                    temp_ind, sofa_subscores[jj]].fillna(5).value_counts()
                temp_vals = (temp_vals / temp_ind.sum()*100).reset_index()
                temp_vals = temp_vals.sort_values(by='index')
                temp_vals = temp_vals.reset_index(drop=True)
                temp_vals['x'] = temp_vals['index'].astype(int).astype(str)
                if (temp_vals['index'] == 5).any():
                    temp_vals.loc[5, 'x'] = 'NA'
                sofa_hist_df[
                    'day_pre_death=' + str(-ii) + ', ' +
                    sofa_subscores[jj] + ', ' + unit_type] = (
                        temp_vals[sofa_subscores[jj]])
                ax[ii, jj].bar(
                    temp_vals['index'], temp_vals[sofa_subscores[jj]],
                    width=1, color=color_list[jj])
                ax[ii, jj].set_ylim([0, 79.9])
                yticks = [0, 20, 40, 60]
                ytick_labels = [str(x) + '%' for x in yticks]
                ax[ii, jj].set_yticks(yticks)
                ax[ii, jj].set_yticklabels(ytick_labels)
                ax[ii, jj].set_xticks(temp_vals['index'])
                ax[ii, jj].set_xticklabels(temp_vals['x'])
        for ii in range(6):
            ax[ii, 0].set_ylabel('Day -' + str(ii))
        for jj in range(6):
            ax[5, jj].set_xlabel(subscore_labels[jj])
        fig.supylabel(
            'SOFA scores for patients with sepsis prior to ' +
            unit_type.upper() + ' death (day 0 = day of death)')
        fig.supxlabel('Score for each SOFA component')
        plt.tight_layout(w_pad=0, h_pad=0)
        plt.subplots_adjust(wspace=0, hspace=0)
        pdf.savefig(fig)
        plt.close()

    figure_name = unit_type + '/' + 'sofa_pre_discharge_' + unit_type + '.pdf'
    with PdfPages(inputs.output_file_path + figure_name) as pdf:
        fig, ax = plt.subplots(
            nrows=6, ncols=6, sharex=True, sharey=True, figsize=[10, 10])
        for ii in range(6):
            for jj in range(6):
                temp_ind = (sofa_hist['day_pre_discharge'] == -ii)
                temp_vals = sofa_hist.loc[
                    temp_ind, sofa_subscores[jj]].fillna(5).value_counts()
                temp_vals = (temp_vals / temp_ind.sum()*100).reset_index()
                temp_vals = temp_vals.sort_values(by='index')
                temp_vals = temp_vals.reset_index(drop=True)
                temp_vals['x'] = temp_vals['index'].astype(int).astype(str)
                if (temp_vals['index'] == 5).any():
                    temp_vals.loc[5, 'x'] = 'NA'
                sofa_hist_df[
                    'day_pre_discharge=' + str(-ii) + ', ' +
                    sofa_subscores[jj] + ', ' + unit_type] = (
                        temp_vals[sofa_subscores[jj]])
                ax[ii, jj].bar(
                    temp_vals['index'], temp_vals[sofa_subscores[jj]],
                    width=1, color=color_list[jj])
                ax[ii, jj].set_ylim([0, 79.9])
                yticks = [0, 20, 40, 60]
                ytick_labels = [str(x) + '%' for x in yticks]
                ax[ii, jj].set_yticks(yticks)
                ax[ii, jj].set_yticklabels(ytick_labels)
                ax[ii, jj].set_xticks(temp_vals['index'])
                ax[ii, jj].set_xticklabels(temp_vals['x'])
        for ii in range(6):
            ax[ii, 0].set_ylabel('Day -' + str(ii))
        for jj in range(6):
            ax[5, jj].set_xlabel(subscore_labels[jj])
        fig.supylabel(
            'SOFA scores for patients with sepsis prior to ' +
            unit_type.upper() +
            ' discharge (day 0 = day of discharge)')
        fig.supxlabel('Score for each SOFA component')
        plt.tight_layout(w_pad=0, h_pad=0)
        plt.subplots_adjust(wspace=0, hspace=0)
        pdf.savefig(fig)
        plt.close()

    figure_name = unit_type + '/' + 'sofa_post_sepsis_' + unit_type + '.pdf'
    with PdfPages(inputs.output_file_path + figure_name) as pdf:
        fig, ax = plt.subplots(
            nrows=6, ncols=6, sharex=True, sharey=True, figsize=[10, 10])
        for ii in range(6):
            for jj in range(6):
                temp_ind = (sofa_hist['day_post_sepsis'] == (5 - ii))
                temp_vals = sofa_hist.loc[
                    temp_ind, sofa_subscores[jj]].fillna(5).value_counts()
                temp_vals = (temp_vals / temp_ind.sum()*100).reset_index()
                temp_vals = temp_vals.sort_values(by='index')
                temp_vals = temp_vals.reset_index(drop=True)
                temp_vals['x'] = temp_vals['index'].astype(int).astype(str)
                if (temp_vals['index'] == 5).any():
                    temp_vals.loc[5, 'x'] = 'NA'
                sofa_hist_df[
                    'day_post_sepsis=' + str(-ii) + ', ' +
                    sofa_subscores[jj] + ', ' + unit_type] = (
                        temp_vals[sofa_subscores[jj]])
                ax[ii, jj].bar(
                    temp_vals['index'], temp_vals[sofa_subscores[jj]],
                    width=1, color=color_list[jj])
                ax[ii, jj].set_ylim([0, 79.9])
                yticks = [0, 20, 40, 60]
                ytick_labels = [str(x) + '%' for x in yticks]
                ax[ii, jj].set_yticks(yticks)
                ax[ii, jj].set_yticklabels(ytick_labels)
                ax[ii, jj].set_xticks(temp_vals['index'])
                ax[ii, jj].set_xticklabels(temp_vals['x'])
        for ii in range(6):
            ax[ii, 0].set_ylabel('Day ' + str(5 - ii))
        for jj in range(6):
            ax[5, jj].set_xlabel(subscore_labels[jj])
        fig.supylabel(
            'SOFA scores by day post ' + unit_type.upper() +
            ' sepsis episode')
        fig.supxlabel('Score for each SOFA component')
        plt.tight_layout(w_pad=0, h_pad=0)
        plt.subplots_adjust(wspace=0, hspace=0)
        pdf.savefig(fig)
        plt.close()

sofa_hist_df = sofa_hist_df.T
sofa_hist_df.iloc[1:] = sofa_hist_df.iloc[1:].fillna(0).apply(
    lambda x: x.round(2))
sofa_hist_df.to_csv(
    inputs.additional_file_path + 'tables/sofa_hist_values.csv')

###############################################################################
# Antibiotic usage (all) per day

for include_prophylactic in [True, False]:
    for unit_type in ['icu', 'mcu']:
        if (unit_type == 'icu') & (include_prophylactic is True):
            abx_usage = abx_usage_icu
        elif (unit_type == 'icu') & (include_prophylactic is False):
            abx_usage = abx_np_usage_icu
        elif (unit_type == 'mcu') & (include_prophylactic is True):
            abx_usage = abx_usage_mcu
        elif (unit_type == 'mcu') & (include_prophylactic is False):
            abx_usage = abx_np_usage_mcu

        abx_daily_usage = abx_usage.groupby(['item_en', 'time']).agg(
                n_item_time=pd.NamedAgg(
                    column='admissionid', aggfunc=lambda x: x.nunique())
            ).reset_index()

        n_patients_time = (
            (descriptors_all['lengthofstay'] // 24).apply(
                lambda x: min((x, MAX_TIME)))).value_counts().reset_index()
        n_patients_time = n_patients_time.sort_values(
            by='index', ascending=False)
        n_patients_time['n_patients'] = (
            n_patients_time['lengthofstay'].cumsum())
        n_patients_time.rename(columns={'index': 'time'}, inplace=True)

        abx_daily_usage = pd.merge(
            abx_daily_usage, n_patients_time, on='time', how='left')

        abx_daily_usage['p_item_time'] = (
            abx_daily_usage['n_item_time'] / abx_daily_usage['n_patients'])

        # We want ICU antibiotics to be ordered in terms of average percentage,
        # and for MCU to match the ICU order
        if (unit_type == 'icu'):
            p_item = abx_daily_usage.groupby('item_en').agg(
                    p_item=pd.NamedAgg(column='p_item_time', aggfunc='mean')
                ).reset_index()
            #
            # item_include_ind = (p_item['p_item'] > 0.02)
            # item_include = p_item.loc[item_include_ind, 'item_en']
            # item_exclude = p_item.loc[(item_include_ind == 0), 'item_en']

            n_courses = abx_usage.groupby('item_en').agg(
                    n_courses=pd.NamedAgg(column='new_course', aggfunc='sum')
                ).reset_index()

            p_item = pd.merge(p_item, n_courses, on='item_en', how='left')

            item_include = p_item.sort_values(
                'n_courses', ascending=False)['item_en'].head(n=10).values
            n_ex = p_item.shape[0] - 10
            item_exclude = p_item.sort_values(
                'n_courses', ascending=False)['item_en'].tail(n=n_ex).values

        abx_daily_usage_other = abx_daily_usage.loc[
            abx_daily_usage['item_en'].isin(item_exclude)]
        abx_daily_usage_other = abx_daily_usage_other.groupby('time').agg(
                n_item_time=pd.NamedAgg(column='n_item_time', aggfunc='sum'),
                p_item_time=pd.NamedAgg(column='p_item_time', aggfunc='sum')
            ).reset_index()
        abx_daily_usage_other['item_en'] = 'Other'

        if (unit_type == 'icu') & (include_prophylactic is True):
            p_item.loc[p_item.shape[0]] = [
                'Other', abx_daily_usage_other['p_item_time'].mean(), np.nan]

        abx_daily_usage = abx_daily_usage.loc[
            abx_daily_usage['item_en'].isin(item_include)]
        abx_daily_usage = pd.merge(
            abx_daily_usage, p_item, on='item_en', how='left')

        abx_daily_usage = abx_daily_usage.sort_values(
            by=['p_item', 'time'], ascending=True)
        abx_daily_usage = pd.concat([abx_daily_usage_other, abx_daily_usage])
        abx_daily_usage.loc[(
                (abx_daily_usage['item_en'] == 'Other')),
            'p_item'] = abx_daily_usage_other['p_item_time'].sum()

        max_p_item_time = abx_daily_usage['p_item_time'].max()

        y_dict = dict(zip(
            abx_daily_usage['item_en'].unique(),
            range(abx_daily_usage['item_en'].nunique())))
        abx_daily_usage['y'] = abx_daily_usage['item_en'].replace(y_dict)
        x_dict = dict(zip(
            abx_daily_usage['time'].unique(),
            range(abx_daily_usage['time'].nunique())))
        abx_daily_usage['x'] = abx_daily_usage['time'].replace(x_dict)
        ytick_labels = list(abx_daily_usage['item_en'].unique())

        xs = abx_daily_usage['x']
        ys = (
            abx_daily_usage['y'] + 0.5 -
            0.5 * abx_daily_usage['p_item_time'] / max_p_item_time)
        widths = [1 for x in abx_daily_usage['x']]
        heights = abx_daily_usage['p_item_time'] / max_p_item_time
        colors = (
            dv.RdYlBu_11.mpl_colors + dv.PRGn_11.mpl_colors +
            dv.BrBG_11.mpl_colors)
        colors = abx_daily_usage['y'].apply(lambda x: colors[x])

        if include_prophylactic is True:
            figure_name = 'abx_relative_usage_'
        else:
            figure_name = 'abx_relative_usage_non_prophylactic_'
        figure_name = unit_type + '/' + figure_name + unit_type + '_v1.pdf'
        with PdfPages(inputs.output_file_path + figure_name) as pdf:
            fig, ax = plt.subplots(figsize=[8, 10])
            # ax.set_aspect('equal')
            ax.set_xlim([min(xs) - 0.1, max(xs) + 1 + 0.1])
            ax.set_ylim([min(ys) - 0.1, max(ys) + 1 + 0.1])

            rect_params = dict(edgecolor='black', lw=1)
            rectangles = [
                Rectangle((x, y), w, h, facecolor=c, **rect_params)
                for x, y, w, h, c in zip(xs, ys, widths, heights, colors)]

            for rectangle in rectangles:
                ax.add_patch(rectangle)

            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.spines['bottom'].set_visible(False)
            ax.spines['left'].set_visible(False)
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_yticks([y + 0.5 for y in y_dict.values()], ytick_labels)
            ax.set_xticks([x + 0.5 for x in x_dict.values()], x_dict.keys())
            ax.set_title(
                'Relative use of antibiotics (%), ' + unit_type.upper())
            ax.set_xlabel('Time (days)')
            plt.tight_layout(h_pad=0, w_pad=0)
            plt.subplots_adjust(wspace=0, hspace=0)
            pdf.savefig(fig)
            plt.close()

###############################################################################
# Sepsis flow chart

n_sepsis_admissions = sepsis.groupby('admissionid').agg(
        n_episodes=pd.NamedAgg(column='new_sepsis_episode', aggfunc='sum')
    ).reset_index()

descriptors_columns = ['admissionid', 'location_IC', 'elective_surgical']
descriptors_columns += ['emergency_surgical', 'emergency_medical']
descriptors_columns += ['sepsis', 'overleden_IC', 'lengthofstay', 'discard']

sepsis_flowchart = pd.merge(
    descriptors_all[descriptors_columns], n_sepsis_admissions,
    on='admissionid', how='left')

sepsis_flowchart.rename(
    columns={'sepsis': 'sepsis_at_admission'}, inplace=True)

sepsis_flowchart['type'] = 'Emergency medical'
sepsis_flowchart.loc[(
    sepsis_flowchart['emergency_surgical']), 'type'] = 'Emergency surgery'
sepsis_flowchart.loc[(
    sepsis_flowchart['elective_surgical']), 'type'] = 'Elective surgery'
sepsis_flowchart['location'] = 'MCU'
sepsis_flowchart.loc[sepsis_flowchart['location_IC'], 'location'] = 'ICU'
sepsis_flowchart['destination'] = 'Discharged'
sepsis_flowchart.loc[(
    sepsis_flowchart['overleden_IC']), 'destination'] = 'ICU/MCU death'
sepsis_flowchart.loc[(
        sepsis_flowchart['lengthofstay'] > 15*24),
    'destination'] = 'Still in ICU/MCU'

sepsis_flowchart = sepsis_flowchart.loc[(sepsis_flowchart['discard'] == 0)]
sepsis_flowchart['n_episodes'] = sepsis_flowchart['n_episodes'].astype(int)
sepsis_flowchart.loc[(sepsis_flowchart['n_episodes'] > 2), 'n_episodes'] = 2
sepsis_flowchart['n_episodes'] = sepsis_flowchart['n_episodes'].astype(str)
sepsis_flowchart.loc[(
    sepsis_flowchart['n_episodes'] == '2'), 'n_episodes'] = '>1'

sepsis_flowchart_columns = ['type', 'location', 'sepsis_at_admission']
sepsis_flowchart_columns += ['n_episodes', 'destination']

sepsis_flowchart = sepsis_flowchart[sepsis_flowchart_columns].groupby(
    sepsis_flowchart_columns).size().reset_index()

sepsis_flowchart.rename(columns={0: 'n'}, inplace=True)

sepsis_flowchart['destination'] = sepsis_flowchart['destination'].replace(
        {'ICU/MCU death': 'ICU death', 'Still in ICU/MCU': 'Still in ICU'})
sepsis_flowchart['sepsis_at_admission'].replace(
    {True: 'Yes', False: 'No'}, inplace=True)

colors = [
    sq.Greens_6.mpl_colors, sq.Oranges_6.mpl_colors, sq.Purples_6.mpl_colors]

sepsis_flowchart_icu = sepsis_flowchart.loc[
    (sepsis_flowchart['location'] == 'ICU')]
sepsis_flowchart_mcu = sepsis_flowchart.loc[
    (sepsis_flowchart['location'] == 'MCU')]


def value_to_n(df):
    df_new = df.copy()
    for column in df.columns:
        df_new[column] = df[column].replace(
            dict(zip(df[column].unique(), range(df[column].nunique()))))
    return df_new


def get_grouped_column(df, column, groupby_columns, fun, as_array=False):
    columns = [column] + groupby_columns
    output = df[columns].groupby(groupby_columns).agg(
        new_col=pd.NamedAgg(column=column, aggfunc=fun))
    output.rename(columns={'new_col': column}, inplace=True)
    if as_array is True:
        output = output.values.reshape(-1)
    else:
        output = output.reset_index()
    return output


def get_df_c(df, columns, groupby_columns):
    columns_all = columns + groupby_columns
    df_cxs_all = get_grouped_column(df, 'cxs', columns_all, 'min')
    df_cwidths_all = get_grouped_column(df, 'cwidths', columns_all, 'sum')
    df_cxs = get_grouped_column(df, 'cxs', columns, 'min')
    df_cwidths = get_grouped_column(df, 'cwidths', columns, 'sum')
    df_new = pd.merge(
        df_cxs_all, df_cwidths_all, on=columns_all, how='left')
    df_new = pd.merge(df_new, df_cxs, on=columns, how='left')
    df_new = pd.merge(df_new, df_cwidths, on=columns, how='left')
    df_new['cxs'] = ((df_new['cxs_x'] - df_new['cxs_y']) / df_new['cwidths_y'])
    df_new['cwidths'] = (df_new['cwidths_x'] / df_new['cwidths_y'])
    df_new = df_new[columns_all + ['cxs', 'cwidths']]
    return df_new


for unit_type in ['icu', 'mcu']:
    if unit_type == 'icu':
        sepsis_flowchart = sepsis_flowchart_icu
    elif unit_type == 'mcu':
        sepsis_flowchart = sepsis_flowchart_mcu

    df = value_to_n(sepsis_flowchart.drop(columns=['n']))
    df = pd.merge(
        df.reset_index(), sepsis_flowchart.reset_index(),
        on='index', how='left', suffixes=['', '_text'])
    df.drop(columns=['index'], inplace=True)
    df.rename(columns={'n_text': 'n'}, inplace=True)
    n_sum = df['n'].sum()

    df['cwidths'] = (df['n'] / n_sum)
    df['cxs'] = [0] + list(df['cwidths'].values.cumsum()[:-1])

    color_columns = ['type', 'sepsis_at_admission', 'n_episodes']
    df['color'] = df[color_columns].apply(
        lambda x: colors[x[0]][3*x[1] + x[2]], axis=1)

    df_te_sort = value_to_n(sepsis_flowchart.drop(columns=['n']))
    df_te_sort.sort_values(by=['type', 'n_episodes'], inplace=True)
    df_te_sort = pd.merge(
        df_te_sort.reset_index(), sepsis_flowchart.reset_index(),
        on='index', how='left', suffixes=['', '_text'])
    df_te_sort.drop(columns=['index'], inplace=True)
    df_te_sort.rename(columns={'n_text': 'n'}, inplace=True)

    df_te_sort['cwidths'] = (df_te_sort['n'] / df_te_sort['n'].sum())
    df_te_sort['cxs'] = [0] + list(df_te_sort['cwidths'].values.cumsum()[:-1])

    color_columns = ['type', 'sepsis_at_admission', 'n_episodes']
    df_te_sort['color'] = df_te_sort[color_columns].apply(
        lambda x: colors[x[0]][3*x[1] + x[2]], axis=1)

    df_td_sort = value_to_n(sepsis_flowchart.drop(columns=['n']))
    df_td_sort.sort_values(by=['type', 'destination'], inplace=True)
    df_td_sort = pd.merge(
        df_td_sort.reset_index(), sepsis_flowchart.reset_index(),
        on='index', how='left', suffixes=['', '_text'])
    df_td_sort.drop(columns=['index'], inplace=True)
    df_td_sort.rename(columns={'n_text': 'n'}, inplace=True)

    df_td_sort['cwidths'] = (df_td_sort['n'] / df_td_sort['n'].sum())
    df_td_sort['cxs'] = [0] + list(df_td_sort['cwidths'].values.cumsum()[:-1])

    color_columns = ['type', 'sepsis_at_admission', 'n_episodes']
    df_td_sort['color'] = df_td_sort[color_columns].apply(
        lambda x: colors[x[0]][3*x[1] + x[2]], axis=1)

    df_t = df.loc[(df['type'].duplicated() == 0)]
    df_t['cwidths'] = get_grouped_column(
        df, 'cwidths', ['type'], 'sum', as_array=True)
    df_t['color'] = df_t['type'].apply(lambda x: colors[x][2])
    df_t['n'] = get_grouped_column(
        df, 'n', ['type'], 'sum', as_array=True)
    df_t.loc[df_t['n'] < 10, 'n'] = '<10'
    df_t['text'] = df_t['type_text'] + '\n(' + df_t['n'].astype(str) + ')'
    n_t = df_t.shape[0]

    df_d = df.loc[(df['destination'].duplicated() == 0)]
    df_d['cwidths'] = get_grouped_column(
        df, 'cwidths', ['destination'], 'sum', as_array=True)
    n_d = df_d.shape[0]

    df_ts = df.loc[(df[['type', 'sepsis_at_admission']].duplicated() == 0)]
    df_ts['cwidths'] = get_grouped_column(
        df, 'cwidths', ['type', 'sepsis_at_admission'], 'sum', as_array=True)
    df_ts['color'] = df_ts[['type', 'sepsis_at_admission']].apply(
        lambda x: colors[x[0]][3*x[1] + 1], axis=1)
    df_ts['n'] = get_grouped_column(
        df, 'n', ['type', 'sepsis_at_admission'], 'sum', as_array=True)
    df_ts.loc[df_ts['n'] < 10, 'n'] = '<10'
    df_ts['text'] = (
        df_ts['sepsis_at_admission_text'].astype(str) + '\n(' +
        df_ts['n'].astype(str) + ')')
    n_ts = df_ts.shape[0]

    df_te = df.loc[(df[['type', 'n_episodes']].duplicated() == 0)]
    df_te['cwidths'] = get_grouped_column(
        df, 'cwidths', ['type', 'n_episodes'], 'sum', as_array=True)
    df_te['color'] = df_te[['type', 'n_episodes']].apply(
        lambda x: colors[x[0]][3 + x[1]], axis=1)
    df_te['n'] = get_grouped_column(
        df, 'n', ['type', 'n_episodes'], 'sum', as_array=True)
    df_te.loc[df_te['n'] < 10, 'n'] = '<10'
    df_te['text'] = (
        df_te['n_episodes_text'] + '\n(' + df_te['n'].astype(str) + ')')
    n_te = df_te.shape[0]

    df_td = df.loc[(df[['type', 'destination']].duplicated() == 0)]
    df_td['cwidths'] = get_grouped_column(
        df, 'cwidths', ['type', 'destination'], 'sum', as_array=True)
    df_td['n'] = get_grouped_column(
        df, 'n', ['type', 'destination'], 'sum', as_array=True)
    df_td.loc[df_td['n'] < 10, 'n'] = '<10'
    df_td['text'] = (
        df_td['destination_text'] + '\n(' + df_td['n'].astype(str) + ')')
    n_td = df_td.shape[0]

    temp_columns = ['type', 'sepsis_at_admission', 'n_episodes']
    df_tse = df.loc[(df[temp_columns].duplicated() == 0)]
    df_tse['cwidths'] = get_grouped_column(
        df, 'cwidths', temp_columns, 'sum', as_array=True)
    n_tse = df_tse.shape[0]

    df_t_s = get_df_c(df, ['type'], ['sepsis_at_admission'])
    df_t_s['color'] = df_t_s[['type', 'sepsis_at_admission']].apply(
            lambda x: colors[x[0]][3*x[1] + 1], axis=1)
    n_t_s = df_t_s.shape[0]

    df_t_e = get_df_c(df_te_sort, ['type'], ['n_episodes'])
    n_t_e = df_t_e.shape[0]

    df_t_se = get_df_c(df, ['type'], ['sepsis_at_admission', 'n_episodes'])
    df_t_se['color'] = df_t_se[temp_columns].apply(
            lambda x: colors[x[0]][3*x[1] + x[2]], axis=1)
    n_t_se = df_t_se.shape[0]

    df_ts_e = get_df_c(df, ['type', 'sepsis_at_admission'], ['n_episodes'])
    df_ts_e['color'] = df_ts_e[temp_columns].apply(
        lambda x: colors[x[0]][3*x[1] + x[2]], axis=1)
    n_ts_e = df_ts_e.shape[0]

    df_te_s = get_df_c(
        df_te_sort, ['type', 'n_episodes'], ['sepsis_at_admission'])
    df_te_s['color'] = df_te_s[temp_columns].apply(
        lambda x: colors[x[0]][3*x[1] + x[2]], axis=1)
    n_te_s = df_te_s.shape[0]

    df_t_es = get_df_c(
        df_te_sort, ['type'], ['n_episodes', 'sepsis_at_admission'])
    df_t_es['color'] = df_t_es[temp_columns].apply(
        lambda x: colors[x[0]][3*x[1] + x[2]], axis=1)
    n_t_es = df_t_es.shape[0]

    df_td_es = get_df_c(
        df_td_sort,
        ['type', 'destination'], ['n_episodes', 'sepsis_at_admission'])
    df_td_es['color'] = df_td_es[temp_columns].apply(
        lambda x: colors[x[0]][3*x[1] + x[2]], axis=1)
    n_td_es = df_td_es.shape[0]

    df_te_sd = get_df_c(
        df_te_sort,
        ['type', 'n_episodes'], ['destination', 'sepsis_at_admission'])
    df_te_sd['color'] = df_te_sd[temp_columns].apply(
        lambda x: colors[x[0]][3*x[1] + x[2]], axis=1)
    n_te_sd = df_te_sd.shape[0]

    pad = 0.02
    n_rows = 5

    # padding to avoid line width obscuring colors
    xs = [4 - pad]
    ys = [2*(n_rows - 1) - pad]
    widths = [5 + 2*pad]
    heights = [1 + 2*pad]
    text_xs = [4 + 2.5]
    text_ys = [2*(n_rows - 1) + 0.5]
    texts = ['All admissions\n(n = ' + df['n'].sum().astype(str) + ')']

    cxs = list(4 + 5*df_tse['cxs'])
    cys = [2*(n_rows - 1) for x in range(n_tse)]
    cwidths = list(5*df_tse['cwidths'])
    cheights = [1 for x in range(n_tse)]
    ccolors = list(df_tse['color'])

    polys = list(zip(
        zip(4 + 5*df_tse['cxs'], [2*(n_rows - 1) for x in range(n_tse)]),
        zip(4 + 5*(df_tse['cxs'] + df_tse['cwidths']),
            [2*(n_rows - 1) for x in range(n_tse)]),
        zip(5*df_t_se['type'] + 3*(df_t_se['cxs'] + df_t_se['cwidths']),
            [2*(n_rows - 2) + 1 for x in range(n_t_se)]),
        zip(5*df_t_se['type'] + 3*df_t_se['cxs'],
            [2*(n_rows - 2) + 1 for x in range(n_t_se)])))
    pcolors = list(df_tse['color'])

    xs += list(5*df_t['type'] - pad)
    ys += [(2*(n_rows - 2) - pad) for x in range(n_t)]
    widths += [(3 + 2*pad) for x in range(n_t)]
    heights += [(1 + 2*pad) for x in range(n_t)]
    text_xs += list(5*df_t['type'] + 1.5)
    text_ys += [2*(n_rows - 2) + 0.5 for x in range(n_t)]
    texts += list(df_t['text'])

    cxs += list(5*df_t_se['type'] + 3*df_t_se['cxs'])
    cys += [2*(n_rows - 2) for x in range(n_t_se)]
    cwidths += list(3*df_t_se['cwidths'])
    cheights += [1 for x in range(n_t_se)]
    ccolors += list(df_t_se['color'])

    polys += list(zip(
        zip(5*df_t_se['type'] + 3*df_t_se['cxs'],
            [2*(n_rows - 2) for x in range(n_t_se)]),
        zip(5*df_t_se['type'] + 3*(df_t_se['cxs'] + df_t_se['cwidths']),
            [2*(n_rows - 2) for x in range(n_tse)]),
        zip(5*df_ts_e['type'] + df_ts_e['sepsis_at_admission'] +
            df_ts_e['cxs'] + df_ts_e['cwidths'] + 0.5 - 2*pad +
            4*pad*df_ts_e['sepsis_at_admission'],
            [2*(n_rows - 3) + 1 for x in range(n_ts_e)]),
        zip(5*df_ts_e['type'] + df_ts_e['sepsis_at_admission'] +
            df_ts_e['cxs'] + 0.5 - 2*pad +
            4*pad*df_ts_e['sepsis_at_admission'],
            [2*(n_rows - 3) + 1 for x in range(n_ts_e)])))
    pcolors += list(df_t_se['color'])

    xs += list(
        5*df_ts['type'] + df_ts['sepsis_at_admission'] + 0.5 -
        3*pad + 4*pad*df_ts['sepsis_at_admission'])
    ys += [(2*(n_rows - 3) - pad) for x in range(n_ts)]
    widths += [(1 + 2*pad) for x in range(n_ts)]
    heights += [(1 + 2*pad) for x in range(n_ts)]
    text_xs += list(
        5*df_ts['type'] + df_ts['sepsis_at_admission'] + 1 -
        2*pad + 4*pad*df_ts['sepsis_at_admission'])
    text_ys += [2*(n_rows - 3) + 0.5 for x in range(n_ts)]
    texts += list(df_ts['text'])

    cxs += list(
        5*df_ts_e['type'] + df_ts_e['sepsis_at_admission'] + df_ts_e['cxs'] +
        0.5 - 2*pad + 4*pad*df_ts_e['sepsis_at_admission'])
    cys += [2*(n_rows - 3) for x in range(n_ts_e)]
    cwidths += list(df_ts_e['cwidths'])
    cheights += [1 for x in range(n_ts_e)]
    ccolors += list(df_ts_e['color'])

    reorder = df_te_s.sort_values(by=['type', 'sepsis_at_admission']).index
    polys += list(zip(
        zip(5*df_ts_e['type'] + df_ts_e['sepsis_at_admission'] +
            df_ts_e['cxs'] + 0.5 - 2*pad +
            4*pad*df_ts_e['sepsis_at_admission'],
            [2*(n_rows - 3) for x in range(n_ts_e)]),
        zip(5*df_ts_e['type'] + df_t_se['sepsis_at_admission'] +
            df_ts_e['cxs'] + df_ts_e['cwidths'] + 0.5 - 2*pad +
            4*pad*df_ts_e['sepsis_at_admission'],
            [2*(n_rows - 3) for x in range(n_ts_e)]),
        zip((5*df_te_s['type'] + df_te_s['n_episodes'] + df_te_s['cxs'] -
            3*pad + 3*pad*df_te_s['n_episodes'] + df_te_s['cwidths'])[reorder],
            [2*(n_rows - 4) + 1 for x in range(n_te_s)]),
        zip((5*df_te_s['type'] + df_te_s['n_episodes'] + df_te_s['cxs'] -
            3*pad + 3*pad*df_te_s['n_episodes'])[reorder],
            [2*(n_rows - 4) + 1 for x in range(n_te_s)])))
    pcolors += list(df_ts_e['color'])

    xs += list(
        5*df_te['type'] + df_te['n_episodes'] -
        4*pad + 3*pad*df_te['n_episodes'])
    ys += [(2*(n_rows - 4) - pad) for x in range(n_te)]
    widths += [(1 + 2*pad) for x in range(n_te)]
    heights += [(1 + 2*pad) for x in range(n_te)]
    text_xs += list(
        5*df_te['type'] + df_te['n_episodes'] + 0.5 -
        4*pad + 3*pad*df_te['n_episodes'])
    text_ys += [2*(n_rows - 4) + 0.5 for x in range(n_te)]
    texts += list(df_te['text'])

    cxs += list(
        5*df_te_s['type'] + df_te_s['n_episodes'] + df_te_s['cxs'] -
        3*pad + 3*pad*df_te_s['n_episodes'])
    cys += [2*(n_rows - 4) for x in range(n_te_s)]
    cwidths += list(df_te_s['cwidths'])
    cheights += [1 for x in range(n_te_s)]
    ccolors += list(df_te_s['color'])

    reorder = df_td_es.sort_values(by=['type', 'n_episodes']).index
    polys += list(zip(
        zip(5*df_te_sd['type'] + df_te_sd['n_episodes'] +
            df_te_sd['cxs'] - 3*pad + 3*pad*df_te_sd['n_episodes'],
            [2*(n_rows - 4) for x in range(n_te_sd)]),
        zip(5*df_te_sd['type'] + df_te_sd['n_episodes'] + df_te_sd['cxs'] -
            3*pad + 3*pad*df_te_sd['n_episodes'] + df_te_sd['cwidths'],
            [2*(n_rows - 4) for x in range(n_te_sd)]),
        zip((5*df_td_es['type'] + 1.5*df_td_es['destination'] - 0.75 +
            1.5*df_td_es['cxs'] + 1.5*df_td_es['cwidths'] +
            4*pad + 3*pad*df_td_es['destination'])[reorder],
            [2*(n_rows - 5) + 1 for x in range(n_td_es)]),
        zip((5*df_td_es['type'] + 1.5*df_td_es['destination'] - 0.75 +
            1.5*df_td_es['cxs'] + 4*pad +
            3*pad*df_td_es['destination'])[reorder],
            [2*(n_rows - 5) + 1 for x in range(n_td_es)])))
    pcolors += list(df_te_sd['color'])

    xs += list(
        5*df_td['type'] + 1.5*df_td['destination'] - 0.75 +
        4*pad + 3*pad*df_td['destination'])
    ys += [2*(n_rows - 5) for x in range(n_td)]
    widths += [1.5 for x in range(n_td)]
    heights += [1 for x in range(n_td)]
    text_xs += list(
        5*df_td['type'] + 1.5*df_td['destination'] +
        4*pad + 3*pad*df_td['destination'])
    text_ys += [2*(n_rows - 5) + 0.5 for x in range(n_td)]
    texts += list(df_td['text'])

    cxs += list(
        5*df_td_es['type'] + 1.5*df_td_es['destination'] - 0.75 +
        1.5*df_td_es['cxs'] + 4*pad + 3*pad*df_td_es['destination'])
    cys += [2*(n_rows - 5) for x in range(n_td_es)]
    cwidths += list(1.5*df_td_es['cwidths'])
    cheights += [1 for x in range(n_td_es)]
    ccolors += list(df_td_es['color'])

    ylims = [y + height for y, height in zip(ys, heights)]
    # ax.set_ylim([min(ys), max(ylims)])
    xlims = [x + width for x, width in zip(xs, widths)]
    # ax.set_xlim([min(xs), max(xlims)])

    text_xs += [min(xlims) - 3] * 5
    text_ys += [6.5, 4.5, 2.5, 0.5, 9.2]
    texts += ['Admission type\n(n)', 'Sepsis at admission\n(n)']
    texts += ['Number of sepsis episodes\n(n)', 'Outcome\n(n)', 'Key']

    text_key_xs = [min(xlims) - 3.9] * 4
    text_key_xs += [min(xlims) - 3.45 + 0.4*x for x in range(5)]
    text_key_ys = [8.7, 8.2, 7.8, 7.4, 8.7, 8.7, 8.7, 8.7, 8.7]
    texts_key = ['Sepsis at admission\nNumber of sepsis episodes']
    texts_key += ['Elective surgery', 'Emergency medical', 'Emergency surgery']
    texts_key += ['No\n0 ', 'No\n1 ', 'No\n>1 ', 'Yes\n1 ', 'Yes\n>1 ']

    key_cxs = list(
        min(xlims) - 3.8 + 0.4*df_tse['n_episodes'] +
        0.8*df_tse['sepsis_at_admission'])
    key_cys = list(8 - 0.4*df_tse['type'])
    key_cwidths = [0.4 for x in range(n_tse)]
    key_cheights = [0.4 for x in range(n_tse)]
    key_ccolors = list(df_tse['color'])

    figure_name = unit_type + '/' + 'sepsis_flowchart_' + unit_type + '.pdf'
    with PdfPages(inputs.output_file_path + figure_name) as pdf:
        fig, ax = plt.subplots(figsize=[12, 8])
        ax.set_aspect('equal')
        ax.set_xlim([min(xs) - 0.1 - 3, max(xlims) + 0.1])
        ax.set_ylim([min(ys) - 0.1, max(ylims) + 0.1])

        c_rectangles = [
            Rectangle((x, y), width, height, facecolor=fc, lw=0, alpha=0.8)
            for x, y, width, height, fc in zip(
                cxs, cys, cwidths, cheights, ccolors)]

        for c_rectangle in c_rectangles:
            ax.add_patch(c_rectangle)

        k_rectangles = [
            Rectangle((x, y), width, height, facecolor=fc, lw=0.5, alpha=0.8)
            for x, y, width, height, fc in zip(
                key_cxs, key_cys, key_cwidths, key_cheights, key_ccolors)]

        for k_rectangle in k_rectangles:
            ax.add_patch(k_rectangle)

        poly_list = [
            Polygon(
                p, facecolor=fc, closed=True,
                edgecolor='black', lw=0.5, alpha=0.6)
            for p, fc in zip(polys, pcolors)]

        for poly in poly_list[::-1]:
            ax.add_patch(poly)

        rectangles = [
            Rectangle((x, y), width, height)
            for x, y, width, height in zip(xs, ys, widths, heights)]

        pc = PatchCollection(
            rectangles, facecolor=(1, 1, 1, 0), edgecolor='black', lw=2)
        ax.add_collection(pc)

        for text_x, text_y, text in zip(text_xs, text_ys, texts):
            ax.text(
                text_x, text_y, text,
                weight='bold',
                horizontalalignment='center', verticalalignment='center')

        for text_x, text_y, text in zip(text_key_xs, text_key_ys, texts_key):
            ax.text(
                text_x, text_y, text,
                weight='bold', fontsize='x-small',
                horizontalalignment='right', verticalalignment='center')

        ax.axis('off')
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close()
