#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  2 17:46:42 2021

@author: willrichardson
"""

# this script creates dataframes of half-hourly outputs from your standard flux processing program
# and saves them to files, as these will be used in intermediate processing steps
# and for creating final program outputs:
    # (1) All half-hourly periods (needed for final outputs)
    # (2) Only best quality half-hour periods:
        # qc flag 0 for CH4 and LE ("0-1-2" system, Mauder & Foken 2004)
        # u star threshold 
        # fetch requirements (wind direction, 90% cumulative flux contribution fetch distance)
        
# Also, save half-hourly sigma_m data into a text file for easy loading later when doing partitioning

#%% import libraries

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator, MultipleLocator
from funcs import filt, get_UnivFunc
# import column header variables from funcs script
from funcs import c1, c2, bound_const, ust, sigma_m, MOST_stab, sigma_q, sigma_c, wq_cov, wc_cov, HH_tstamp_beg, datetime_sep, tstamp_cols, HH_output_header, HH_output_skiprow, NaN_vals

#%% get relevant environment variables

# get working directory as base path
base_path = os.getcwd()

# filepaths
hh_main_path = os.environ['HH_output']
hh_anc_path = os.environ['HH_stats']
site_yr = os.environ['site_yr']
run_ID = os.environ['run_ID']

# make paths for output files
out_path = base_path + '/ref_data/Output/%s/%s_Interm_RefDf' %(site_yr, site_yr)
plot_path = base_path + '/ref_data/Output/%s/%s_' %(site_yr, site_yr)
# half-hourly sigma_m and sigma_m/u*
sigmaCH4_fname = base_path + '/ref_data/Input/%s/%s_sigmaCH4.txt' %(site_yr, site_yr)
sigmaCH4ust_fname = base_path + '/ref_data/Input/%s/%s_sigmaCH4ust.txt' %(site_yr, site_yr)

# filtering parameters
ust_thresh = np.float64(os.environ['ust_thresh'])
winddir_min = np.float64(os.environ['WD_min'])
winddir_max = np.float64(os.environ['WD_max'])
if os.environ['fetch_x'] == 'False':
    fetch_max = False
else:
    fetch_max = np.float64(os.environ['fetch_x'])

#%% load files

# header and skiprows specifications in the below two lines may need to be modififed depending on your file header
HH_main = pd.read_csv(hh_main_path, header=HH_output_header, skiprows=HH_output_skiprow, na_values=NaN_vals)
# Make timestamp, convert to beginning of time period if necessary
if datetime_sep == True:
    HH_main['Timestamp'] = pd.to_datetime(HH_main[[tstamp_cols[0],tstamp_cols[1]]].apply(lambda x: ' '.join(x), axis=1))
    HH_main = HH_main.set_index('Timestamp')
else:
    HH_main = HH_main.set_index(tstamp_cols[0])
    HH_main.index.name = 'Timestamp'
    
# Shift timestamps back 30 minutes to match convention of other workflow steps
if HH_tstamp_beg == False:
    HH_main.index = HH_main.index - pd.Timedelta(30, 'm')
    HH_main['DOY'] = HH_main['DOY'] - (0.5/24)

# repeat steps for ancillary dataframe if one exists
if hh_anc_path != None:
    HH_anc = pd.read_csv(hh_anc_path, header='infer', na_values=[-9999, -999900])
    # Make timestamp, convert to beginning of time period if necessary
    if datetime_sep == True:
        HH_anc['Timestamp'] = pd.to_datetime(HH_anc[[tstamp_cols[0],tstamp_cols[1]]].apply(lambda x: ' '.join(x), axis=1))
        HH_anc = HH_anc.set_index('Timestamp')
    else:
        HH_anc = HH_anc.set_index(tstamp_cols[0])
        HH_anc.index.name = 'Timestamp'
    
    # Shift timestamps back 30 minutes to match convention of other workflow steps
    if HH_tstamp_beg == False:
        HH_anc.index = HH_anc.index - pd.Timedelta(30, 'm')    
        
    # get full dataframe
    HH_full = HH_main.join(HH_anc, how='left', rsuffix='_stats')

else:
    HH_full = HH_main

## calculate sigma_m/u* for an ebullition threshold quantity option
HH_full.loc[:, 'sigmaCH4_ust'] = HH_full.loc[:, sigma_m]/HH_full.loc[:, ust]

#%% save full dataframe

HH_full.to_csv(out_path + '_allHH.csv')

#%% do filtering
HH_qPart_filt, HH_cPart_filt, HH_ch4_filt = filt(HH_full, ust_thresh, winddir_min, winddir_max, fetch_max)
    
#%% plot and save distribution of u* values in each filtered dataframe

fig, ax = plt.subplots(nrows=1, ncols=2, sharex=True, sharey=True, figsize=(5,2.2))

ax[0].hist(HH_qPart_filt[ust], bins=40, color='mediumaquamarine', ec='k', lw=0.5)
ax0r = ax[0].twinx()
n_0, bins_0, patches_0 = ax0r.hist(HH_qPart_filt[ust], bins=1000, density=True, cumulative=True, histtype='step', color='goldenrod', edgecolor='goldenrod', 
                                   lw=0.9, zorder=4)
patches_0[0].set_xy(patches_0[0].get_xy()[:-1])
ax0r.tick_params(axis='y', which='major', labelright=False)
ax0r.set_ylim((0,1.005))

ax[1].hist(HH_cPart_filt[ust], bins=40, color='mediumaquamarine', ec='k', lw=0.5)
ax1r = ax[1].twinx()
n_1, bins_1, patches_1 = ax1r.hist(HH_cPart_filt['u*'], bins=1000, density=True, cumulative=True, histtype='step', color='goldenrod', edgecolor='goldenrod', 
                                   lw=0.9, zorder=4)
patches_1[0].set_xy(patches_1[0].get_xy()[:-1])
ax1r.tick_params(which='major', axis='y', labelsize=7)
ax1r.set_ylim((0,1.005))

ax[0].set_title('best quality periods\nfor partitioning based on $q$, n = %i' %HH_qPart_filt.shape[0], size=7)
ax[1].set_title('best quality periods\nfor partitioning based on $c$, n = %i' %HH_cPart_filt.shape[0], size=7)

for i in np.arange(ax.size):
    ax.reshape(-1)[i].set_xlabel('$u_*$ (m s$^{-1}$)', size=8)
    ax.reshape(-1)[i].xaxis.set_major_locator(MultipleLocator(0.1))    
    ax.reshape(-1)[i].xaxis.set_minor_locator(AutoMinorLocator())
    ax.reshape(-1)[i].yaxis.set_minor_locator(AutoMinorLocator())
    ax.reshape(-1)[i].grid(which='major', axis='both', alpha=0.3)
    ax.reshape(-1)[i].tick_params(which='major', axis='both', labelsize=7)
    ax.reshape(-1)[i].set_axisbelow(True)
        
ax[0].set_ylabel('Count', size=8)
ax1r.set_ylabel('CDF', size=8)
ax[1].tick_params(axis='y', which='major', labelleft=False)

fig.suptitle('%s $u_*$ distribution' %site_yr, size=8, y=1.1)
fig.subplots_adjust(wspace=0.08)
plt.savefig(plot_path + 'ustar_distribution.png', dpi=300, bbox_inches='tight')

#%% save half-hourly sigma_m data (and sigma_m/ust) into a text file for easy loading later

filenames = [x.strftime('%Y%m%d') + '_' + x.strftime('%H%M') for x in HH_full.index]
sigmaCH4_df = HH_full.copy(); sigmaCH4_ust_df = HH_full.copy()
sigmaCH4_df.insert(0, 'file_name', filenames); sigmaCH4_ust_df.insert(0, 'file_name', filenames)
# drop everything but quantity of interest
sigmaCH4_out = sigmaCH4_df[['file_name', sigma_m]]
sigmaCH4ust_out = sigmaCH4_ust_df[['file_name', 'sigmaCH4_ust']]
# save files
sigmaCH4_out.to_csv(sigmaCH4_fname, sep='\t', header=False, index=False)
sigmaCH4ust_out.to_csv(sigmaCH4ust_fname, sep='\t', header=False, index=False)
# 
#%% filter for non-local processes based on flux variance similarity in reference scalars

# dimensionless flux-variance quantity
HH_qPart_filt['q_MOST_scale'] = HH_qPart_filt[sigma_q]/(np.abs(-1*HH_qPart_filt[wq_cov]/HH_qPart_filt[ust]))
HH_cPart_filt['c_MOST_scale'] = HH_cPart_filt[sigma_c]/(np.abs(-1*HH_cPart_filt[wc_cov]/HH_cPart_filt[ust]))

# stability range for universal function
x_universal = np.arange(-4.9, 0.53, 0.001)

# get universal functions (De Bruin et al 1993)
univ_func = np.empty(len(x_universal), dtype=np.float64)
for i in np.arange(len(x_universal)):
    if x_universal[i] > 0:
        univ_func[i] = 2
    else:
        univ_func[i] = c1*(1 - c2*x_universal[i])**(-1/3)

# get boundaries for filtering
univ_func_low = univ_func/bound_const
univ_func_up = univ_func*bound_const

#%% filter for points falling outside of bounds

MOST_phi_q = get_UnivFunc(HH_qPart_filt)
HH_qPart_filt = HH_qPart_filt.join(MOST_phi_q)
HH_qNonLocal = HH_qPart_filt.loc[(HH_qPart_filt['q_MOST_scale'] < HH_qPart_filt['MOST_phi_lower']) | (HH_qPart_filt['q_MOST_scale'] > HH_qPart_filt['MOST_phi_upper'])]
HH_qLocal = HH_qPart_filt.loc[(HH_qPart_filt['q_MOST_scale'] >= HH_qPart_filt['MOST_phi_lower']) & (HH_qPart_filt['q_MOST_scale'] <= HH_qPart_filt['MOST_phi_upper'])]

MOST_phi_c = get_UnivFunc(HH_cPart_filt)
HH_cPart_filt = HH_cPart_filt.join(MOST_phi_c)
HH_cNonLocal = HH_cPart_filt.loc[(HH_cPart_filt['c_MOST_scale'] < HH_cPart_filt['MOST_phi_lower']) | (HH_cPart_filt['c_MOST_scale'] > HH_cPart_filt['MOST_phi_upper'])]
HH_cLocal = HH_cPart_filt.loc[(HH_cPart_filt['c_MOST_scale'] >= HH_cPart_filt['MOST_phi_lower']) & (HH_cPart_filt['c_MOST_scale'] <= HH_cPart_filt['MOST_phi_upper'])]

#%% visualize relationships
fig, ax = plt.subplots(nrows=1, ncols=2, sharex=True, sharey=True, figsize=(5,2.2))

# universal funcs
for i in np.arange(ax.size):
    ax[i].plot(-1*x_universal, univ_func, '-', c='red', lw=1.0)
    ax[i].plot(-1*x_universal, univ_func_low, '--', c='red', lw=1.0)
    ax[i].plot(-1*x_universal, univ_func_up, '--', c='red', lw=1.0)

# plot data
ax[0].scatter(-1*HH_qPart_filt[MOST_stab], HH_qPart_filt['q_MOST_scale'], facecolors='none',
              edgecolors='k', alpha=0.3, s=10, linewidths=0.7, label='n = %i'%HH_qPart_filt.shape[0])
ax[0].scatter(-1*HH_qNonLocal[MOST_stab], HH_qNonLocal['q_MOST_scale'], facecolors='none',
              edgecolors='r', alpha=0.8, s=10, linewidths=0.7, label='n = %i removed\n(non-local processes)'%HH_qNonLocal.shape[0])

ax[1].scatter(-1*HH_cPart_filt[MOST_stab], HH_cPart_filt['c_MOST_scale'], facecolors='none',
              edgecolors='k', alpha=0.3, s=10, linewidths=0.7, label='n = %i'%HH_cPart_filt.shape[0])
ax[1].scatter(-1*HH_cNonLocal[MOST_stab], HH_cNonLocal['c_MOST_scale'], facecolors='none',
              edgecolors='r', alpha=0.8, s=10, linewidths=0.7, label='n = %i removed'%HH_cNonLocal.shape[0])

ax[0].set_ylabel('$\sigma_q$/|$q_*$|', size=9)
ax[1].set_ylabel('$\sigma_c$/|$c_*$|', size=9, labelpad=8)

ax[0].set_xscale('log')
ax[0].set_yscale('log')
ax[1].set_yscale('log')

for i in np.arange(ax.size):
    ax.reshape(-1)[i].set_xlabel('- $\\frac{(z-d)}{L}$', size=11)
    ax.reshape(-1)[i].grid(which='major', axis='both', alpha=0.5)
    ax.reshape(-1)[i].tick_params(which='major', axis='both', labelsize=7)
    ax.reshape(-1)[i].set_axisbelow(True)
    ax.reshape(-1)[i].legend(loc='lower left', bbox_to_anchor=(0.01, 1.02), fontsize=7)
    
fig.subplots_adjust(hspace=0.25)

ax[0].set_ylim((0.3,50))

plt.savefig(plot_path + 'MOST_FV_plot.png', dpi=400, bbox_inches='tight')

#%% save filtered dfs

HH_qPart_filt.to_csv(out_path + '_qPart_filtHH.csv')
HH_cPart_filt.to_csv(out_path + '_cPart_filtHH.csv')
HH_qLocal.to_csv(out_path + '_qPart_filtHH_NLfilt.csv')
HH_cLocal.to_csv(out_path + '_cPart_filtHH_NLfilt.csv')
