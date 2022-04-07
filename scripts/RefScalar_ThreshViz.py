#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 22 12:36:11 2021

@author: willrichardson
"""

# visualize reference scalar flux thresholds for adequate similarity, calculated in
# RefScalar_Thresh.R

#%% import libraries

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
from funcs import getBinAvg
from funcs import LE_name

#%% environment variables

base_path = os.getcwd()
site_yr = os.environ['site_yr']
run_ID = os.environ['run_ID']
user_LEthresh = np.float64(os.environ['LE_thresh'])
user_Fco2thresh = np.float64(os.environ['Fco2_thresh'])
NL_filt = bool(os.environ['NL_filt'])

#%% load master dataframe with all stats & qc'd frames; filter for low CH4 flux and negative ebullition
# master
master_df = pd.read_csv(base_path + '/ref_data/Output/%s/%s/%s_Interm_RefDf_allHH_PostProcStats.csv'%(site_yr, run_ID, site_yr),
                        index_col=0, infer_datetime_format=True, parse_dates=True)

if NL_filt == True:
    qPart_filt = pd.read_csv(base_path + '/ref_data/Output/%s/%s_Interm_RefDf_qPart_filtHH_NLfilt.csv' %(site_yr, site_yr),
                          index_col='Timestamp', infer_datetime_format=True, parse_dates=True)
    
    cPart_filt = pd.read_csv(base_path + '/ref_data/Output/%s/%s_Interm_RefDf_cPart_filtHH_NLfilt.csv' %(site_yr, site_yr),
                          index_col='Timestamp', infer_datetime_format=True, parse_dates=True)
else:
    qPart_filt = pd.read_csv(base_path + '/ref_data/Output/%s/%s_Interm_RefDf_qPart_filtHH.csv' %(site_yr, site_yr),
                          index_col='Timestamp', infer_datetime_format=True, parse_dates=True)
    
    cPart_filt = pd.read_csv(base_path + '/ref_data/Output/%s/%s_Interm_RefDf_cPart_filtHH.csv' %(site_yr, site_yr),
                          index_col='Timestamp', infer_datetime_format=True, parse_dates=True) 

# q; use indices of previously filtered dataframe to keep desired rows from master_df
qPart_df = master_df[master_df.index.isin(qPart_filt.index)]
# repeat for c
cPart_df = master_df[master_df.index.isin(cPart_filt.index)]

# filter for low CH4 flux and negative ebullition
qPart_df = qPart_df.loc[qPart_df['CH4_tot'] >= 0.01]; cPart_df = cPart_df.loc[cPart_df['CH4_tot'] >= 0.01]
qPart_df = qPart_df.loc[qPart_df['frac_eb_q'] >= 0.0]; cPart_df = cPart_df.loc[cPart_df['frac_eb_c'] >= 0.0] 

#%% load breakpoint parameters

brkpt_path = base_path + '/ref_data/Output/%s/%s/%s_FracEb_RefScalFlux_brkpts.csv'%(site_yr, run_ID, site_yr)
if os.path.exists(brkpt_path):
    brkpt_df = pd.read_csv(brkpt_path)

#%% get bin averages for variable pairs

bin_stats = getBinAvg([qPart_df, cPart_df], [LE_name, 'co2_flux_mag'], ['frac_eb_q', 'frac_eb_c'], 30)

#%% visualize results

fig, ax = plt.subplots(nrows=1, ncols=2, sharey=True, figsize=(5,2.2))

ax[0].scatter(qPart_df[LE_name], qPart_df['frac_eb_q'], s=10, facecolors='none', edgecolors='k',
              linewidths=0.6, alpha=0.5)
ax[0].scatter(bin_stats['0'].iloc[:,0], bin_stats['0'].iloc[:,2], s=15, c='deepskyblue', edgecolors='k',
                              linewidths=0.6, zorder=5)
ax[0].axvline(user_LEthresh, c='forestgreen', ls='--', zorder=3, lw=1.0, label='%.1f W m$^{-2}$ (user-supplied)' %user_LEthresh)
ax[0].set_xlabel('$LE$ (W m$^{-2}$)', size=8, labelpad=2)
ax[0].set_ylabel('$F_{eb}$/$F_{tot}$', size=8, labelpad=2)
ax[0].legend(loc='lower left', bbox_to_anchor=(-0.03, 1.005), fontsize=7)

ax[1].scatter(cPart_df['co2_flux_mag'], cPart_df['frac_eb_c'], s=8, facecolors='none', edgecolors='k',
              linewidths=0.6, alpha=0.5)
ax[1].scatter(bin_stats['1'].iloc[:,0], bin_stats['1'].iloc[:,2], s=15, c='deepskyblue', edgecolors='k',
                              linewidths=0.6, zorder=5)
ax[1].axvline(user_Fco2thresh, c='forestgreen', ls='--', zorder=3, lw=1.0, label='%.1f $\mu$mol m$^{-2}$ s$^{-1}$ (user-supplied)' %user_Fco2thresh)
ax[1].set_xlabel('|$F_{CO_2}$| ($\mu$mol m$^{-2}$ s$^{-1}$)', size=8, labelpad=2)
ax[1].legend(loc='lower left', bbox_to_anchor=(-0.07, 1.005), fontsize=7)

if os.path.exists(brkpt_path):
    ax[0].axvline(np.float64(brkpt_df['FracEbq_LE_brkpt']), c='gold', ls='--', zorder=3, lw=1.0, label='%.1f $\pm$ %.2f W m$^{-2}$\n(Broken-line regression)' %(np.float64(brkpt_df['FracEbq_LE_brkpt']), np.float64(brkpt_df['FracEbq_LE_brkpt_SE'])))
    ax[1].axvline(np.float64(brkpt_df['FracEbc_Fco2_brkpt']), c='gold', ls='--', zorder=3, lw=1.0, label='%.1f $\pm$ %.2f $\mu$mol m$^{-2}$ s$^{-1}$\n(Broken-line regression)' %(np.float64(brkpt_df['FracEbc_Fco2_brkpt']), np.float64(brkpt_df['FracEbc_Fco2_brkpt_SE'])))

for i in np.arange(ax.size):
    ax.reshape(-1)[i].xaxis.set_minor_locator(AutoMinorLocator())
    ax.reshape(-1)[i].yaxis.set_minor_locator(AutoMinorLocator())
    ax.reshape(-1)[i].tick_params(axis='both', which='major', labelsize=7)
    ax.reshape(-1)[i].grid(which='major', axis='both', alpha=0.3, zorder=0)
    ax.reshape(-1)[i].set_axisbelow(True)

ax[0].set_ylim((0,1))
xmax0 = np.percentile(qPart_df[LE_name], 99); xmax1 = np.percentile(cPart_df['co2_flux_mag'], 99)
ax[0].set_xlim(-5,xmax0); ax[1].set_xlim(-1,xmax1)

plt.savefig(base_path + '/ref_data/Output/%s/%s/%s_RefScalarThreshPlot.png' %(site_yr, run_ID, site_yr),
            dpi=300, bbox_inches='tight')
