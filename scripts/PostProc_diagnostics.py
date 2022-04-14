#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  9 17:10:51 2022

@author: willrichardson
"""
# this script takes the full partitioning program output and creates some initial figures
# which are useful for determining fn_LB and other program parameters/settings

#%% import libraries

import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, AutoMinorLocator, LogLocator, NullFormatter
import matplotlib.lines as Lines

#%% load data, declare other parameters

site_yr = os.environ['site_yr']
base_path = os.getcwd()
run_ID = os.environ['run_ID']
fn_LB = np.float64(os.environ['LB'])

save_path = base_path + '/ref_data/Output/' + site_yr 
qPart_path = base_path + '/ref_data/Output/' + site_yr + '/' + run_ID + '/' + site_yr + '_FullOutput_qPart_filtHH_NLfilt.csv'
cPart_path = base_path + '/ref_data/Output/' + site_yr + '/' + run_ID + '/' + site_yr + '_FullOutput_cPart_filtHH_NLfilt.csv'

qPart_df = pd.read_csv(qPart_path, index_col='Timestamp', parse_dates=True, infer_datetime_format=True)
cPart_df = pd.read_csv(cPart_path, index_col='Timestamp', parse_dates=True, infer_datetime_format=True)

# make dataframe where both reference scalars are good
qc_both_tstamps = qPart_df.join(cPart_df, how='inner', rsuffix='_cPart').index
RefScal_filt_df = qPart_df.loc[qc_both_tstamps]

#%% helper functions and routines for making figures 

# function for getting bin averages
def BinAvg_2(df, fn_str, var_str, nBins):
    # get variables 
    if var_str in ['mq_coh1', 'mc_coh1', 'mT_coh1']:
        # print('get rid of Non-ebullitive coherences')
        bad_label = 'NonEb_' + var_str
        var2d_all = df.filter(like=var_str).drop(list(df.filter(like=bad_label).columns), axis=1)
    else:
        var2d_all = df.filter(like=var_str)
    
    fn2d_all = df.filter(like=fn_str)
    
    # for coherences, drop 2 coarsest scales/lowest frequencies
    if 'coh' in var_str:
        fn2d_all = fn2d_all.drop(labels=[fn2d_all.columns[0], fn2d_all.columns[1]], axis=1)
    
    # concatenate these in a new frame
    var_all = pd.DataFrame(columns=['fn', var_str])
    for i in np.arange(var2d_all.shape[0]):
        var_hh = var2d_all.iloc[i,:]
        fn_hh = fn2d_all.iloc[i,:]
        intermediate_df = pd.DataFrame([fn_hh.values, var_hh.values]).swapaxes(0,1).rename({0: 'fn', 1: var_str}, axis=1)
        var_all = var_all.append(intermediate_df)
    
    # if variable is qc, mc, or Tm coherence, take the absolute value of the coherences 
    if 'mc_coh' or 'mT_coh' in var_str:
        var_all[var_str] = np.abs(var_all[var_str])
    
    # sort new frame by ascending frequency
    fn_sort = var_all.sort_values(by='fn')
    
    sort_split = np.array_split(fn_sort, nBins, axis=0)
    
    # get bin stats
    bin_means = pd.DataFrame(columns=fn_sort.columns)
    bin_std = pd.DataFrame(columns=fn_sort.columns)
    for frame in sort_split:
        
        # prevent extreme values from impacting stats
        if 'NormCov' in var_str:
            thresh = np.percentile(frame.loc[:,var_str], 98)
            frame = frame.loc[frame[var_str] < thresh]  
            
        mean = frame.mean(axis=0)
        bin_means = pd.concat([bin_means, mean], ignore_index=True, sort=True)
        std = frame.std(axis=0)
        bin_std = pd.concat([bin_std, std], ignore_index=True, sort=True)
    
    return(bin_means, bin_std)

#%% Examine how fn_LB value impacts inclusion of Wx scales

# desired range of fn_LBs to examine
LBs = np.arange(0.000, 0.0305, 0.0005)

ts0p43_out_fnzd = []; ts0p85_out_fnzd = []; ts1p7_out_fnzd = []; ts3p4_out_fnzd = []; ts6p8_out_fnzd = []; ts13p7_out_fnzd = []; 
for i in np.arange(len(LBs)):

    fnzd_j9_ngt_LB_new = qPart_df.loc[qPart_df['fn_j9_ts25.6'] >= LBs[i]]
    ts0p43_out_fnzd.append((fnzd_j9_ngt_LB_new.shape[0]/qPart_df.shape[0])*100) 
    
    fnzd_j10_ngt_LB_new = qPart_df.loc[qPart_df['fn_j10_ts51.2'] >= LBs[i]]
    ts0p85_out_fnzd.append((fnzd_j10_ngt_LB_new.shape[0]/qPart_df.shape[0])*100)     
    
    fnzd_j11_ngt_LB_new = qPart_df.loc[qPart_df['fn_j11_ts102.4'] >= LBs[i]]
    ts1p7_out_fnzd.append((fnzd_j11_ngt_LB_new.shape[0]/qPart_df.shape[0])*100)    
    
    fnzd_j12_ngt_LB_new = qPart_df.loc[qPart_df['fn_j12_ts204.8'] >= LBs[i]]
    ts3p4_out_fnzd.append((fnzd_j12_ngt_LB_new.shape[0]/qPart_df.shape[0])*100)
    
    fnzd_j13_ngt_LB_new = qPart_df.loc[qPart_df['fn_j13_ts409.6'] >= LBs[i]]
    ts6p8_out_fnzd.append((fnzd_j13_ngt_LB_new.shape[0]/qPart_df.shape[0])*100)
    
    fnzd_j14_ngt_LB_new = qPart_df.loc[qPart_df['fn_j14_ts819.2'] >= LBs[i]]
    ts13p7_out_fnzd.append((fnzd_j14_ngt_LB_new.shape[0]/qPart_df.shape[0])*100)

ts0p43_out_fnzd = np.asarray(ts0p43_out_fnzd)
ts0p85_out_fnzd = np.asarray(ts0p85_out_fnzd)
ts1p7_out_fnzd = np.asarray(ts1p7_out_fnzd)
ts3p4_out_fnzd = np.asarray(ts3p4_out_fnzd)
ts6p8_out_fnzd = np.asarray(ts6p8_out_fnzd)
ts13p7_out_fnzd = np.asarray(ts13p7_out_fnzd)

#%% Get spectral stats, plot

nBin = 30
# get bin averages
wm_cov_means, wm_cov_std = BinAvg_2(qPart_df, 'fn_j', 'wm_NormCov', nBin)
mq_coh_means, mq_coh_std = BinAvg_2(qPart_df, 'fn_j', 'mq_coh1', nBin)
mc_coh_means, mc_coh_std = BinAvg_2(cPart_df, 'fn_j', 'mc_coh1', nBin)
qc_coh_means, qc_coh_std = BinAvg_2(RefScal_filt_df, 'fn_j', 'qc_coh1', nBin)

# column lists for easy filtering later
fn_cols = list(qPart_df.filter(like='fn_j').columns); fn_coh_cols = fn_cols[2:]

# fnLB and fnUB values to plot as vlines
fnLB_plot = 0.003
fnUB_plot = 1.0

#%% plot

fig, ax = plt.subplots(nrows=2, ncols=2, sharex=True, sharey=False, figsize=(6.8,5.5))

# plot half-hourly values
ax[0,0].scatter(qPart_df[fn_cols], qPart_df.filter(like='wm_Norm'),
                facecolors='none', edgecolors='k', alpha=0.2, zorder=0, s=11, linewidths=0.8)
ax[0,1].scatter(RefScal_filt_df[fn_coh_cols], RefScal_filt_df.filter(like='qc_coh1'),
                facecolors='none', edgecolors='k', alpha=0.2, zorder=0, s=11, linewidths=0.8)
ax[1,0].scatter(qPart_df[fn_coh_cols], qPart_df.filter(like='mq_coh1'),
                facecolors='none', edgecolors='k', alpha=0.2, zorder=0, s=11, linewidths=0.8)
ax[1,1].scatter(cPart_df[fn_coh_cols], cPart_df.filter(like='mc_coh1'),
                facecolors='none', edgecolors='k', alpha=0.2, zorder=0, s=11, linewidths=0.8)

# add binned values
ax[0,0].scatter(wm_cov_means['fn'], wm_cov_means['wm_NormCov'], facecolors='grey', edgecolors='grey', zorder=2, s=13)
ax[0,0].errorbar(wm_cov_means['fn'], wm_cov_means['wm_NormCov'], yerr=wm_cov_std['wm_NormCov'], 
                  xerr=None, ecolor='grey', elinewidth=1, capsize=1, ls='none')
ax[0,1].scatter(qc_coh_means['fn'], -1*qc_coh_means['qc_coh1'], facecolors='lightseagreen', edgecolors='lightseagreen', zorder=2, s=13)
ax[0,1].errorbar(qc_coh_means['fn'], -1*qc_coh_means['qc_coh1'], yerr=qc_coh_std['qc_coh1'], 
                  xerr=None, ecolor='lightseagreen', elinewidth=1, capsize=1, ls='none')
ax[1,0].scatter(mq_coh_means['fn'], mq_coh_means['mq_coh1'], facecolors='royalblue', edgecolors='royalblue', zorder=2, s=13)
ax[1,0].errorbar(mq_coh_means['fn'], mq_coh_means['mq_coh1'], yerr=mq_coh_std['mq_coh1'], 
                  xerr=None, ecolor='royalblue', elinewidth=1, capsize=1, ls='none')
ax[1,1].scatter(mc_coh_means['fn'], -1*mc_coh_means['mc_coh1'], facecolors='limegreen', edgecolors='limegreen', zorder=2, s=13)
ax[1,1].errorbar(mc_coh_means['fn'], -1*mc_coh_means['mc_coh1'], yerr=mc_coh_std['mc_coh1'], 
                  xerr=None, ecolor='limegreen', elinewidth=1, capsize=1, ls='none')

# format and apply labels
ax[0,0].set_xscale('log'); ax[0,0].set_xlim((1.3e-4, 90))
ax[0,0].set_xticks([10**-3, 10**-2, 10**-1, 10**0, 10**1], minor=False)
ax[0,0].xaxis.set_minor_locator(LogLocator(subs=(0.1,)))

ax[1,0].set_xlabel(r'$\frac{z_m - z_d}{t_s*u}$', size=11, labelpad=2)
ax[1,1].set_xlabel(r'$\frac{z_m - z_d}{t_s*u}$', size=11, labelpad=2)

ax[0,0].set_ylabel('normalized wm covariance', size=10)
ax[0,1].set_ylabel('qc coherence', size=10)
ax[1,0].set_ylabel('mq coherence', size=10)
ax[1,1].set_ylabel('mc coherence', size=10)

ax[0,0].set_ylim((-0.2, 0.4)); ax[0,0].yaxis.set_major_locator(MultipleLocator(0.1)); ax[0,0].yaxis.set_minor_locator(AutoMinorLocator())

for i in np.arange(ax.size):
    if i != 0:
        ax.reshape(-1)[i].set_ylim((-1.05, 1.05))
        ax.reshape(-1)[i].yaxis.set_major_locator(MultipleLocator(0.25))
        ax.reshape(-1)[i].yaxis.set_minor_locator(AutoMinorLocator())
    
    ax.reshape(-1)[i].grid(which='major', axis='both', alpha=0.45)
    ax.reshape(-1)[i].grid(which='minor', axis='both', alpha=0.15)
    ax.reshape(-1)[i].tick_params(axis='both', which='major', labelsize=7)

    x_minor = LogLocator(base = 10.0, subs = np.arange(1.0, 10.0) * 0.1, numticks=10)
    ax.reshape(-1)[i].xaxis.set_minor_locator(x_minor)
    ax.reshape(-1)[i].yaxis.set_minor_formatter(NullFormatter())

    ax.reshape(-1)[i].axvline(x=fnLB_plot, ls='--', color='gold', lw=1.1, zorder=1, label='f$_n$ default boundaries')
    ax.reshape(-1)[i].axvline(x=fnUB_plot, ls='--', color='gold', lw=1.1, zorder=1)
    ax.reshape(-1)[i].axhline(y=0, ls=':', color='k', zorder=0, lw=1.1)
    ax.reshape(-1)[i].set_axisbelow(True)

ax[0,0].legend(loc='lower center', bbox_to_anchor=(0.5,1.01), fontsize=8)

fig.subplots_adjust(hspace=0.08, wspace=0.32)
fig.suptitle('%s Spectral Plots'%site_yr, size=10, y=0.98)
plt.savefig(save_path + '/SpectralPlots.png', dpi=400, bbox_inches='tight')

#%% Make figure visualizing inclusion of Wx scales across fnLB values

fig, ax = plt.subplots(nrows=1, ncols=1, sharey=True, figsize=(4.3,3))

ax.plot(LBs, ts0p43_out_fnzd, '-', c='darkviolet', label='25.6 s', lw=1.2)
ax.plot(LBs, ts0p85_out_fnzd, '-', c='mediumblue', label='51.2 s', lw=1.2)
ax.plot(LBs, ts1p7_out_fnzd, '-', c='g', label='1.7 min', lw=1.2)
ax.plot(LBs, ts3p4_out_fnzd, '-', c='gold', label='3.4 min', lw=1.2)
ax.plot(LBs, ts6p8_out_fnzd, '-', c='orange',  label='6.8 min', lw=1.2)
ax.plot(LBs, ts13p7_out_fnzd, '-', c='crimson', label='13.7 min', lw=1.2)

ax.set_xlabel('$f_{n,LB}$', size=9, labelpad=3)
ax.set_ylabel('% of periods including $W_x$ in partitioning', size=9)
ax.axvline(0.003, ls=':', c='darkred', zorder=-1, lw=1.2)
ax.legend(loc='lower center', bbox_to_anchor=(0.5, 0.99), ncol=3, fontsize=8, handlelength=1, columnspacing=1)
ax.grid(axis='both', which='major', alpha=0.3)

ax.set_axisbelow(True)
ax.yaxis.set_minor_locator(MultipleLocator(5))
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.grid(axis='y', which='minor', alpha=0.15)
ax.set_xlim((0.0,0.0302))
ax.set_ylim((0,101))

ax.tick_params(which='major', axis='both', labelsize=7, pad=2)
ax.tick_params(which='major', axis='x', labelrotation=15, pad=1, length=4)

# extend vertical line signifying default fn_LB below subplot and annotate it
vline = Lines.Line2D([0.203,0.203], [0.02,0.12], lw=1.2, transform=fig.transFigure, figure=fig, ls=':', c='darkred')
fig.lines.extend([vline])
ax.annotate('default value', (0.09,-0.158), xycoords='axes fraction', size=8, c='darkred', ha='center')

# plot empirically selected value on top (uncomment below lines if you want to visualize an updated guess)
# ax.axvline(fn_LB, ls=':', c='darkturquoise', zorder=-1, lw=1.2)
# ax.annotate('empirical\nvalue', (0.0175,66), xycoords='data', size=8, c='darkturquoise', ha='center')

fig.suptitle('%s'%site_yr, size=10, y=1.08)

plt.savefig(save_path + '/fnLB_WxScale_inclusion.png', dpi=400, bbox_inches='tight')

