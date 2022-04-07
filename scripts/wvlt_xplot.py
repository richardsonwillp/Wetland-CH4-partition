#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 14 16:37:59 2022

@author: willrichardson
"""

# this script is written to visualize best quality periods in wavelet domain,
# plotting methane against both reference scalars (q and c) and plotting the
# reference scalars against each other

# run this AFTER completing at least one partitioning run (necessary for annotations)

#%% user variables to change
#   these should generally match the settings specified in part_master.sh 

# site year identifier  
site_yr = 'Way3-Summer-2015'
# partitioning run identifier
run_ID = 'Run2'
# time zones
tz_logging = 'America/Chicago'
tz_local = 'America/Chicago'
# filepaths to your partitioning program, desired location of figures
prog_loc = '/Volumes/Backup-Will/ch4flux_partition_wpr'
fig_path = '/Users/willrichardson/Desktop'

# partitioning  parameters
# interval between data points (see readme for determination based on your acquisition frequency)
DT = 0.05
# lower frequency bound, whatever was used in the partitioning run of interest
fn_LB = 0.015
# upper frequency bound, whatever was used in the partitioning run of interest
fn_UB = 1.0
# measurement height
zm = 2.20
# ebullition threshold sepcifics, for selecting ebullition threshold used in this run
EbThresh_type = 'NonEb_LinReg'
# NOTE this variable should resemble the same variable in part_master.sh, but replace 'm' with CH4 and remove all underscores
EbThresh_qty = 'sigmaCH4'

# specify column names of other variables not listed at top of 'funcs.py'; i.e., mean wind speed, sensible heat flux, air temp in celsius, relative humidity
ws_mean = 'wind_speed'
sens_heat = 'H'
air_temp_C = 'T_air_C'
rel_hum = 'RH'

#%% import libraries
# import sys
import os
# sys.path.insert(0, '/Users/willrichardson/opt/anaconda3/lib/python3.8/site-packages')
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.ticker import AutoMinorLocator
# get column names from funcs.py
from funcs import sigma_m, MOST_stab, LE_name, Fco2_name, WD, ust

#%% Load data, get any other derived quantities

# reference data frames for ancillary information on each period
qPart_df = pd.read_csv(prog_loc + '/ref_data/Output/%s/%s/%s_FullOutput_qPart_filtHH_NLfilt.csv' %(site_yr, run_ID, site_yr),
                       index_col='Timestamp', infer_datetime_format=True, parse_dates=True)
cPart_df = pd.read_csv(prog_loc + '/ref_data/Output/%s/%s/%s_FullOutput_cPart_filtHH_NLfilt.csv' %(site_yr, run_ID, site_yr),
                       index_col='Timestamp', infer_datetime_format=True, parse_dates=True)
full_HH_df = pd.read_csv(prog_loc + '/ref_data/Output/%s/%s/%s_FullOutput_allHH.csv' %(site_yr, run_ID, site_yr),
                       index_col='Timestamp', infer_datetime_format=True, parse_dates=True)
# canopy height data
zc_df = pd.read_csv(prog_loc + 'ref_data/Input/%s/%s_zc.txt' %(site_yr, site_yr),
                    sep='\t', index_col=0, names=['Timestamp', 'zc'], parse_dates=True, infer_datetime_format=True)

# get ebullition threshold widths for this run
EbThresh_widths = pd.read_csv(prog_loc + 'ref_data/Output/%s_RMSDx-sigmaCH4_slopes_EbThreshWidth_fnLB_%s_fnUB_%s.csv' %(site_yr, str(fn_LB), str(fn_UB)),
                              index_col=0)

## get strings for indexing columns of EbThresh_widths dataframe
q_col = 'RMSDq_%s_slope' %EbThresh_qty
c_col = 'RMSDc_%s_slope' %EbThresh_qty

RMSDq_slope = EbThresh_widths.loc[EbThresh_type, q_col]
RMSDc_slope = EbThresh_widths.loc[EbThresh_type, c_col]

# get wavelet transform scale information for fn calcs
nj = int(np.rint(np.log2((30*60)/DT)))
N_wvlt = 2**nj


#%% function for creating a harmonized dataframe
#   i.e., keep estimate based on q when it meets all quality requirements,
#   take estimate based on c if it meets all quality requirements but q does not

def Hrmz_RefScal(full_df, q_df, c_df):
    
    hrmz_df = pd.DataFrame(index=full_df.index, columns=['CH4_eb_hrmnz', 'frac_eb_hrmnz', 'CH4_diff_hrmnz', 'Ebhrmz_rerr_FtotNorm'],
                           dtype=np.float64)
    
    # iterate through periods and make harmonized frame
    for tstamp in hrmz_df.index:
        if tstamp in q_df.index:
            hrmz_df.loc[tstamp, 'CH4_eb_hrmnz'] = q_df.loc[tstamp, 'CH4_eb_q']
            hrmz_df.loc[tstamp, 'frac_eb_hrmnz'] = q_df.loc[tstamp, 'frac_eb_q']
            hrmz_df.loc[tstamp, 'CH4_diff_hrmnz'] = q_df.loc[tstamp, 'CH4_diff_q']
            hrmz_df.loc[tstamp, 'Ebhrmz_rerr_FtotNorm'] = q_df.loc[tstamp, 'Ebq_rerr_FtotNorm']
        elif (tstamp not in q_df.index) and (tstamp in c_df.index):
            hrmz_df.loc[tstamp, 'CH4_eb_hrmnz'] = c_df.loc[tstamp, 'CH4_eb_c']
            hrmz_df.loc[tstamp, 'frac_eb_hrmnz'] = c_df.loc[tstamp, 'frac_eb_c']
            hrmz_df.loc[tstamp, 'CH4_diff_hrmnz'] = c_df.loc[tstamp, 'CH4_diff_c']
            hrmz_df.loc[tstamp, 'Ebhrmz_rerr_FtotNorm'] = c_df.loc[tstamp, 'Ebc_rerr_FtotNorm']  
        else:
            hrmz_df.loc[tstamp, 'CH4_eb_hrmnz'] = np.nan
            hrmz_df.loc[tstamp, 'frac_eb_hrmnz'] = np.nan
            hrmz_df.loc[tstamp, 'CH4_diff_hrmnz'] = np.nan
            hrmz_df.loc[tstamp, 'Ebhrmz_rerr_FtotNorm'] = np.nan

    hrmz_frame = full_df.join(hrmz_df, on='Timestamp')
    # get periods with no good value, drop them
    no_RefScal_idx = hrmz_frame[hrmz_frame['CH4_eb_hrmnz'].isna()].index
    hrmz_frame = hrmz_frame.drop(no_RefScal_idx, axis=0)    
    # drop negative ebullitive fluxes from frame
    hrmz_frame = hrmz_frame.drop(hrmz_frame.loc[hrmz_frame['CH4_eb_hrmnz'] <= 0].index)

    return(hrmz_frame)

#%% main function for making Wx scatterplots

def main_WxPlot(tstamp, hrmz_df):
    
    plt.ioff()
    
    # convert timestamp to format of file naming convention
    ## check to see if logging and local time zones are different; if so, read in files based on 'raw_file' datetimes
    if tz_logging != tz_local:
        day = hrmz_df.loc[tstamp, 'raw_file'][0:8]; datetime = hrmz_df.loc[tstamp, 'raw_file']
    else:
        day = tstamp.strftime('%Y%m%d'); datetime = tstamp.strftime('%Y%m%d_%H%M')
    
    # load data
    wvlt_df = pd.read_table(prog_loc + '/wvlet/%s/%s/wvlet-%s.dat' %(site_yr, day, datetime), sep='\s+',
                              names=['Index_wvlt', 'u', 'w', 'T', 'q', 'c', 'm'], delim_whitespace=False, 
                              skiprows=1)

    # get canopy height
    zc = np.float64(zc_df.loc[tstamp.floor('D')])
    # calculate displacement height as 0.66*canopy height
    d = 0.67*zc #[m]
    # calculate frequency for filtering
    u = np.mean(wvlt_df['u'])
    wvlt_df['j'] = (N_wvlt/wvlt_df['Index_wvlt'])*0.05
    wvlt_df['fn'] = (zm - d)/(wvlt_df['j']*u)
    
    # filter out low frequency components for partitioning
    wvlt_df_filt = wvlt_df[wvlt_df['fn'] > fn_LB]

    fig, ax = plt.subplots(nrows=1, ncols=3, sharey=False, figsize=(10,2.7))
    
    fig.suptitle('%s  %s\n' %(site_yr, tstamp.strftime('%d-%b %H:%M')) + 'F$_{CH_4, tot}$  = %.2f $\\mu$mol m$^{-2}$ s$^{-1}$  | LE = %.0f W m$^{-2}$  | F$_{CO_2}$ = %.1f $\\mu$mol m$^{-2}$ s$^{-1}$  | H = %.0f W m$^{-2}$  | $u_*$ = %.2f m s$^{-1}$ \n$\\sigma_m$ = %.4f mmol m$^{-3}$ | $\\bar{u}$ = %.2f m/s  | T$_{air}$ = %.1f $^{\\circ}$C  | RH = %.0f %%  |  wind direction = %.0f$^{\\circ}$  | $\\frac{z - d}{L}$ = %.2f'
                 %(hrmz_df.loc[tstamp, 'CH4_tot'], hrmz_df.loc[tstamp, LE_name], hrmz_df.loc[tstamp, Fco2_name], hrmz_df.loc[tstamp, sens_heat], hrmz_df.loc[tstamp, ust], hrmz_df.loc[tstamp, sigma_m], hrmz_df.loc[tstamp, ws_mean], hrmz_df.loc[tstamp, air_temp_C], hrmz_df.loc[tstamp, rel_hum], hrmz_df.loc[tstamp, WD], hrmz_df.loc[tstamp, MOST_stab]), size=8, y=1.18)
    
    ## Get axes limits and arrays for lines of best fit
    ymax01 = np.percentile(wvlt_df_filt['m'], 99.98); y_lim01 = (-1*ymax01, ymax01)
    ### mq
    xmax0 = np.percentile(wvlt_df_filt['q'], 99.97); x_lim0 = (-1*xmax0, xmax0)
    x_bf0 = np.array(x_lim0)
    y_bf0 = hrmz_df.loc[tstamp, 'slope_mq']*x_bf0
    ### mc
    xmax1 = np.percentile(wvlt_df_filt['c'], 99.97); x_lim1 = (-1*xmax1, xmax1)
    x_bf1 = np.array(x_lim1)
    y_bf1 = hrmz_df.loc[tstamp, 'slope_mc']*x_bf1    
    ### qc
    xmax2 = np.percentile(wvlt_df_filt['q'], 99.97); x_lim2 = (-1*xmax2, xmax2)
    ymax2 = np.percentile(wvlt_df_filt['c'], 99.98); y_lim2 = (-1*ymax2, ymax2)
    x_bf2 = np.array(x_lim2)
    y_bf2 = hrmz_df.loc[tstamp, 'slope_cq']*x_bf2  
    
    ## Create normalization for colorbar
    plot_norm = mpl.colors.LogNorm(vmin=fn_LB, vmax=200)
    ## color of ebullition threshold
    thresh_clr = 'r'
    
    ## mq plotting and annotations
    ax[0].scatter(wvlt_df_filt['q'], wvlt_df_filt['m'], s=10, c=wvlt_df_filt['fn'], cmap='plasma',
                  norm=plot_norm, edgecolor='k', linewidths=0.5, alpha=0.8)
    
    ### line of best fit 
    ax[0].plot(x_bf0, y_bf0, ls='-', c='dimgrey',
               lw=1.0)
    
    ### plot thresholds
    ax[0].plot(x_bf0, y_bf0 + 3*RMSDq_slope*hrmz_df.loc[tstamp, sigma_m], ls='--', c=thresh_clr[0], lw=1.0)
    ax[0].plot(x_bf0, y_bf0 - 3*RMSDq_slope*hrmz_df.loc[tstamp, sigma_m], ls='--', c=thresh_clr[0], lw=1.0)    
    title_0 = '$F_{eb,q}$/$F_{tot}$ = %.2f'%hrmz_df.loc[tstamp, 'frac_eb_q']
    qfilt = tstamp in qPart_df.index
    if qfilt == False:
        title_0 = title_0 + ' | $^{*}$flagged for bad quality'
    ax[0].set_title(title_0, size=8, c='k', loc='left')     
    
    ax[0].set_ylabel('$W_m$', size=8, labelpad=2); ax[0].set_xlabel('$W_q$', size=8)
    
    ## mc plotting and annotations
    ax[1].scatter(wvlt_df_filt['c'], wvlt_df_filt['m'], s=10, c=wvlt_df_filt['fn'], cmap='plasma',
                  norm=plot_norm, edgecolor='k', linewidths=0.5, alpha=0.8)
    
    ### line of best fit 
    ax[1].plot(x_bf1, y_bf1, ls='-', c='dimgrey',
               lw=1.0)
    
    ### plot thresholds
    ax[1].plot(x_bf1, y_bf1 + 3*RMSDc_slope*hrmz_df.loc[tstamp, sigma_m], ls='--', c=thresh_clr[0], lw=1.0)
    ax[1].plot(x_bf1, y_bf1 - 3*RMSDc_slope*hrmz_df.loc[tstamp, sigma_m], ls='--', c=thresh_clr[0], lw=1.0)    
    title_1 = '$F_{eb,c}$/$F_{tot}$ = %.2f'%hrmz_df.loc[tstamp, 'frac_eb_c']
    cfilt = tstamp in cPart_df.index
    if cfilt == False:
        title_1 = title_1 + ' | $^{*}$flagged for bad quality'
    ax[1].set_title(title_1, size=8, c='k', loc='left')
    
    ax[1].set_xlabel('$W_c$', size=8); ax[1].set_ylabel('$W_m$', size=8, labelpad=1)
    
    ## qc plotting and annotations
    ax[2].scatter(wvlt_df_filt['q'], wvlt_df_filt['c'], s=10, c=wvlt_df_filt['fn'], cmap='plasma',
                  norm=plot_norm, edgecolor='k', linewidths=0.5, alpha=0.8)
    
    ### line of best fit 
    ax[2].plot(x_bf2, y_bf2, ls='-', c='dimgrey',
               lw=1.0)
    
    title_2 = '$W_c$-$W_q$  R$^2$ = %.2f'%hrmz_df.loc[tstamp, 'R2_cq']
    ax[2].set_title(title_2, size=8, c='k', loc='center')
    
    ax[2].set_xlabel('$W_q$', size=8); ax[2].set_ylabel('$W_c$', size=8, labelpad=-3)
    
    ## set axes limits
    ax[0].set_xlim(x_lim0); ax[0].set_ylim(y_lim01)
    ax[1].set_xlim(x_lim1); ax[1].set_ylim(y_lim01)
    ax[2].set_xlim(x_lim2); ax[2].set_ylim(y_lim2)
    
    # scientifc notation for Wm axes
    ax[0].ticklabel_format(axis='y', style='scientific', scilimits=(0,0), useMathText=True); ax[0].yaxis.offsetText.set_fontsize(7); ax[0].yaxis.offsetText.set_position((-0.2,1.0))
    ax[1].ticklabel_format(axis='y', style='scientific', scilimits=(0,0), useMathText=True); ax[1].yaxis.offsetText.set_fontsize(7); ax[1].yaxis.offsetText.set_position((-0.2,1.0))
    
    
    ## modify plot elements for each subplot
    for p in np.arange(ax.size):
        ax.reshape(-1)[p].tick_params(which='major', axis='both', labelsize=7)
        ax.reshape(-1)[p].xaxis.set_minor_locator(AutoMinorLocator())
        ax.reshape(-1)[p].yaxis.set_minor_locator(AutoMinorLocator())
        ax.reshape(-1)[p].grid(which='major', axis='both', alpha=0.3)
        ax.reshape(-1)[p].set_axisbelow(True)
        # plot ref lines
        ax.reshape(-1)[p].axvline(0, linestyle='dotted', color='k', lw=1, zorder=0)
        ax.reshape(-1)[p].axhline(0, linestyle='dotted', color='k', lw=1, zorder=0)
    
    ## add color bar
    scm = plt.cm.ScalarMappable(cmap='plasma', norm=plot_norm)
    # scm.set_array([])
    fig.subplots_adjust(right=0.85)
    cbar_ax = fig.add_axes([0.855, 0.11, 0.01, 0.77])
    cbar_ax.tick_params(labelsize=7, length=5)
    cbar_ax.tick_params(length=3, which='minor')
    cbar_ax.minorticks_on()
    cbar = fig.colorbar(scm, cax=cbar_ax)
    cbar.set_label('f$_n$ (-)', size=8, labelpad=4)
    
    fig.subplots_adjust(wspace=0.28)
    
    # Make subfolder for each day if one doesn't already exist
    subf = '%s' %day
    subf_path = fig_path + '/' + subf
    if os.path.isdir(subf_path):
        pass
    else:
        os.mkdir(subf_path)
    
    plt.savefig(subf_path + '/%s.png' %datetime,  dpi=400, bbox_inches='tight')
    plt.clf()

    print('%s done' %datetime)
    return()

#%% get harmonized data frame for this site year

hrmz_df_siteyr = Hrmz_RefScal(full_HH_df, qPart_df, cPart_df)

#%% make a test plot first to ensure things look correct before iterating over all QC'd periods

main_WxPlot(pd.to_datetime('2015-07-05 13:00'), hrmz_df_siteyr)

#%%  now iterate through all QC'd periods

for n in np.arange(hrmz_df_siteyr.shape[0]):
    main_WxPlot(hrmz_df_siteyr.index[n], hrmz_df_siteyr)


