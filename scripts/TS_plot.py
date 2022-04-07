#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  4 13:22:30 2021

@author: willrichardson
"""
# this script makes times series plots of the high frequency data for all observation periods
## NaNs in the time series of each variable are plotted as red circles; the total
## for each variable is marked in a legend on each subplot
    
#%% user defined variables

# time zone info
logging_tz = 'America/Chicago'
local_tz = 'America/Chicago'
# boolean for whether timestmap in filenames refer to beginning or end of averaging period
tstamp_beg = False

# filepaths
## NOTE: these high frequency data should be those prior to NaN removal
HiFreq_data_dir = '/Volumes/Backup-Will/DATA_2015_2016/Way3/Summer_2015/EddyPro/Output/eddypro_raw_datasets/level_6'
plot_dir = '/Volumes/Backup-Will/DATA_2015_2016/Way3/Summer_2015/test_TS'

# specifics on raw data files
## file delimiter
delim = ''
## number of header lines in files
header_rows = 10
## specify where the timestamp information is in the filename (counting starts at 0)
year_bnd = [0,4]
month_bnd = [4,6]
day_bnd = [6,8]
hour_bnd = [9,11]
minute_bnd = [11,13]
## value representing NaNs
NaN_val = -9999.0
## a list of the variables in the order their columns appear in the raw data; use the following variable options
### 3D wind vector ('u', 'v', 'w'), sonic temp ('Ts'), co2 ('c'), h2o ('q'), ch4 ('m'), air temp ('Ta'), air pressure ('P')
var_cols = ['u', 'v', 'w', 'Ts', 'c', 'q', 'm', 'Ta', 'P']
## specify units of the respective columns (use LaTeX notation for pretty axes labels)
col_units = ['m s$^{-1}$', 'm s$^{-1}$', 'm s$^{-1}$', 'K', 'mmol m$^{-3}$',
             'mmol m$^{-3}$', 'mmol m$^{-3}$', 'K', 'Pa']

#%% import libraries
# import sys
# sys.path.insert(0, '/Users/willrichardson/opt/anaconda3/lib/python3.8/site-packages')
import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

#%% other variables used by program
# array of frequencies for plotting in minutes
freqs = np.array([10, 20]) # array of possible frequencies
# plotting variables
pres_lab_fs = 10.5
pres_ticklab_fs = 8
pres_labpad = 10
# font sizes, line thicknesses, etc.
fs = 12 # points
lab_fs = 8
lw = 0.4 # points

#%% bottom cell: main routine

plt.ioff()

# list files in directory
files = os.listdir(HiFreq_data_dir)
# get rid of any weird temp files
files = list(filter(lambda x: x[0] != '.', files))
# sort files chronologically
files_sorted = sorted(files, key = lambda x: int(x[minute_bnd[0]:minute_bnd[1]]))
files_sorted = sorted(files_sorted, key = lambda x: int(x[hour_bnd[0]:hour_bnd[1]]))
files_sorted = sorted(files_sorted, key = lambda x: int(x[day_bnd[0]:day_bnd[1]]))
files_sorted = sorted(files_sorted, key = lambda x: int(x[month_bnd[0]:month_bnd[1]]))
files_sorted = sorted(files_sorted, key = lambda x: int(x[year_bnd[0]:year_bnd[1]]))

# Counter variable for tracking progress of code
n_file_iter = 1
n_files = len(files_sorted)

# iterate for all files!
for fname in files_sorted:
    # load data
    data = np.genfromtxt(HiFreq_data_dir + '/' + fname, delimiter=delim, skip_header=header_rows,
                         missing_values=NaN_val, filling_values=np.nan)
    # genfromtxt doesn't always catch NaNs, replace them manually
    data[data == NaN_val] = np.nan
    # get timestamp, adjust if not referring to beginning of averaging period
    tstamp_str = fname[year_bnd[0]:year_bnd[1]] + fname[month_bnd[0]:month_bnd[1]] + fname[day_bnd[0]:day_bnd[1]] + '-' + fname[hour_bnd[0]:hour_bnd[1]] + fname[minute_bnd[0]:minute_bnd[1]]
    datetime = pd.to_datetime(tstamp_str, format='%Y%m%d-%H%M')
    if tstamp_beg == False:
        datetime = datetime - pd.Timedelta(30, 'm')
    
    tstamp_str_out = datetime.strftime('%Y%m%d_%H%M')
    
    print('Plotting period %i of %i... %s' %(n_file_iter, n_files, tstamp_str_out))

    # Count NaNs; get extra plotting variables
    NaN_count = np.empty(data.shape[1]); means = np.empty(data.shape[1])
    NaN_locs = {}
    for col in np.arange(data.shape[1]):
        NaN_count[col] = np.sum(np.isnan(data[:,col]))
        means[col] = np.nanmean(data[:,col])
        NaN_loc = np.argwhere(np.isnan(data[:,col]))
        NaN_locs.update({col: NaN_loc})
    
    # make a variable for easier time series plotting
    T = data.shape[0]/(60*freqs) # length of data period if each possible frequency were right
    freq = freqs[np.argmin(np.absolute(30 - T))] # set frequency as value which yields a length closest to half-hour
    minutes = np.arange(0, data.shape[0], 1)/(freq*60)
    
    # make plot
    fig, ax = plt.subplots(nrows = len(var_cols), ncols=1, sharex=True)
    if logging_tz != local_tz:
        true_tstamp = datetime.tz_localize(logging_tz).tz_convert(local_tz).tz_localize(None).strftime('%Y-%m-%d %H:%M')
        fig.suptitle('%s (%s LT)' %(datetime.strftime('%Y-%m-%d %H:%M'), true_tstamp),
                     y=0.9, size=pres_lab_fs+1)
    else:
        fig.suptitle(datetime.strftime('%Y-%m-%d %H:%M'),
                 y=0.9, size=pres_lab_fs+1)
    
    ## Subplots 
    for i in np.arange(len(var_cols)):
        
        # plot data
        ax[i].plot(minutes, data[:,i], '-', c='k', linewidth=lw, alpha=0.9, zorder=2)

        # plot NaNs
        if NaN_count[i] > 0:
            flat_locs = np.squeeze(NaN_locs.get(i), axis=-1)
            means_plotting = np.ones(len(flat_locs))*means[i]
            ax[i].plot(minutes[flat_locs], means_plotting, 'ro', ms=6, fillstyle='none', mew=0.5, zorder=2, label=r'# NaN = %i' %NaN_count[i])
            ax[i].legend(loc=0, fontsize=fs-5, handletextpad=0.2, framealpha=0.3)
        
        # lables and formatting
        ylab = var_cols[i] + '(%s)' %col_units[i]
        ax[i].set_ylabel(ylab, size=pres_lab_fs, labelpad=pres_labpad)
        ax[i].tick_params('y', labelsize=pres_ticklab_fs)
        ax[i].ticklabel_format(axis='y', style='plain', useOffset=False)
        # put a major xtick every 5 minutes and a minor xtick at each minute
        ax[i].xaxis.set_major_locator(MultipleLocator(5))
        ax[i].xaxis.set_minor_locator(MultipleLocator(1))
        ax[i].grid(which='major', axis='both', alpha=0.5)
        ax[i].grid(which='minor', axis='both', alpha=0.2)
        
        # add x label to bottom supblot
        if i == (len(var_cols) - 1):
            ax[i].set_xlabel('Time (min)', size=pres_lab_fs)
            ax[i].tick_params('x', labelsize=pres_ticklab_fs)
            ax[i].set_xlim((-0.5,30.5))
    
    # Manually set size
    fig.set_size_inches(7, 12)
    
    # Make subfolder for each day if one doesn't already exist
    subf = tstamp_str_out[0:8] 
    subf_path = plot_dir + '/' + subf
    if os.path.isdir(subf_path):
        pass
    else:
        os.mkdir(subf_path)
    
    # Save fig; time of day is the filename 
    plt.savefig(subf_path + '/%s.png' %tstamp_str_out, dpi=300, bbox_inches='tight')
    plt.clf()
    plt.close(fig)
    
    # keep track of which file you're on
    n_file_iter = n_file_iter + 1
