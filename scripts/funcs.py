#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  9 15:28:30 2022

@author: willrichardson
"""

# this script houses all functions used in python-based processing steps of the partitioning program
# as well as variables holding strings with column names for variables of interest

#%% variables specifying column names for data selection, other characteristics of input files, etc.

# info on timestamps in raw data files and ancillary data files
raw_tstamp_beg = True # whether timestamps in raw data filenames refer to beginning or end of averaging period
HH_tstamp_beg = False # whether timestamps in ancillary dataset refer to beginning or end of averaging period
datetime_sep = True # whether date and time are separated in two columns or alread combined into one timestamp
tstamp_cols = ['date', 'time']
HH_output_header = 1
HH_output_skiprow = [2]
NaN_vals = [-9999, -999900]

# stationarity and integral turbulence characteristics tests columns
FCH4_SS = 'FCH4_SS_TEST'
FH2O_SS = 'FH2O_SS_TEST'
FCO2_SS = 'FC_SS_TEST'
w_ITC = 'W_ITC_TEST'
# other ancillary variables
sigma_m = 'st_dev(ch4)'
MOST_stab = '(z-d)/L'
WD = 'wind_dir'
ust = 'u*'
fetch_name = 'x_90%'
sigma_q = 'st_dev(h2o)'
sigma_c = 'st_dev(co2)'
wq_cov = 'cov(w/h2o)'
wc_cov = 'cov(w/co2)'
LE_name = 'LE'
Fco2_name = 'co2_flux'

# whether or not to convert units of some ancillary data in post-processing
anc_unit_conv = True
# optional names of environmental vars. if unit conversion desired (K to degC, Pa to kPa)
air_temp_K = 'air_temperature'
T_dew = 'Tdew'
P_atm = 'air_pressure'
VP = 'e'
VP_sat = 'es'
VPD_name = 'VPD'

# specify which reference scalar threshold to use ('user-supplied' overrides the the breakpoints calculated using broken-line regression)
RefScal_Type = 'user-supplied'

#%%  constants for flux-variance based non-local process filtering (DeBruin et al 1993)
c1 = 2.9; c2 = 20; bound_const = 1.7

#%% import libraries

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

#%% function for filtering full dataframe to based on quality requirements

def filt(df, ust_min, wd_min, wd_max, x90):
    
    qPart = df.copy(); cPart = df.copy(); m_frame = df.copy()
    
    # filter frames for (1) CH4 flux stationarity, (2) ref. scalar flux stationarity, and (3) ITC w test
    qPart_filt = qPart.loc[(qPart[FCH4_SS] <= 2) & (qPart[FH2O_SS] <= 2) & (qPart[w_ITC] <= 2)]
    cPart_filt = cPart.loc[(cPart[FCH4_SS] <= 2) & (cPart[FCO2_SS] <= 2) & (cPart[w_ITC] <= 2)]
    m_filt = m_frame.loc[(m_frame[FCH4_SS] <= 2) & (m_frame[w_ITC] <= 2)]
    
    print('%.i periods where CH4 and H2O fluxes meet stationarity and well-developed turbulence requirements' %qPart_filt.shape[0])
    print('%.i periods where CH4 and CO2 fluxes meet stationarity and well-developed turbulence requirements' %cPart_filt.shape[0])
    
    # filter for fetch requirements
    qPart_filt = qPart_filt.loc[(qPart_filt[WD] >= wd_min) & (qPart_filt[WD] <= wd_max)]
    cPart_filt = cPart_filt.loc[(cPart_filt[WD] >= wd_min) & (cPart_filt[WD] <= wd_max)]
    m_filt = m_filt.loc[(m_filt[WD] >= wd_min) & (m_filt[WD] <= wd_max)]
    
    if x90 != False:
        qPart_filt = qPart_filt.loc[qPart_filt[fetch_name] <= x90]
        cPart_filt = cPart_filt.loc[cPart_filt[fetch_name] <= x90]
        m_filt = m_filt.loc[m_filt[fetch_name] <= x90]

        print('%.i with best quality CH4 & H2O fluxes and from desired source' %qPart_filt.shape[0])
        print('%.i with best quality CH4 & CO2 fluxes and from desired source\n' %cPart_filt.shape[0])

    else:
        print('%.i with best quality CH4 & H2O fluxes and from desired source' %qPart_filt.shape[0])
        print('%.i with best quality CH4 & CO2 fluxes and from desired source\n' %cPart_filt.shape[0])

    if ust_min > 0:
        qPart_filt = qPart_filt.loc[(qPart_filt[ust] >= ust_min)]
        cPart_filt = cPart_filt.loc[(cPart_filt[ust] >= ust_min)]

        print('%.i with best quality CH4 and H2O fluxes, from desired source, and above u* threshold' %qPart_filt.shape[0])
        print('%.i with best quality CH4 and CO2 fluxes, from desired source, and above u* threshold\n' %cPart_filt.shape[0])

    return(qPart_filt, cPart_filt, m_filt) 

#%% function for getting the univ. func and its boundaries for a given stability

def get_UnivFunc(df):
    univ_out = np.empty(df.shape[0], dtype=np.float64)
    for i in np.arange(df.shape[0]):
        if df.loc[df.index[i], MOST_stab] > 0:
            univ_out[i] = 2
        else:
            univ_out[i] = c1*(1 - c2*df.loc[df.index[i], MOST_stab])**(-1/3)
    
    univ_out_low = univ_out/bound_const
    univ_out_high = univ_out*bound_const

    univ_vals = pd.DataFrame(np.stack([univ_out, univ_out_low, univ_out_high], axis=1), index=df.index, columns=['MOST_phi', 'MOST_phi_lower', 'MOST_phi_upper'])
    return(univ_vals)

#%% Helper functions for doing regressions

#   Does linear regression (forced through origin), returns results and points defining line of best fit
def LR_xy_NC(df, x_var, y_var):
    
    reg = df[[y_var, x_var]].dropna() 
    
    # rename variables that would throw off formulas
    if '(' in x_var:
        x_var_new = x_var.replace('(', '_'); x_var_new = x_var_new.replace(')', '')
        reg = reg.rename(columns = {x_var: x_var_new})
        x_var = x_var_new
        
    if '(' in y_var:
        y_var_new = y_var.replace('(', '_'); y_var_new = y_var_new.replace(')', '')
        reg = reg.rename(columns = {y_var: y_var_new})        
        y_var = y_var_new
    
    mod = smf.ols(formula='%s ~ 0 + %s' %(y_var, x_var), data=reg)
    res = mod.fit()
    x_plot = np.linspace(np.min(reg[x_var]), np.max(reg[x_var]), reg.shape[0])
    y_plot = res.params[0]*x_plot
    bf = np.stack((x_plot, y_plot), axis=1)
    
    return(res, bf)

# Same as above but doing quantile regression across the user-defined range of quantiles
def QuantReg_xy_NC(df, x_var, y_var, percentile_arr):
    
    reg = df[[y_var, x_var]].dropna()
    
    # rename variables that would throw off formulas
    if '(' in x_var:
        x_var_new = x_var.replace('(', '_'); x_var_new = x_var_new.replace(')', '')
        reg = reg.rename(columns = {x_var: x_var_new})
        x_var = x_var_new
        
    if '(' in y_var:
        y_var_new = y_var.replace('(', '_'); y_var_new = y_var_new.replace(')', '')
        reg = reg.rename(columns = {y_var: y_var_new})        
        y_var = y_var_new
    
    mod = smf.quantreg(formula='%s ~ 0 + %s' %(y_var, x_var), data=reg)
    res_all = {}
    bf_all = {}
    for i in np.arange(len(percentile_arr)):
        
        res = mod.fit(q=percentile_arr[i])
        out_key = np.around(percentile_arr[i], 2)
        # print(out_key)
        # print(len(str(out_key)))
        if len(str(out_key)) == 3:
            # print(True)
            out_key = str(out_key).ljust(4, '0')
        else:
            out_key = str(out_key)
        res_all[out_key] = res
        
        x_plot = np.linspace(np.min(reg[x_var]), np.max(reg[x_var]), reg.shape[0])
        y_plot = res.params[0]*x_plot
        bf = np.stack((x_plot, y_plot), axis=1)
        bf_all[out_key] = bf

    return(res_all, bf_all)

#%% helper functions for plotting

# QR lines of best fit
def ax_Add_bf_lines(ax, bf_dict_list, res_dict_list, color_list):
        
    lwidth = 0.9; lstyle = 'solid'
    for i in np.arange(ax.size):
        bf_dict = bf_dict_list[i]; res_dict = res_dict_list[i]
        c_index = 0
        for percentile in bf_dict.keys():
            ax.reshape(-1)[i].plot(bf_dict[percentile][:,0], bf_dict[percentile][:,1], c=color_list[c_index], lw=lwidth, ls=lstyle,
                       label='Q = %s: y = %.3fx' %(percentile, res_dict[percentile].params[0]))
            c_index = c_index + 1
    return(ax)

# adding another line of best fit for non-ebulltive periods 
def ax_AddLine(ax, slopes, line_label, linecolor):
    clr = linecolor; lstyle = 'dashed'; lwidth = 0.9
    for i in np.arange(ax.size):
        x_plot = np.array([0, ax[i].get_xlim()[1]])
        ax.reshape(-1)[i].plot(x_plot, slopes[i]*x_plot, ls=lstyle, c=clr, lw=lwidth, label='%s: y = %.3fx' %(line_label, slopes[i]))

    return(ax)

# plotting non-ebullitive periods on top in a different color
def ax_AddPoints_NoC(ax, dfs, y_var, x_var, EdgeColor):
    lw=0.6; ms = 8; face = 'none'
    
    if type(x_var) == str:
        for i in np.arange(ax.size):
            ax.reshape(-1)[i].scatter(dfs[i][x_var], dfs[i][y_var[i]], s=ms, facecolors=face, edgecolors=EdgeColor,
                          linewidths=lw, zorder=6)
    else:
        for i in np.arange(ax.size):
            ax.reshape(-1)[i].scatter(dfs[i][x_var[i]], dfs[i][y_var[i]], s=ms, facecolors=face, edgecolors=EdgeColor,
                          linewidths=lw, zorder=6)    
    return(ax)

#%% functions for sorting files from partitioning output, reading them in and making a dataframe

# sort files that end in a string with YYYYMMDD chronologically
def sort_ByDate_DMY(filelist):
    filelist = list(filter(lambda x: x[0] != '.', filelist))
    filelist = sorted(filelist, key = lambda x: int(x[-6:-4]))
    filelist = sorted(filelist, key = lambda x: int(x[-8:-6]))
    filelist = sorted(filelist, key = lambda x: int(x[-12:-8]))
    
    return(filelist)

# read flux and random error data in, return dataframe and raw file names
def read_conc(flux_filelist, rerr_filelist):
    
    flux_headers = ['Timestamp', 'CH4_tot', 'CH4_eb_T', 'CH4_eb_q', 'CH4_eb_c', 'LFC_CH4'] # first col is string, all others in units of mmol*m^-2*s^-1
    rerr_headers = ['Timestamp', 'CH4_diffT_err', 'CH4_ebT_err', 'CH4_diffq_err', 'CH4_ebq_err','CH4_diffc_err', 'CH4_ebc_err'] # first col is string, all others in units of mmol*m^-2*s^-1
    comb_headers = flux_headers[1:] + rerr_headers[1:]
    
    df = pd.DataFrame(columns=comb_headers)

    for f in np.arange(len(flux_filelist)):
        
        # Read flux estimates, make timestamps the index
        flux_frame = pd.read_table(flux_filelist[f], sep='\s+', names=flux_headers, skiprows=1, delim_whitespace=False)
        flux_frame['Timestamp'] = pd.to_datetime(flux_frame['Timestamp'], format='wvlet-%Y%m%d_%H%M.dat')
        flux_frame = flux_frame.set_index('Timestamp')

        # Read random error estimates, make timestamps the index
        rerr_frame = pd.read_table(rerr_filelist[f], sep='\s+', names=rerr_headers, skiprows=1, delim_whitespace=False)
        rerr_frame['Timestamp'] = pd.to_datetime(rerr_frame['Timestamp'], format='wvlet-%Y%m%d_%H%M.dat')
        rerr_frame = rerr_frame.set_index('Timestamp')
        
        # Merge frames
        comb_frame = flux_frame.merge(rerr_frame, on='Timestamp')
        
        # append to full df
        df = df.append(comb_frame)
        df.index.name = 'Timestamp'
        # check to see if raw data filenames refer to beginning or end of averaging period
        if raw_tstamp_beg == False:
            df.index = df.index - pd.Timedelta(30, 'm')
    
    raw_files = pd.Series(data=df.index.strftime('%Y%m%d_%H%M'), index=df.index, dtype=str, name='raw_file')
            
    return(raw_files, df)

#%% function for getting bin averages

def getBinAvg(dfs, x_vars, y_vars, nBins):
    
    out_df_all = {}
    
    for i in np.arange(len(x_vars)):
        x_var = x_vars[i]; y_var = y_vars[i]
        df = dfs[i]
        df_sort = df[[x_var, y_var]]
        df_sort = df_sort.sort_values(x_var)
        
        df_split = np.array_split(df_sort, nBins)
        
        x_mean = []; x_se = []; y_mean = []; y_se = []
        for j in np.arange(len(df_split)):
            x_mean.append(np.nanmean(df_split[j][x_var]))
            x_se.append(np.nanstd(df_split[j][x_var]))
            y_mean.append(np.nanmean(df_split[j][y_var]))
            y_se.append(np.nanstd(df_split[j][y_var]))

        x_mean = np.asarray(x_mean); x_se = np.asarray(x_se); y_mean = np.asarray(y_mean); y_se = np.asarray(y_se)
        
        bin_stats = pd.DataFrame(np.stack([x_mean, x_se, y_mean, y_se], axis=1), columns=[x_var + '_mean', x_var + '_sd', y_var + '_mean', y_var + '_sd'])
        
        out_df_all[str(i)] = bin_stats
        
    return(out_df_all)
