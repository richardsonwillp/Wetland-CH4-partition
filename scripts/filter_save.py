#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 22 14:39:31 2021

@author: willrichardson
"""
# add reference scalar flux flags, filter and save final copies of dataframes
# filter final frames for: CH4 flux magnitude, negative ebullition, and reference scalar flux magnitude

#%% import libraries
import os
import pandas as pd
import numpy as np
import datetime
from funcs import LE_name, RefScal_Type

#%% get relevant environment variables
site_yr = os.environ['site_yr']
run_ID = os.environ['run_ID']
fn_LB = os.environ['LB']
fn_UB = os.environ['UB']
tz_local = os.environ['tz_local']
tz_logging = os.environ['tz_logging']
NL_filt = bool(os.environ['NL_filt'])
FCH4_thresh = np.float64(os.environ['FCH4_min'])

# get working directory
base_path = os.getcwd()

#%% load master dataframe, intermediately filtered frames for each reference scalar, refernece scalar threshold values

# master, intermediately filtered frames
master_df = pd.read_csv(base_path + '/ref_data/Output/%s/%s/%s_Interm_RefDf_allHH_PostProcStats.csv' %(site_yr, run_ID, site_yr),
                        index_col=0, infer_datetime_format=True, parse_dates=True)


if NL_filt == True:
    qPart_filt = pd.read_csv(base_path + '/ref_data/Output/%s/%s_Interm_RefDf_qPart_filtHH_NLfilt.csv' %(site_yr, site_yr),
                          index_col='Timestamp', infer_datetime_format=True, parse_dates=True)
    
    cPart_filt = pd.read_csv(base_path + '/ref_data/Output/%s/%s_Interm_RefDf_cPart_filtHH_NLfilt.csv' %(site_yr, site_yr),
                          index_col='Timestamp', infer_datetime_format=True, parse_dates=True)

    q_fulloutput_fname = base_path + '/ref_data/Output/%s/%s/%s_FullOutput_qPart_filtHH_NLfilt.csv' %(site_yr, run_ID, site_yr)
    c_fulloutput_fname = base_path + '/ref_data/Output/%s/%s/%s_FullOutput_cPart_filtHH_NLfilt.csv' %(site_yr, run_ID, site_yr)

else:
    qPart_filt = pd.read_csv(base_path + '/ref_data/Output/%s/%s_Interm_RefDf_qPart_filtHH.csv' %(site_yr, site_yr),
                          index_col='Timestamp', infer_datetime_format=True, parse_dates=True)
    
    cPart_filt = pd.read_csv(base_path + '/ref_data/Output/%s/%s_Interm_RefDf_cPart_filtHH.csv' %(site_yr, site_yr),
                          index_col='Timestamp', infer_datetime_format=True, parse_dates=True) 
    q_fulloutput_fname = base_path + '/ref_data/Output/%s/%s/%s_FullOutput_qPart_filtHH.csv' %(site_yr, run_ID, site_yr)
    c_fulloutput_fname = base_path + '/ref_data/Output/%s/%s/%s_FullOutput_cPart_filtHH.csv' %(site_yr, run_ID, site_yr)

# reference scalar thresholds
brkpt_path = base_path + '/ref_data/Output/%s/%s/%s_FracEb_RefScalFlux_brkpts.csv' %(site_yr, run_ID, site_yr)
if os.path.exists(brkpt_path):
    brkpt_df = pd.read_csv(brkpt_path)

if RefScal_Type == 'user-supplied' or os.path.exists(brkpt_path) == False:
    LE_thresh = np.float64(os.environ['LE_thresh'])
    Fco2_thresh = np.float64(os.environ['Fco2_thresh'])
else:    
    LE_thresh = np.float64(brkpt_df['FracEbq_LE_brkpt']) + np.float64(brkpt_df['FracEbq_LE_brkpt_SE'])
    Fco2_thresh = np.float64(brkpt_df['FracEbc_Fco2_brkpt']) + np.float64(brkpt_df['FracEbc_Fco2_brkpt_SE'])

# add reference scalar flux flags (0 is good, 1 is flagged as being too low)
master_df['LE_flag'] = master_df[LE_name] < LE_thresh
master_df['Fco2_flag'] = master_df['co2_flux_mag'] < Fco2_thresh

#%% get intermediately filtered frames with all stats/flags added

qPart_df = master_df[master_df.index.isin(qPart_filt.index)]
cPart_df = master_df[master_df.index.isin(cPart_filt.index)]

#%% filter for conditions

qPart_df = qPart_df.loc[(qPart_df['CH4_tot'] >= FCH4_thresh) & (qPart_df['frac_eb_q'] >= 0.0) & (qPart_df['LE_flag'] == 0)]
cPart_df = cPart_df.loc[(cPart_df['CH4_tot'] >= FCH4_thresh) & (cPart_df['frac_eb_c'] >= 0.0) & (cPart_df['Fco2_flag'] == 0)]

# print final number of best quality periods
print('%i periods meeting all quality requirements for partitioning based on q' %qPart_df.shape[0])
print('%i periods meeting all quality requirements for partitioning based on c' %cPart_df.shape[0])

#%% check to see if data were logged in their local timezone; if not, convert time stamps
#   (logged timestamp is preserved in the 'raw_file' column of the dataframes)

if tz_logging != tz_local:
    master_df.index = master_df.index.tz_localize(tz_logging).tz_convert(tz_local).tz_localize(None)
    qPart_df.index = qPart_df.index.tz_localize(tz_logging).tz_convert(tz_local).tz_localize(None)
    cPart_df.index = cPart_df.index.tz_localize(tz_logging).tz_convert(tz_local).tz_localize(None)

#%% save frames

master_df.to_csv(base_path + '/ref_data/Output/%s/%s/%s_FullOutput_allHH.csv' %(site_yr, run_ID, site_yr), index_label='Timestamp')
qPart_df.to_csv(q_fulloutput_fname, index_label='Timestamp')
cPart_df.to_csv(c_fulloutput_fname, index_label='Timestamp')

#%% print current time to terminal and write to InfoFile for tracking purposes 

end_time = datetime.datetime.now().strftime('%Y-%m-%d %H:%M')
end_time_out = 'Run end time: %s'%end_time

print(end_time_out)

InfoFile_name = base_path + '/ref_data/Output/%s/%s_%s_PartRunInfo.txt' %(os.environ['site_yr'], os.environ['site_yr'], os.environ['run_ID'])

InfoFile = open(InfoFile_name, mode='a')
InfoFile.write('\n###################################\n' + end_time_out)
if RefScal_Type != 'user-supplied' and os.path.exists(brkpt_path) == False:
    InfoFile.write('\n ***Fitted reference scalar thresholds did not converge - user-supplied used instead')
InfoFile.close()

