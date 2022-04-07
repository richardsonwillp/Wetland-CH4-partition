#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  9 14:13:23 2021

@author: willrichardson
"""
# This script calculates other relevant stats and quantities of interest for each half-hour period
# includes 30-minute quantites and spectral quantites in each period

# leaving RMSDs as local variables; pass them into this script

#%% import libraries
import sys
sys.path.insert(0, '/Users/willrichardson/opt/anaconda3/lib/python3.8/site-packages')
import os
import glob
import pandas as pd
import numpy as np
from funcs import sort_ByDate_DMY, read_conc
from funcs import anc_unit_conv, Fco2_name, sigma_m

#%% get current working directory, relevant environment variables
args = sys.argv

RMSDq_slope = np.float64(args[1])
RMSDc_slope = np.float64(args[2])

base_path = os.getcwd()
site_yr = os.environ['site_yr']
run_ID = os.environ['run_ID']
fn_LB = os.environ['LB']
fn_UB = os.environ['UB']

DT = np.float64(os.environ['DT'])
nj = int(np.rint(np.log2((30*60)/DT)))
zm = np.float64(os.environ['Zm'])
N_wvlt = 2**nj

# make output file name
output_fname = base_path + '/ref_data/Output/%s/%s/%s_Interm_RefDf_allHH_PostProcStats.csv'%(site_yr, run_ID, site_yr)

#%% load reference dfs, outputs from partitioning program, other ancillary info

df_ref = pd.read_csv(base_path + '/ref_data/Output/%s/%s_Interm_RefDf_allHH.csv' %(site_yr, site_yr),
                      index_col='Timestamp', infer_datetime_format=True, parse_dates=True)

RLM_df = pd.read_csv(base_path + '/ref_data/Output/%s/%s_RLMfitRMSDs_fnLB_%s_fnUB_%s.csv' %(site_yr, site_yr, fn_LB, fn_UB),
                     index_col=0, infer_datetime_format=True, parse_dates=True)

zc_df = pd.read_csv(os.environ['zc_file'], sep='\t', index_col=0, names=['Timestamp', 'zc'], parse_dates=True, infer_datetime_format=True)

# partitioning program outputs
## load and sort partitioned fluxes and random error estimates
flux = glob.glob(base_path + '/flux/%s/flux*.txt'%site_yr)
rerr = glob.glob(base_path + '/flux/%s/rerror*.txt'%site_yr)
flux = sort_ByDate_DMY(flux); rerr = sort_ByDate_DMY(rerr) 
raw_files, part_df = read_conc(flux, rerr)
## partitioning program yields fluxes in mmol m-2 s-1; convert to umol m-2 s-1
part_df = part_df*1000

#%% join into a unified dataframe (use 'inner' join so as to only keep rows where raw data made it all the way to the partitioning stage)

df_master = df_ref.join([part_df, RLM_df], how='inner')

#%% 30 minute stats

# add raw filenames as a column in the data frame
df_master['raw_file'] = raw_files
# fraction of total CH4 flux that is ebullition
df_master['frac_eb_q'] = df_master['CH4_eb_q']/df_master['CH4_tot']
df_master['frac_eb_c'] = df_master['CH4_eb_c']/df_master['CH4_tot']
df_master['frac_eb_T'] = df_master['CH4_eb_T']/df_master['CH4_tot'] 
# diffusive fluxes
df_master['CH4_diff_q'] = df_master['CH4_tot'] - df_master['CH4_eb_q']
df_master['CH4_diff_c'] = df_master['CH4_tot'] - df_master['CH4_eb_c']
df_master['CH4_diff_T'] = df_master['CH4_tot'] - df_master['CH4_eb_T']

# Some unit sonversions if desired
if anc_unit_conv == True:
    from funcs import air_temp_K, T_dew, P_atm, VP, VP_sat, VPD_name
    df_master['T_air_C'] = df_master[air_temp_K] - 273.15; df_master[T_dew] = df_master[T_dew] - 273.15 #[K to C]
    df_master[VP] = df_master[VP]/1000; df_master[VP_sat] = df_master[VP_sat]/1000 #[Pa to kPa]
    df_master[VPD_name] = df_master[VPD_name]/1000; df_master['P'] = df_master[P_atm]/1000 # [Pa to kPa]
    df_master.drop([P_atm, air_temp_K], axis=1, inplace=True)

# co2 flux magnitude (for reference scalar thresholding)
df_master['co2_flux_mag'] = np.absolute(df_master[Fco2_name])

# normalized random error stats
df_master['Ebq_rerr_FebNorm'] = df_master['CH4_ebq_err']/df_master['CH4_eb_q']
df_master['Ebq_rerr_FtotNorm'] = df_master['CH4_ebq_err']/df_master['CH4_tot']    
df_master['Diffq_rerr_FdiffNorm'] = df_master['CH4_diffq_err']/df_master['CH4_diff_q']
df_master['Diffq_rerr_FtotNorm'] = df_master['CH4_diffq_err']/df_master['CH4_tot']
df_master['Ebc_rerr_FebNorm'] = df_master['CH4_ebc_err']/df_master['CH4_eb_c']
df_master['Ebc_rerr_FtotNorm'] = df_master['CH4_ebc_err']/df_master['CH4_tot']    
df_master['Diffc_rerr_FdiffNorm'] = df_master['CH4_diffc_err']/df_master['CH4_diff_c']
df_master['Diffc_rerr_FtotNorm'] = df_master['CH4_diffc_err']/df_master['CH4_tot']

#%% function for spectral stats on each period

def SpectralStats(tstamp):
    
    # convert timestamp to format of file naming convention
    day = tstamp.strftime('%Y%m%d'); datetime = tstamp.strftime('%Y%m%d_%H%M')
    # load data; first row is coarse-grained mean, skip it
    wvlt_df = pd.read_table(base_path + '/wvlet/%s/'%site_yr + day + '/' + 'wvlet-' + datetime + '.dat', sep='\s+',
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
    wvlt_df_filt = wvlt_df[wvlt_df['fn'] > np.float64(fn_LB)]
    # add ebullition flags so that stats on diffusive/ebullitive fluxes can be calculated
    wvlt_df_filt.loc[:, 'thresh_q_upp'] = df_master.loc[tstamp, 'slope_mq']*wvlt_df_filt.loc[:, 'q'] + (3*(RMSDq_slope*df_master.loc[tstamp, sigma_m]))
    wvlt_df_filt.loc[:, 'thresh_q_low'] = df_master.loc[tstamp, 'slope_mq']*wvlt_df_filt.loc[:, 'q'] - (3*(RMSDq_slope*df_master.loc[tstamp, sigma_m]))
    wvlt_df_filt.loc[:, 'thresh_c_upp'] = df_master.loc[tstamp, 'slope_mc']*wvlt_df_filt.loc[:, 'c'] + (3*(RMSDc_slope*df_master.loc[tstamp, sigma_m]))
    wvlt_df_filt.loc[:, 'thresh_c_low'] = df_master.loc[tstamp, 'slope_mc']*wvlt_df_filt.loc[:, 'c'] - (3*(RMSDc_slope*df_master.loc[tstamp, sigma_m]))

    wvlt_df_filt.loc[:, 'Ebq'] = (wvlt_df_filt.loc[:, 'm'] > wvlt_df_filt.loc[:, 'thresh_q_upp']) | (wvlt_df_filt.loc[:, 'm'] < wvlt_df_filt.loc[:, 'thresh_q_low'])
    wvlt_df_filt.loc[:, 'Ebc'] = (wvlt_df_filt.loc[:, 'm'] > wvlt_df_filt.loc[:, 'thresh_c_upp']) | (wvlt_df_filt.loc[:, 'm'] < wvlt_df_filt.loc[:, 'thresh_c_low'])

    # do spectral stats
    ## result objects
    ### master output for funtion
    out = pd.DataFrame(index=[tstamp])
    ### frames for scale-wise stats
    #### frequency values
    fn_res = pd.DataFrame(index=[tstamp])
    #### variances
    mvar_res = pd.DataFrame(index=[tstamp]); qvar_res = pd.DataFrame(index=[tstamp])
    cvar_res = pd.DataFrame(index=[tstamp]); Tvar_res = pd.DataFrame(index=[tstamp]); wvar_res = pd.DataFrame(index=[tstamp])
    #### ebullitive covariances
    Ebq_cov_res = pd.DataFrame(index=[tstamp]); Ebc_cov_res = pd.DataFrame(index=[tstamp]) 
    #### covariances
    wmcov_res = pd.DataFrame(index=[tstamp]); wqcov_res = pd.DataFrame(index=[tstamp]); wccov_res = pd.DataFrame(index=[tstamp]); wTcov_res = pd.DataFrame(index=[tstamp])
    mqcoh_res = pd.DataFrame(index=[tstamp]); mccoh_res = pd.DataFrame(index=[tstamp]); qccoh_res = pd.DataFrame(index=[tstamp]);
    ## variables for 30-min. values calculated from scale-wise stats
    m_var = 0; q_var = 0; c_var = 0; T_var = 0; w_var = 0;
    m_ebq_var = 0; m_ebc_var = 0
    wm_tot_cov = 0; wm_ebq_cov = 0; wm_ebc_cov = 0
    wq_tot_cov = 0; wc_tot_cov = 0; wT_tot_cov = 0
    
    if wvlt_df.isnull().values.any() == False and wvlt_df.shape[0] > 30000:
    
        for i in np.arange(nj):
            # select data at this scale
            scale_index = 2**i        
            # all scales of transform
            scale_data = wvlt_df.loc[wvlt_df['Index_wvlt'] == scale_index]
            # dataframe filtered for coefficients to be included in partitioning
            scale_data_filt = wvlt_df_filt.loc[wvlt_df_filt['Index_wvlt'] == scale_index]
                
            # get labelling string of this scale for results; make scale-wise variable labels
            scal_lab = 'j' + str(15-i) + '_ts' + str(np.float64(scale_data['j'].iloc[0]))
            fn_scal_lab = 'fn_' + scal_lab
            fn_scale = np.float64(np.unique(scale_data['fn'])); fn_res[fn_scal_lab] = fn_scale
            
            ## variances
            mlab = 'm_NormVar_' + scal_lab; qlab = 'q_NormVar_' + scal_lab; clab = 'c_NormVar_' + scal_lab;
            Tlab = 'Ts_NormVar_' + scal_lab; wlab = 'w_NormVar_' + scal_lab; 
            ## covariances
            wm_lab = 'wm_NormCov_' + scal_lab; wq_lab = 'wq_NormCov_' + scal_lab; wc_lab = 'wc_NormCov_' + scal_lab;
            wT_lab = 'wT_NormCov_' + scal_lab;
            Ebq_lab = 'Ebq_NormCov_' + scal_lab; Ebc_lab = 'Ebc_NormCov_' + scal_lab
            ## coherences
            mqcoh_lab1 = 'mq_coh1_' + scal_lab; mccoh_lab1 = 'mc_coh1_' + scal_lab; qccoh_lab1 = 'qc_coh1_' + scal_lab
    
            # variances and covariances (variance calculated in this iterative process, std for m calculated from square root of global variance at the end)
            mvar_scale = np.sum(scale_data['m']**2)/N_wvlt; mvar_res.loc[tstamp, mlab] = mvar_scale; m_var = m_var + mvar_scale
            qvar_scale = np.sum(scale_data['q']**2)/N_wvlt; qvar_res.loc[tstamp, qlab] = qvar_scale; q_var = q_var + qvar_scale
            cvar_scale = np.sum(scale_data['c']**2)/N_wvlt; cvar_res.loc[tstamp, clab] = cvar_scale; c_var = c_var + cvar_scale
            Tvar_scale = np.sum(scale_data['T']**2)/N_wvlt; Tvar_res.loc[tstamp, Tlab] = Tvar_scale; T_var = T_var + Tvar_scale
            wvar_scale = np.sum(scale_data['w']**2)/N_wvlt; wvar_res.loc[tstamp, wlab] = wvar_scale; w_var = w_var + wvar_scale
            
            wmcov_scale = np.sum(scale_data['w']*scale_data['m'])/N_wvlt; wmcov_res.loc[tstamp, wm_lab] = wmcov_scale; wm_tot_cov = wm_tot_cov + wmcov_scale
            wqcov_scale = np.sum(scale_data['w']*scale_data['q'])/N_wvlt; wqcov_res.loc[tstamp, wq_lab] = wqcov_scale; wq_tot_cov = wq_tot_cov + wqcov_scale
            wccov_scale = np.sum(scale_data['w']*scale_data['c'])/N_wvlt; wccov_res.loc[tstamp, wc_lab] = wccov_scale; wc_tot_cov = wc_tot_cov + wccov_scale
            wTcov_scale = np.sum(scale_data['w']*scale_data['T'])/N_wvlt; wTcov_res.loc[tstamp, wT_lab] = wTcov_scale; wT_tot_cov = wT_tot_cov + wTcov_scale
            
            if scale_data_filt.empty != True:
                
                # variances
                m_ebq_var = m_ebq_var + np.sum(scale_data_filt['Ebq']*scale_data_filt['m']**2)/N_wvlt
                m_ebc_var = m_ebc_var + np.sum(scale_data_filt['Ebc']*scale_data_filt['m']**2)/N_wvlt
                
                # covariances
                Ebq_cov_scale = np.sum(scale_data_filt['Ebq']*scale_data_filt['w']*scale_data_filt['m'])/N_wvlt
                Ebq_cov_res[Ebq_lab] = Ebq_cov_scale
                wm_ebq_cov = wm_ebq_cov + Ebq_cov_scale
                Ebc_cov_scale = np.sum(scale_data_filt['Ebc']*scale_data_filt['w']*scale_data_filt['m'])/N_wvlt
                Ebc_cov_res[Ebc_lab] = Ebc_cov_scale
                wm_ebc_cov = wm_ebc_cov + Ebc_cov_scale
                    
            # calculate coherences (must double check to make sure this scale has 3 or more data points)
            if scale_data.shape[0] >= 3:
                # (1) just using an out of box correlation coefficient calculation
                mqcoh_res[mqcoh_lab1] = np.corrcoef(scale_data['q'], scale_data['m'])[0][1]
                mccoh_res[mccoh_lab1] = np.corrcoef(scale_data['c'], scale_data['m'])[0][1]
                qccoh_res[qccoh_lab1] = np.corrcoef(scale_data['q'], scale_data['c'])[0][1]
    
    # convert variances to stdevs; normalize scale-wise variances
    m_std = np.sqrt(m_var); m_ebq_std = np.sqrt(m_ebq_var); m_diffq_std = m_std -  m_ebq_std
    m_ebc_std = np.sqrt(m_ebc_var); m_diffc_std = m_std - m_ebc_std
    mvar_res = mvar_res/m_var; qvar_res = qvar_res/q_var; cvar_res = cvar_res/c_var; Tvar_res = Tvar_res/T_var; wvar_res = wvar_res/w_var;
    
    # convert w-m covariances from mmol m-2 s-1 to umol m-2 s-1; normalize scale-wise covariances
    wmcov_res = wmcov_res/wm_tot_cov; Ebq_cov_res = Ebq_cov_res/wm_ebq_cov; Ebc_cov_res = Ebc_cov_res/wm_ebc_cov
    wm_tot_cov = wm_tot_cov*1000; wm_ebq_cov = wm_ebq_cov*1000; wm_diffq_cov = wm_tot_cov - wm_ebq_cov
    wm_ebc_cov = wm_ebc_cov*1000; wm_diffc_cov = wm_tot_cov - wm_ebc_cov
    wqcov_res = wqcov_res/wq_tot_cov; wccov_res = wccov_res/wc_tot_cov; wTcov_res = wTcov_res/wT_tot_cov
        
    # save all stats to output frame
    ## stanalone values not already in a preliminary dataframe
    cols_write = ['var_m', 'stdev_m', 'stdev_m_ebq', 'stdev_m_diffq', 'stdev_m_ebc', 'stdev_m_diffc', 'cov_wm', 'cov_wm_ebq',
                'cov_wm_diffq', 'cov_wm_ebc','cov_wm_diffc', 'var_q', 'var_c', 'var_T', 'var_w',
                'cov_wq', 'cov_wc', 'cov_wT']
    output_arr = np.array([m_var, m_std, m_ebq_std, m_diffq_std, m_ebc_std, m_diffc_std, wm_tot_cov, wm_ebq_cov, wm_diffq_cov,
                           wm_ebc_cov, wm_diffc_cov, q_var, c_var, T_var, w_var, wq_tot_cov, wc_tot_cov, wT_tot_cov])
    output_arr = np.expand_dims(output_arr, 0)
    
    out_interm = pd.DataFrame(output_arr, index=[tstamp], columns=cols_write)
    out = out.join([fn_res, out_interm, mvar_res, qvar_res, cvar_res, Tvar_res, wvar_res, wmcov_res, wqcov_res, wccov_res,
                    wTcov_res, mqcoh_res, mccoh_res, qccoh_res, Ebq_cov_res, Ebc_cov_res])
    
    print('%s spectral stats done'%datetime)
    
    return(out)


#%% join spectral stats with master frame, save output

for n in np.arange(df_master.shape[0]):
    spec_out = SpectralStats(df_master.index[n])
    if n == 0:
        out_df = spec_out
    else:
        out_df = out_df.append(spec_out)

#%% save new reference df with all of these stats

df_master = df_master.join(out_df, how='left')
df_master.to_csv(output_fname, index_label='Timestamp')

