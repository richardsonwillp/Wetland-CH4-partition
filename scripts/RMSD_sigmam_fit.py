#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  4 16:12:15 2021

@author: willrichardson
"""
# This script fits relationships between RMSD and sigma_m for various selections of periods
# in order to estimate an ebullition threshold width:
    #(1) fit linear regression using only your manually ID'd non-ebullitive periods
    #(2) fit a selection of quantile regressions using all qc'd periods
    
# should be passed arguments for:
    # (1) qc'd dataframe filepath, (2) RMSD/RLM output filename, (3) text file with list of non-ebullitive timestamps
    # (1) site-year, (2) path to non-ebulltive periods (3) fn_LB, (4) fn_UB, 

#%% user variables
percentiles = [0.05, 0.15, 0.25, 0.33, 0.40, 0.50, 0.75, 0.95]
QR_c_list = ['firebrick', 'orangered', 'orange', 'gold', 'yellowgreen', 'mediumspringgreen', 'royalblue', 'darkblue']

#%% import libraries
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
from matplotlib import patches, lines
from funcs import LR_xy_NC, QuantReg_xy_NC, ax_Add_bf_lines, ax_AddLine, ax_AddPoints_NoC
from funcs import sigma_m

#%% get relevant environment variables
site_yr = os.environ['site_yr']
run_ID = os.environ['run_ID']
fn_LB = os.environ['LB']
fn_UB = os.environ['UB']
NL_filt = bool(os.environ['NL_filt'])

# get working directory
base_path = os.getcwd()

#%% Load data
if NL_filt == True:
    q_df_filt = pd.read_csv(base_path + '/ref_data/Output/%s/%s_Interm_RefDf_qPart_filtHH_NLfilt.csv' %(site_yr, site_yr),
                          index_col='Timestamp', infer_datetime_format=True, parse_dates=True)
    
    c_df_filt = pd.read_csv(base_path + '/ref_data/Output/%s/%s_Interm_RefDf_cPart_filtHH_NLfilt.csv' %(site_yr, site_yr),
                          index_col='Timestamp', infer_datetime_format=True, parse_dates=True)
else:
    q_df_filt = pd.read_csv(base_path + '/ref_data/Output/%s/%s_Interm_RefDf_qPart_filtHH.csv' %(site_yr, site_yr),
                          index_col='Timestamp', infer_datetime_format=True, parse_dates=True)
    
    c_df_filt = pd.read_csv(base_path + '/ref_data/Output/%s/%s_Interm_RefDf_cPart_filtHH.csv' %(site_yr, site_yr),
                          index_col='Timestamp', infer_datetime_format=True, parse_dates=True)    

RMSD_df = pd.read_csv(base_path + '/ref_data/Output/%s/%s_RLMfitRMSDs_fnLB_%s_fnUB_%s.csv' %(site_yr, site_yr, fn_LB, fn_UB),
                      index_col=0, infer_datetime_format=True, parse_dates=True)
# join these two into one dataframe
q_df_main = q_df_filt.join(RMSD_df)
c_df_main = c_df_filt.join(RMSD_df)

# non-ebullitive periods
NonEb_periods = pd.read_csv(os.environ['NonEb_periods'], sep='\t', index_col=0, header=None, infer_datetime_format=True, parse_dates=True)
qNonEb_df = NonEb_periods.join(q_df_main, how='left')
cNonEb_df = NonEb_periods.join(c_df_main, how='left')
# check for any rows with all NaNs (i.e., user-selectd non-eb period didn't meet all quality requirements), drop
qNonEb_df.dropna(axis=0, how='all'); cNonEb_df.dropna(axis=0, how='all')

#%% Make output filenames

# ebullition threshold parameters
EbThresh_fname = base_path + '/ref_data/Output/%s/%s_RMSDx-sigmaCH4_slopes_EbThreshWidth_fnLB_%s_fnUB_%s.csv' %(site_yr, site_yr, fn_LB, fn_UB)
# figure showing RMSD - simga_m relationships and ebullition threshold widths
fig_name = base_path + '/ref_data/Output/%s/%s/%s_RMSDx-sigmaCH4_plot_fnLB_%s_fnUB_%s' %(site_yr, run_ID, site_yr, fn_LB, fn_UB)
fig_name1 = base_path + '/ref_data/Output/%s/%s/%s_RMSDx-sigmaCH4_ust_plot_fnLB_%s_fnUB_%s' %(site_yr, run_ID, site_yr, fn_LB, fn_UB)


#%% do regressions

# Linear regression between RMSD - sigma(m), all periods
## i.e.,  mean response of all periods
RMSDq_sdm_res, RMSDq_sdm_bf = LR_xy_NC(q_df_main, sigma_m, 'RMSD_mq')
RMSDc_sdm_res, RMSDc_sdm_bf = LR_xy_NC(c_df_main, sigma_m, 'RMSD_mc')
RMSDT_sdm_res, RMSDT_sdm_bf = LR_xy_NC(q_df_main, sigma_m, 'RMSD_mT')
print('RMSDq - sigmaCH4 R2, all periods: %.2f; RMSDc - sigmaCH4 R2, all periods: %.2f; RMSDT - sigmaCH4 R2, all periods: %.2f\n' %(RMSDq_sdm_res.rsquared, RMSDc_sdm_res.rsquared, RMSDT_sdm_res.rsquared))

# Linear regression between RMSD - sigma(m), 'non-ebullitive' periods
## i.e.,  mean response of non-ebullitive periods
NonEb_RMSDq_sdm_res, NonEb_RMSDq_sdm_bf = LR_xy_NC(qNonEb_df, sigma_m, 'RMSD_mq')
NonEb_RMSDc_sdm_res, NonEb_RMSDc_sdm_bf = LR_xy_NC(cNonEb_df, sigma_m, 'RMSD_mc')
NonEb_RMSDT_sdm_res, NonEb_RMSDT_sdm_bf = LR_xy_NC(qNonEb_df, sigma_m, 'RMSD_mT')
print('RMSDq - sigmaCH4, non-ebullitive periods: y = %.2fx, R2 = %.2f\n' %(NonEb_RMSDq_sdm_res.params[0], NonEb_RMSDq_sdm_res.rsquared)) 
print('RMSDc - sigmaCH4, non-ebullitive periods: y = %.2fx, R2 = %.2f\n' %(NonEb_RMSDc_sdm_res.params[0], NonEb_RMSDc_sdm_res.rsquared)) 
print('RMSDT - sigmaCH4, non-ebullitive periods: y = %.2fx, R2 = %.2f\n' %(NonEb_RMSDT_sdm_res.params[0], NonEb_RMSDT_sdm_res.rsquared)) 

#% Quantile regression between RMSD - sigma(m), all periods
##  get response of all quantiles from 0.05 to 0.95
QR_RMSDq_sdm_res, QR_RMSDq_sdm_bf = QuantReg_xy_NC(q_df_main, sigma_m, 'RMSD_mq',
                                                      percentiles)
QR_RMSDc_sdm_res, QR_RMSDc_sdm_bf = QuantReg_xy_NC(c_df_main, sigma_m, 'RMSD_mc',
                                                      percentiles)
QR_RMSDT_sdm_res, QR_RMSDT_sdm_bf = QuantReg_xy_NC(q_df_main, sigma_m, 'RMSD_mT',
                                                      percentiles)

# Repeat for sigma_m/u* as the x variable
RMSDq_sdm_ust_res, RMSDq_sdm_ust_bf = LR_xy_NC(q_df_main, 'sigmaCH4_ust', 'RMSD_mq')
RMSDc_sdm_ust_res, RMSDc_sdm_ust_bf = LR_xy_NC(c_df_main, 'sigmaCH4_ust', 'RMSD_mc')
RMSDT_sdm_ust_res, RMSDT_sdm_ust_bf = LR_xy_NC(q_df_main, 'sigmaCH4_ust', 'RMSD_mT')
print('RMSDq - sigmaCH4 R2, all periods: %.2f; RMSDc - sigmaCH4 R2, all periods: %.2f; RMSDT - sigmaCH4 R2, all periods: %.2f\n' %(RMSDq_sdm_ust_res.rsquared, RMSDc_sdm_ust_res.rsquared, RMSDT_sdm_ust_res.rsquared))

NonEb_RMSDq_sdm_ust_res, NonEb_RMSDq_sdm_ust_bf = LR_xy_NC(qNonEb_df, 'sigmaCH4_ust', 'RMSD_mq')
NonEb_RMSDc_sdm_ust_res, NonEb_RMSDc_sdm_ust_bf = LR_xy_NC(cNonEb_df, 'sigmaCH4_ust', 'RMSD_mc')
NonEb_RMSDT_sdm_ust_res, NonEb_RMSDT_sdm_ust_bf = LR_xy_NC(qNonEb_df, 'sigmaCH4_ust', 'RMSD_mT')
print('RMSDq - sigmaCH4, non-ebullitive periods: y = %.2fx, R2 = %.2f\n' %(NonEb_RMSDq_sdm_ust_res.params[0], NonEb_RMSDq_sdm_ust_res.rsquared)) 
print('RMSDc - sigmaCH4, non-ebullitive periods: y = %.2fx, R2 = %.2f\n' %(NonEb_RMSDc_sdm_ust_res.params[0], NonEb_RMSDc_sdm_ust_res.rsquared)) 
print('RMSDT - sigmaCH4, non-ebullitive periods: y = %.2fx, R2 = %.2f\n' %(NonEb_RMSDT_sdm_ust_res.params[0], NonEb_RMSDT_sdm_ust_res.rsquared)) 

QR_RMSDq_sdm_ust_res, QR_RMSDq_sdm_ust_bf = QuantReg_xy_NC(q_df_main, 'sigmaCH4_ust', 'RMSD_mq',
                                                      percentiles)
QR_RMSDc_sdm_ust_res, QR_RMSDc_sdm_ust_bf = QuantReg_xy_NC(c_df_main, 'sigmaCH4_ust', 'RMSD_mc',
                                                      percentiles)
QR_RMSDT_sdm_ust_res, QR_RMSDT_sdm_ust_bf = QuantReg_xy_NC(q_df_main, 'sigmaCH4_ust', 'RMSD_mT',
                                                      percentiles)

#%% make figure to visualize results of fitting

fig, RMSD_sdm_ax = plt.subplots(nrows=1, ncols=2, sharex=True, sharey=True, figsize=(4.7,2))

# do plotting
## q
RMSD_sdm_ax[0].scatter(q_df_main[sigma_m], q_df_main['RMSD_mq'], s=8, facecolors='none', edgecolors='k',
              linewidths=0.7, alpha=0.5)
RMSD_sdm_ax[0].set_ylabel('$RMSD_q$', size=9, labelpad=1)
## c 
RMSD_sdm_ax[1].scatter(c_df_main[sigma_m], c_df_main['RMSD_mc'], s=8, facecolors='none', edgecolors='k',
              linewidths=0.7, alpha=0.5)
RMSD_sdm_ax[1].set_ylabel('$RMSD_c$', size=9, labelpad=7)

# add lines of best fit from QR
RMSD_sdm_ax = ax_Add_bf_lines(RMSD_sdm_ax, [QR_RMSDq_sdm_bf, QR_RMSDc_sdm_bf],
                                  [QR_RMSDq_sdm_res, QR_RMSDc_sdm_res],
                                  QR_c_list)
# add mean response of non-ebullitive periods
RMSD_sdm_ax = ax_AddLine(RMSD_sdm_ax, [NonEb_RMSDq_sdm_res.params[0], NonEb_RMSDc_sdm_res.params[0]], 'non-ebullitive', 'fuchsia')
# add stats for mean response of all periods
RMSD_sdm_ax = ax_AddLine(RMSD_sdm_ax, [RMSDq_sdm_res.params[0], RMSDc_sdm_res.params[0]], 'all', 'grey')

# plot NonEb periods on top in a different color
RMSD_sdm_ax = ax_AddPoints_NoC(RMSD_sdm_ax, [qNonEb_df, cNonEb_df], ['RMSD_mq', 'RMSD_mc'], sigma_m, 'r')

# add R2 annotation in title (for all points and non-eb points)
RMSD_sdm_ax[0].annotate('R$^2$ = %.2f (all points) | n = %i'%(RMSDq_sdm_res.rsquared, RMSDq_sdm_res.nobs), (0.03,1.11), xycoords='axes fraction', size=8)
RMSD_sdm_ax[0].annotate('R$^2$ = %.2f (non-eb) | n = %i'%(NonEb_RMSDq_sdm_res.rsquared, NonEb_RMSDq_sdm_res.nobs), (0.03,1.02), xycoords='axes fraction', size=8, c='fuchsia')
RMSD_sdm_ax[1].annotate('R$^2$ = %.2f (all points) | n = %i'%(RMSDc_sdm_res.rsquared, RMSDc_sdm_res.nobs), (0.02,1.11), xycoords='axes fraction', size=8)
RMSD_sdm_ax[1].annotate('R$^2$ = %.2f (non-eb) | n = %i'%(NonEb_RMSDc_sdm_res.rsquared, NonEb_RMSDc_sdm_res.nobs), (0.02,1.02), xycoords='axes fraction', size=8, c='fuchsia')

y_vars = ['RMSD_mq', 'RMSD_mc']
subplot_labs=['(a)', '(b)']
for i in np.arange(RMSD_sdm_ax.size):
    
    # tick formatting
    RMSD_sdm_ax.reshape(-1)[i].xaxis.set_minor_locator(AutoMinorLocator())
    RMSD_sdm_ax.reshape(-1)[i].yaxis.set_minor_locator(AutoMinorLocator())
    RMSD_sdm_ax.reshape(-1)[i].tick_params(axis='both', which='major', labelsize=7.5, pad=2.5)
    RMSD_sdm_ax.reshape(-1)[i].ticklabel_format(axis='both', style='scientific', scilimits=(0,0), useMathText=True)
    RMSD_sdm_ax.reshape(-1)[i].yaxis.offsetText.set_fontsize(7.5); RMSD_sdm_ax.reshape(-1)[i].xaxis.offsetText.set_fontsize(7.5)
    t = RMSD_sdm_ax.reshape(-1)[i].yaxis.get_offset_text(); t.set_x(-0.2)
    t_x = RMSD_sdm_ax.reshape(-1)[i].xaxis.get_offset_text(); t_x.set_x(1.12)

    # grids
    RMSD_sdm_ax.reshape(-1)[i].grid(which='major', axis='both', alpha=0.4, zorder=0)
    RMSD_sdm_ax.reshape(-1)[i].set_axisbelow(True)

    # other labels
    RMSD_sdm_ax.reshape(-1)[i].set_xlabel('$\sigma_m$ (mmol m$^{-3}$)', size=9, labelpad=0)

# make fancy legends:
    # (1) get portion to the left of colon in each label in RMSDq plot; make legend in far left under plots with handles and these labels
    # (2) get equation portion for each; make centered annotation for each equation underneath respective subplots
q_handles, q_labs = RMSD_sdm_ax[0].get_legend_handles_labels()
c_handles, c_labs = RMSD_sdm_ax[1].get_legend_handles_labels()
quant_labs = []; q_eqns = []; c_eqns = []
for label in q_labs:
    quant_labs.append(label.split(':')[0])
    q_eqns.append(label.split(':')[1][1:])
    
for label in c_labs:
    c_eqns.append(label.split(':')[1][1:])
    
leg = RMSD_sdm_ax[0].legend(q_handles, quant_labs, loc='upper left', bbox_to_anchor=(0.76,-0.2), ncol=1, fontsize=8, handlelength=1.0, handletextpad=0.8, labelspacing=0.3, edgecolor='grey', framealpha=0.0, markerfirst=False)

x_loc = [-40,-40,-40,-40,-40,-40,-40,-40,-1,-100]
x_count = 0
for t in leg.get_texts():
    t.set_position((x_loc[x_count],0))
    x_count = x_count + 1

x_eqn_loc = 0.52
y_loc = [-0.32, -0.41, -0.50, -0.59, -0.68, -0.775, -0.865, -0.955, -1.045, -1.137]
y_idx = 0
for eqn in q_eqns:
    RMSD_sdm_ax[0].annotate(eqn, xy=(x_eqn_loc, y_loc[y_idx]), xycoords='axes fraction', ha='center', size=8)
    y_idx = y_idx + 1

x_eqn_loc = 0.52
y_idx = 0
for eqn in c_eqns:
    RMSD_sdm_ax[1].annotate(eqn, xy=(x_eqn_loc, y_loc[y_idx]), xycoords='axes fraction', ha='center', size=8)
    y_idx = y_idx + 1
    
leg_outline = patches.Rectangle((0.22,-0.79), 0.6, 0.715, ec='grey', fc='none', lw=0.5) 
col_line1 =  lines.Line2D([0.395,0.395], [-0.79,-0.79+0.715], lw=0.5, transform=fig.transFigure, figure=fig, ls='-', c='grey')
col_line2 =  lines.Line2D([0.64,0.64], [-0.79,-0.79+0.715], lw=0.5, transform=fig.transFigure, figure=fig, ls='-', c='grey')
fig.lines.extend([col_line1, col_line2])

y_loc_fig = [-0.155, -0.225, -0.295, -0.365, -0.435, -0.505, -0.575, -0.645, -0.715]
for row_y in y_loc_fig:
    line = lines.Line2D([0.22,0.22+0.6], [row_y,row_y], lw=0.5, transform=fig.transFigure, ls='-', c='grey')
    fig.lines.extend([line])

fig.add_artist(leg_outline)
fig.subplots_adjust(wspace=0.2)

# save figure
plt.savefig(fig_name + '.png', dpi=300, bbox_inches='tight')

# make a zoomed in version and save it also 
xmax = np.nanpercentile(q_df_main[sigma_m], 99); ymax = np.nanpercentile(q_df_main['RMSD_mq'], 99)
RMSD_sdm_ax[0].set_xlim((-0.0001,xmax)) 
RMSD_sdm_ax[0].set_ylim((-0.0001,ymax)) 
plt.savefig(fig_name + '_zoom.png', dpi=300, bbox_inches='tight')

#%% repeat for sigma_m/u* as quantinty 
#   (extend number of decimal points in slopes to 3 in figure)

fig, RMSD_sdm_ax = plt.subplots(nrows=1, ncols=2, sharex=True, sharey=True, figsize=(4.7,2))

# do plotting
## q
RMSD_sdm_ax[0].scatter(q_df_main['sigmaCH4_ust'], q_df_main['RMSD_mq'], s=8, facecolors='none', edgecolors='k',
              linewidths=0.7, alpha=0.5)
RMSD_sdm_ax[0].set_ylabel('$RMSD_q$', size=9, labelpad=1)
## c 
RMSD_sdm_ax[1].scatter(c_df_main['sigmaCH4_ust'], c_df_main['RMSD_mc'], s=8, facecolors='none', edgecolors='k',
              linewidths=0.7, alpha=0.5)
RMSD_sdm_ax[1].set_ylabel('$RMSD_c$', size=9, labelpad=7)

# add lines of best fit from QR
RMSD_sdm_ax = ax_Add_bf_lines(RMSD_sdm_ax, [QR_RMSDq_sdm_ust_bf, QR_RMSDc_sdm_ust_bf],
                                  [QR_RMSDq_sdm_ust_res, QR_RMSDc_sdm_ust_res],
                                  QR_c_list)
# add mean response of non-ebullitive periods
RMSD_sdm_ax = ax_AddLine(RMSD_sdm_ax, [NonEb_RMSDq_sdm_ust_res.params[0], NonEb_RMSDc_sdm_ust_res.params[0]], 'non-ebullitive', 'fuchsia')
# add stats for mean response of all periods
RMSD_sdm_ax = ax_AddLine(RMSD_sdm_ax, [RMSDq_sdm_ust_res.params[0], RMSDc_sdm_ust_res.params[0]], 'all', 'grey')

# plot NonEb periods on top in a different color
RMSD_sdm_ax = ax_AddPoints_NoC(RMSD_sdm_ax, [qNonEb_df, cNonEb_df], ['RMSD_mq', 'RMSD_mc'], 'sigmaCH4_ust', 'r')

# add R2 annotation in title (for all points and non-eb points)
RMSD_sdm_ax[0].annotate('R$^2$ = %.2f (all points) | n = %i'%(RMSDq_sdm_ust_res.rsquared, RMSDq_sdm_ust_res.nobs), (0.03,1.11), xycoords='axes fraction', size=8)
RMSD_sdm_ax[0].annotate('R$^2$ = %.2f (non-eb)'%NonEb_RMSDq_sdm_ust_res.rsquared, (0.03,1.02), xycoords='axes fraction', size=8, c='fuchsia')
RMSD_sdm_ax[1].annotate('R$^2$ = %.2f (all points) | n = %i'%(RMSDc_sdm_ust_res.rsquared, RMSDc_sdm_res.nobs), (0.02,1.11), xycoords='axes fraction', size=8)
RMSD_sdm_ax[1].annotate('R$^2$ = %.2f (non-eb) | n = %i'%(NonEb_RMSDc_sdm_ust_res.rsquared, NonEb_RMSDc_sdm_ust_res.nobs), (0.02,1.02), xycoords='axes fraction', size=8, c='fuchsia')

y_vars = ['RMSD_mq', 'RMSD_mc']
subplot_labs=['(a)', '(b)']
for i in np.arange(RMSD_sdm_ax.size):
    
    # tick formatting
    RMSD_sdm_ax.reshape(-1)[i].xaxis.set_minor_locator(AutoMinorLocator())
    RMSD_sdm_ax.reshape(-1)[i].yaxis.set_minor_locator(AutoMinorLocator())
    RMSD_sdm_ax.reshape(-1)[i].tick_params(axis='both', which='major', labelsize=7.5, pad=2.5)
    RMSD_sdm_ax.reshape(-1)[i].ticklabel_format(axis='both', style='scientific', scilimits=(0,0), useMathText=True)
    RMSD_sdm_ax.reshape(-1)[i].yaxis.offsetText.set_fontsize(7.5); RMSD_sdm_ax.reshape(-1)[i].xaxis.offsetText.set_fontsize(7.5)
    t = RMSD_sdm_ax.reshape(-1)[i].yaxis.get_offset_text(); t.set_x(-0.2)
    t_x = RMSD_sdm_ax.reshape(-1)[i].xaxis.get_offset_text(); t_x.set_x(1.12)

    # grids
    RMSD_sdm_ax.reshape(-1)[i].grid(which='major', axis='both', alpha=0.4, zorder=0)
    RMSD_sdm_ax.reshape(-1)[i].set_axisbelow(True)

    # other labels
    RMSD_sdm_ax.reshape(-1)[i].set_xlabel('$\sigma_m$/$u_*$', size=9, labelpad=0)

# make fancy legends:
    # (1) get portion to the left of colon in each label in RMSDq plot; make legend in far left under plots with handles and these labels
    # (2) get equation portion for each; make centered annotation for each equation underneath respective subplots
q_handles, q_labs = RMSD_sdm_ax[0].get_legend_handles_labels()
c_handles, c_labs = RMSD_sdm_ax[1].get_legend_handles_labels()
quant_labs = []; q_eqns = []; c_eqns = []
for label in q_labs:
    quant_labs.append(label.split(':')[0])
    q_eqns.append(label.split(':')[1][1:])
    
for label in c_labs:
    c_eqns.append(label.split(':')[1][1:])
    
leg = RMSD_sdm_ax[0].legend(q_handles, quant_labs, loc='upper left', bbox_to_anchor=(0.76,-0.2), ncol=1, fontsize=8, handlelength=1.0, handletextpad=0.8, labelspacing=0.3, edgecolor='grey', framealpha=0.0, markerfirst=False)

x_loc = [-40,-40,-40,-40,-40,-40,-40,-40,-1,-100]
x_count = 0
for t in leg.get_texts():
    t.set_position((x_loc[x_count],0))
    x_count = x_count + 1

x_eqn_loc = 0.52
y_loc = [-0.32, -0.41, -0.50, -0.59, -0.68, -0.775, -0.865, -0.955, -1.045, -1.137]
y_idx = 0
for eqn in q_eqns:
    RMSD_sdm_ax[0].annotate(eqn, xy=(x_eqn_loc, y_loc[y_idx]), xycoords='axes fraction', ha='center', size=8)
    y_idx = y_idx + 1

x_eqn_loc = 0.52
y_idx = 0
for eqn in c_eqns:
    RMSD_sdm_ax[1].annotate(eqn, xy=(x_eqn_loc, y_loc[y_idx]), xycoords='axes fraction', ha='center', size=8)
    y_idx = y_idx + 1
    
leg_outline = patches.Rectangle((0.22,-0.79), 0.6, 0.715, ec='grey', fc='none', lw=0.5) 
col_line1 =  lines.Line2D([0.395,0.395], [-0.79,-0.79+0.715], lw=0.5, transform=fig.transFigure, figure=fig, ls='-', c='grey')
col_line2 =  lines.Line2D([0.64,0.64], [-0.79,-0.79+0.715], lw=0.5, transform=fig.transFigure, figure=fig, ls='-', c='grey')
fig.lines.extend([col_line1, col_line2])

y_loc_fig = [-0.155, -0.225, -0.295, -0.365, -0.435, -0.505, -0.575, -0.645, -0.715]
for row_y in y_loc_fig:
    line = lines.Line2D([0.22,0.22+0.6], [row_y,row_y], lw=0.5, transform=fig.transFigure, ls='-', c='grey')
    fig.lines.extend([line])

fig.add_artist(leg_outline)
fig.subplots_adjust(wspace=0.2)

# save figure
plt.savefig(fig_name1 + '.png', dpi=300, bbox_inches='tight')

# make a zoomed in version and save it also 
xmax = np.nanpercentile(q_df_main['sigmaCH4_ust'], 99); ymax = np.nanpercentile(q_df_main['RMSD_mq'], 99)
RMSD_sdm_ax[0].set_xlim((-0.0001,xmax)) 
RMSD_sdm_ax[0].set_ylim((-0.0001,ymax)) 
plt.savefig(fig_name1 + '_zoom.png', dpi=300, bbox_inches='tight')


#%% save slopes to a text file for reference, with appropriate labels

# make index rows (slope type); columns will be 'RMSDq_sigmaCH4_slope' and 'RMSDc_sigmaCH4_slope'
row_ind = ['NonEb_LinReg', 'All_LinReg']; QuantReg_ind = []
for quant in list(np.rint(np.array(percentiles)*100).astype(int)):
    QuantReg_ind.append('QuantReg_q' + str(quant))

row_ind = row_ind + QuantReg_ind

# make empty dataframe
EbThreshWidths_df = pd.DataFrame(index=row_ind, columns=['RMSDq_sigmaCH4_slope', 'RMSDc_sigmaCH4_slope', 'RMSDT_sigmaCH4_slope',
                                                         'RMSDq_sigmaCH4ust_slope', 'RMSDc_sigmaCH4ust_slope', 'RMSDT_sigmaCH4ust_slope'],
                                 dtype=np.float64)
# fill it
EbThreshWidths_df.loc['NonEb_LinReg', 'RMSDq_sigmaCH4_slope'] = NonEb_RMSDq_sdm_res.params[0]
EbThreshWidths_df.loc['NonEb_LinReg', 'RMSDc_sigmaCH4_slope'] = NonEb_RMSDc_sdm_res.params[0]
EbThreshWidths_df.loc['NonEb_LinReg', 'RMSDT_sigmaCH4_slope'] = NonEb_RMSDT_sdm_res.params[0]
EbThreshWidths_df.loc['NonEb_LinReg', 'RMSDq_sigmaCH4ust_slope'] = NonEb_RMSDq_sdm_ust_res.params[0]
EbThreshWidths_df.loc['NonEb_LinReg', 'RMSDc_sigmaCH4ust_slope'] = NonEb_RMSDc_sdm_ust_res.params[0]
EbThreshWidths_df.loc['NonEb_LinReg', 'RMSDT_sigmaCH4ust_slope'] = NonEb_RMSDT_sdm_ust_res.params[0]
EbThreshWidths_df.loc['All_LinReg', 'RMSDq_sigmaCH4_slope'] = RMSDq_sdm_res.params[0]
EbThreshWidths_df.loc['All_LinReg', 'RMSDc_sigmaCH4_slope'] = RMSDc_sdm_res.params[0]
EbThreshWidths_df.loc['All_LinReg', 'RMSDT_sigmaCH4_slope'] = RMSDT_sdm_res.params[0]
EbThreshWidths_df.loc['All_LinReg', 'RMSDq_sigmaCH4ust_slope'] = RMSDq_sdm_ust_res.params[0]
EbThreshWidths_df.loc['All_LinReg', 'RMSDc_sigmaCH4ust_slope'] = RMSDc_sdm_ust_res.params[0]
EbThreshWidths_df.loc['All_LinReg', 'RMSDT_sigmaCH4ust_slope'] = RMSDT_sdm_ust_res.params[0]

for lab in QuantReg_ind:
    perc = str(int(lab[10:])/100)
    if len(perc) < 4:
        perc = perc + '0'
    EbThreshWidths_df.loc[lab, 'RMSDq_sigmaCH4_slope'] = QR_RMSDq_sdm_res[perc].params[0]
    EbThreshWidths_df.loc[lab, 'RMSDc_sigmaCH4_slope'] = QR_RMSDc_sdm_res[perc].params[0]
    EbThreshWidths_df.loc[lab, 'RMSDT_sigmaCH4_slope'] = QR_RMSDT_sdm_res[perc].params[0]    

    EbThreshWidths_df.loc[lab, 'RMSDq_sigmaCH4ust_slope'] = QR_RMSDq_sdm_ust_res[perc].params[0]
    EbThreshWidths_df.loc[lab, 'RMSDc_sigmaCH4ust_slope'] = QR_RMSDc_sdm_ust_res[perc].params[0]
    EbThreshWidths_df.loc[lab, 'RMSDT_sigmaCH4ust_slope'] = QR_RMSDT_sdm_ust_res[perc].params[0]    

# save it
EbThreshWidths_df.to_csv(EbThresh_fname, index_label='Threshold_type', float_format='%.6f')

