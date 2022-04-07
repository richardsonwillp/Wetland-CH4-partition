#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  9 15:30:57 2021

@author: willrichardson
"""

# run this at beginning of partitioning to save a record of the settings used in each partitioning run

import os
import datetime
print('Run start time: %s'%datetime.datetime.now().strftime('%Y-%m-%d %H:%M'))
print('\nSite-year ID: %s' %os.environ['site_yr'])
print('Partitioning run ID: %s' %os.environ['run_ID'])
print('Time zone of the site: %s' %os.environ['tz_local'])
print('Time zone in which data were logged: %s' %os.environ['tz_logging'])
print('\nFILE PATH INFORMATION')
print('Output from your typical flux processing pipeline located in ' + os.environ['HH_output'])
print('Other stats needed located in ' + os.environ['HH_stats'])
print('Non-ebullitive periods listed in ' + os.environ['NonEb_periods'])
print('Canopy height data located in ' + os.environ['zc_file'])
print('\nSITE PARAMETERS')
print('Measurement height = %s m' %os.environ['Zm'])
print('eddy covariance system acquisition rate = %s s' %os.environ['DT'])
print('user-supplied u_star threshold = %s m s-1' %os.environ['ust_thresh'])
print('minimum allowed wind direction = %s degrees' %os.environ['WD_min'])
print('maximum allowed wind direction = %s degrees' %os.environ['WD_max'])
print('fetch expands %s m from tower' %os.environ['fetch_x'])
print('minimum CH4 flux to ensure adequate scalar similarity for partitioning: %s umol m-2 s-1' %os.environ['FCH4_min'])
print('user-supplied minimum LE to ensure adequate scalar similarity for partitioning: %s W m-2' %os.environ['LE_thresh'])
print('user-supplied minimum Fco2 to ensure adequate scalar similarity for partitioning: %s umol m-2 s-1' %os.environ['Fco2_thresh'])
print('filtering for non-local processes? %s' %os.environ['NL_filt'])
print('\nPARTITIONING PARAMETERS')
print('lower frequency bound = %s' %os.environ['LB'])
print('upper frequency bound = %s' %os.environ['UB'])
print('Ebullition threshold type: %s' %os.environ['EbThresh_type'])
print('Ebullition threshold scaling quantity: %s' %os.environ['EbThresh_qty'])

base_path = os.getcwd()
save_fname = base_path + '/ref_data/Output/%s/%s_%s_PartRunInfo.txt' %(os.environ['site_yr'], os.environ['site_yr'], os.environ['run_ID'])

part_settings = ['Run start time: %s'%datetime.datetime.now().strftime('%Y-%m-%d %H:%M'),
                 'Site-year ID: %s' %os.environ['site_yr'],
                 'Partitioning run ID: %s' %os.environ['run_ID'],
                 'Time zone of the site: %s' %os.environ['tz_local'],
                 'Time zone in which data were logged: %s' %os.environ['tz_logging'],
                 'FILE PATH INFORMATION',
                 'Output from your typical flux processing pipeline located in ' + os.environ['HH_output'],
                 'Other stats needed located in ' + os.environ['HH_stats'],
                 'Non-ebullitive periods listed in ' + os.environ['NonEb_periods'],
                 'Canopy height data located in ' + os.environ['zc_file'],
                 'SITE PARAMETERS',
                 'Measurement height = %s m' %os.environ['Zm'],
                 'eddy covariance system acquisition rate = %s s' %os.environ['DT'],
                 'user-supplied u_star threshold = %s m s-1' %os.environ['ust_thresh'],
                 'minimum allowed wind direction = %s degrees' %os.environ['WD_min'],
                 'maximum allowed wind direction = %s degrees' %os.environ['WD_max'],
                 'fetch expands %s m from tower' %os.environ['fetch_x'],
                 'minimum CH4 flux to ensure adequate scalar similarity for partitioning: %s umol m-2 s-1' %os.environ['FCH4_min'],
                 'user-supplied minimum LE to ensure adequate scalar similarity for partitioning: %s W m-2' %os.environ['LE_thresh'],
                 'user-supplied minimum Fco2 to ensure adequate scalar similarity for partitioning: %s umol m-2 s-1' %os.environ['Fco2_thresh'],
                 'filtering for non-local processes? %s' %os.environ['NL_filt'],
                 'PARTITIONING PARAMETERS',
                 'lower frequency bound = %s' %os.environ['LB'],
                 'upper frequency bound = %s' %os.environ['UB'],
                 'Ebullition threshold type: %s' %os.environ['EbThresh_type'],
                 'Ebullition threshold scaling quantity: %s' %os.environ['EbThresh_qty']]

outfile = open(save_fname, mode='w')
for text in part_settings:
    outfile.write(text + '\n')

outfile.close()