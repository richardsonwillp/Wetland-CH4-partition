#! /usr/local/bin/tcsh -f

# Script for running the modified partitioning program
## combines what was formerly the iteration script and the partitioning script,
## further post-processing steps added

##############################################################################################
#					USER MODIFIED VARIABLES			    	     #
##############################################################################################
# currently setup to process on: Way4 Summer 2015 with empirically determined fnLB values

# WORKING DIRECTORY (base path of where program is located)
set WORKDIR = /Volumes/Backup-Will/Wetland-CH4-partition
# make site-year ID for filenames (avoid using underscores)
setenv site_yr Way3-Summer-2015
# Identifier for this partitioning run (useful if doing multiple partitioning runs on one site-year of data)
setenv run_ID Run_test
# Time zones in which data were logged and where site is located
setenv tz_local America/Chicago
setenv tz_logging America/Chicago
# indicator variable for whether or not this site-year has been partitioned before (for isolating qc'd periods, making reference dfs)
set first_run = True

# set a variable named 'DATA' as list of days to run partitioning for; will iterate over this list, running partitioning for all files on each day
set DATA = ( 20150603 20150605 20150606 20150607 20150609 20150610 20150611 20150612 20150615 20150616 20150617 20150618 20150619 20150620 20150621 20150622 20150623 20150626 20150628 20150629 20150630 20150701 20150702 20150704 20150705 20150706 20150707 20150708 20150709 20150710 20150711 20150712 20150713 20150714 20150715 20150717 20150718 20150719 20150724 20150726 20150728 20150803 20150804 20150805 )

# ANCILLARY DATA (use absolute paths)
## 30 minute time series output from EddyPro and any other stats
setenv HH_output /Volumes/Backup-Will/Wetland-CH4-partition/ref_data/Input/Way3-Summer-2015/eddypro_adv_full_output_2021-11-09T211349_adv.csv
setenv HH_stats /Volumes/Backup-Will/Wetland-CH4-partition/ref_data/Input/Way3-Summer-2015/eddypro_Run1_st6_plusQCflags.csv

## text file with list of non-ebullitive periods
setenv NonEb_periods /Volumes/Backup-Will/Wetland-CH4-partition/ref_data/Input/Way3-Summer-2015/Way3-Summer-2015_NonEbPeriods.txt
## canopy height data (daily interval)
setenv zc_file /Volumes/Backup-Will/Wetland-CH4-partition/ref_data/Input/Way3-Summer-2015/Way3-Summer-2015_zc.txt

# PARTITIONING PARAMETERS
## measurement height
setenv Zm 2.20 #[m]

## items related to normalized frequency (bounds, acquisition frequency)
setenv LB 0.003
setenv UB 1.0
setenv DT 0.05 #[s]

## ebullition threshold type, quantity for estimating threshold
### default type set to 'NonEb_LinReg'; modify to other quantiles if desired, see indexes in '/ref_data/Output/site-yr/site-yr_RMSDx-sigmaCH4_slopes_EbThreshWidth*.csv'
### default quantity set to 'sigma_m'; can also use sigma_m normalized by u* by setting to 'sigma_m_ust' 
setenv EbThresh_type NonEb_LinReg
setenv EbThresh_qty sigma_m

# FILTERING CRITERIA
## ustar threshold (plausibility of value is examined in output)
setenv ust_thresh 0.0 # [m s-1]
## wind direction (degrees CW from N; if desired fetch surrounds tower in every direction, set to 0 and 360)
setenv WD_min 95 # [degrees]
setenv WD_max 265 # [degrees]
## approximate distance at which desired fetch stops; if footprint doesn't ever extend past desired area, set to False
setenv fetch_x 350 #[m]
# lowest CH4 flux magnitude at which scalar similarity is maintained (approx 0.01 umol m-2 s-1 in Richardson et al. (2022))
setenv FCH4_min 0.01
# reference scalar thresholds for adequate scalar similarity (visually identified, seems more reliable)
setenv LE_thresh 27 # [W m-2]
setenv Fco2_thresh 6 #[umol m-2 s-1]
# whether or not to filter for non-local processes ()
setenv NL_filt True

##########################################################################################################
# navigate to working directory
cd $WORKDIR 

##########################################################################################################
# 			Make filtered and non-filtered reference dataframes                        	 #
# 		   if this is the first time this site-year has been partitioned 		   	 #
##########################################################################################################

if ( $first_run == True ) then
	echo "first time partitioning this site-year, make reference dataframes"
	mkdir ref_data/Output/$site_yr	
	# make reference frames!
	./scripts/create_reference_frames.py
else
	echo "site-year previously partitioned"
endif

# make subdirectory for reference data outputs
mkdir ref_data/Output/$site_yr/$run_ID

# make variable holding path to sigma_CH4 data
set sigmaCH4_file = `echo "$WORKDIR""/ref_data/Input/""$site_yr""/""$site_yr""_sigmaCH4.txt"`
set sigmaCH4ust_file = `echo "$WORKDIR""/ref_data/Input/""$site_yr""/""$site_yr""_sigmaCH4ust.txt"`

# record settings used in this partitioning run
./scripts/write_part_sttngs.py

########################################################################################################## 
# 	check to see if RMSDs have been fitted for this site-year fn parameter combination
# 		if not, run R script to do so and save in /ref_data/Output
# 	
#       then fit relationships between RMSD and sigmaCH4 to get ebullition threshold values
##########################################################################################################

# get other relevant variables
## RMSD filepath to check for previous partitioning with this fn parameter set
set RMSD_fname = `echo "$WORKDIR""/ref_data/Output/""$site_yr""/""$site_yr""_RLMfitRMSDs_""fnLB_""$LB""_""fnUB_""$UB"".csv"`

if ( ! -f "$RMSD_fname" ) then
	echo "do RMSD fitting"
	Rscript scripts/RLM_fitting.R
	./scripts/RMSD_sigmam_fit.py
else 
	echo "RMSD fitting previously done for these fn boundaries"
endif

## get ebullition thresholds for partitioning
set EbThresh_file = `echo "$WORKDIR""/ref_data/Output/""$site_yr""/""$site_yr""_RMSDx-sigmaCH4_slopes_EbThreshWidth_fnLB_""$LB""_""fnUB_""$UB"".csv"`

## m-ref RMSDs
if ( $EbThresh_qty == sigma_m ) then
	set SD_QCH = `grep $EbThresh_type $EbThresh_file | cut -f2 -d ','`
	set SD_CCH = `grep $EbThresh_type $EbThresh_file | cut -f3 -d ','`
	set SD_TCH = `grep $EbThresh_type $EbThresh_file | cut -f4 -d ','`
	# specify file containing quantity for ebullition threshold
	set RMSD_qty_file = $sigmaCH4_file
	
else if ( $EbThresh_qty == sigma_m_ust ) then
	set SD_QCH = `grep $EbThresh_type $EbThresh_file | cut -f5 -d ','`
	set SD_CCH = `grep $EbThresh_type $EbThresh_file | cut -f6 -d ','`
	set SD_TCH = `grep $EbThresh_type $EbThresh_file | cut -f7 -d ','`
	# specify file containing quantity for ebullition threshold
	set RMSD_qty_file = $sigmaCH4ust_file

endif

## q-T RMSD; default value from original program; deprecated, ignore
set SD_QT = 0.0876

#########################################################################################################
#				 	EXECUTE PARTITIONING 						#
#########################################################################################################

foreach day ( $DATA )
	echo $day
 	set Zc = `grep $day $zc_file | awk '{print $2}'`
  	echo $Zc

  	set DIR = $day
  	set FILE = `ls wvlet/$site_yr/$DIR/*.dat | awk -F/ '{print $4}'`
	echo "      parameter ( zm = " $Zm ")" > fort/param2
	echo "      parameter ( lb = " $LB ", ub = " $UB ")" > fort/param3
	echo "      parameter ( sd_tch = " $SD_TCH ", sd_qch = " $SD_QCH ")" > fort/param4
	echo "      parameter ( sd_qt = " $SD_QT ")" > fort/param5
	echo "      parameter ( sd_cch = " $SD_CCH ")" > fort/param7

	# put the day's canopy height into a temporary file; do the same for DT so it can be read into fortran programs
	echo $Zc > zc.txt
	echo $DT > DT.txt	
	
	# make directories, outer and inner
	if ( -d fig/$site_yr ) then
	else
 		mkdir fig/$site_yr
	endif
	if ( -d fig/$site_yr/$DIR ) then
	else
 		mkdir fig/$site_yr/$DIR
	endif
	
	if ( -d txt/$site_yr ) then
	else
  		mkdir txt/$site_yr
	endif
	if ( -d txt/$site_yr/$DIR ) then
	else
  		mkdir txt/$site_yr/$DIR
	endif

	if ( -d filt/$site_yr ) then
	else
 		mkdir filt/$site_yr
	endif
	if ( -d filt/$site_yr/$DIR ) then
	else
 		mkdir filt/$site_yr/$DIR
	endif

	if ( -d slope/$site_yr ) then
	else
 		mkdir slope/$site_yr
	endif
	if ( -d flux/$site_yr ) then
	else
 		mkdir flux/$site_yr
	endif
	
	echo "file_name total_CH4_flux bubble_CH4_flux_based_on_T-CH4 bubble_CH4_flux_based_on_H2O-CH4 bubble_CH4_flux_based_on_CO2-CH4 low-frequency_component_CH4_flux" > flux.txt
	echo "file_name diffusive_flux_based_on_T-CH4 bubble_CH4_flux_based_on_T-CH4 diffusive_flux_based_on_H2O-CH4 bubble_CH4_flux_based_on_H2O-CH4 diffusive_flux_based_on_CO2-CH4 bubble_CH4_flux_based_on_CO2-CH4" > rerror.txt
	echo "file_name T=H2O*slope CH4=T*slope CH4=H2O*slope CH4=CO2*slope" > slope.txt
	
	# iterate overall files on this day
	foreach file ( $FILE )
		echo $file
		cp wvlet/$site_yr/$DIR/$file temp1.txt
 		set NUM = `wc temp1.txt | awk '{print $1}'`
 		echo "      parameter ( ni = " $NUM ")" > fort/param6
 		echo "u  w  T  q  CO2  CH4" > temp2.txt
  
		# get sigma(CH4) for this period 
  		set tstamp = `echo $file | awk '{print substr($0,7,13)}'`
  		set sigma_CH4 =  `grep $tstamp $RMSD_qty_file | awk '{print $2}'`
		echo $sigma_CH4
  		echo $sigma_CH4 > sigma_CH4.txt
  
  		more temp1.txt >> temp2.txt
  		gfortran fort/pickup.f -o pickup.out
  		./pickup.out
  		Rscript scripts/robust.R
  		mv wt.jpeg fig/$site_yr/$DIR/$file.jpeg
  		mv wt2.jpeg fig/$site_yr/$DIR/$file-2.jpeg
  		mv temp3.txt txt/$site_yr/$DIR/$file
  		echo $file > date.txt
  		paste date.txt out1.txt out2.txt out3.txt out4.txt >> slope.txt
  		gfortran fort/filter.f -o filter.out
  		./filter.out
  		paste date.txt part.txt >> flux.txt
  		gfortran fort/ranerror.f -o ranerror.out
  		./ranerror.out
  		paste date.txt error.txt >> rerror.txt
  		mv filter.txt filt/$site_yr/$DIR/$file
	end

	mv slope.txt slope/$site_yr/slope$DIR.txt
	mv flux.txt flux/$site_yr/flux$DIR.txt
	mv rerror.txt flux/$site_yr/rerror$DIR.txt

	rm filter.out pickup.out ranerror.out date.txt out?.txt part.txt error.txt winddata temp?.txt zc.txt sigma_CH4.txt DT.txt

end

# do post-processing stats
./scripts/post_proc_stats.py $SD_QCH $SD_CCH

# fit reference scalar thresholds and visualize results
Rscript scripts/RefScalar_Thresh.R
./scripts/RefScalar_ThreshViz.py

# save different levels of filtered output
./scripts/filter_save.py

# run post-processing diagnostic plots if this is the first time the site-year has been processed
if ( $first_run == True ) then
	./scripts/PostProc_diagnostics.py
endif
