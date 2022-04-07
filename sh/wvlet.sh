#! /usr/local/bin/tcsh -f

# This shell script conducts the wavelet transformation step of the partitioning program
# Current state: modified for Way4 Summer 2015

# Enter the following variables. 
# WORKDIR: this variable should be set to the full path to your ch4flux_partition- directory, 
# e.g., /home/hiwata/ch4flux_partition-v2.0
set WORKDIR = /Volumes/Backup-Will/Wetland-CH4-partition
# make site-year ID for filenames and subdirectories (avoid using underscores)
set site_yr = Way3-Summer-2015
# DATA: data directories. Multiple data directories can be set for DATA variable,
# e.g., ( 20160823 20160907 20160928 ).
set DATA = ( 20150603 20150605 20150606 20150607 20150609 20150610 20150611 20150612 20150615 20150616 20150617 20150618 20150619 20150620 20150621 20150622 20150623 20150626 20150628 20150629 20150630 20150701 20150702 20150704 20150705 20150706 20150707 20150708 20150709 20150710 20150711 20150712 20150713 20150714 20150715 20150717 20150718 20150719 20150724 20150726 20150728 20150803 20150804 20150805 )
# DT: time interval of time-series data in second (e.g., 10Hz --> 0.1, 20Hz --> 0.05)
set DT = 0.05

cd $WORKDIR
foreach data ( $DATA )
  echo $data
##Get a list of data
  set DATFIL = `ls data/$site_yr/$data/*.dat | awk -F/ '{print $4}'`
  if ( -d wvlet/$site_yr/ ) then
  else
    mkdir -p wvlet/$site_yr/
  endif

  if ( -d wvlet/$site_yr/$data/ ) then
  else
    mkdir -p wvlet/$site_yr/$data/
  endif
##Calculation of each data run
  foreach datfil ( $DATFIL )
    echo $datfil
    cp data/$site_yr/$data/$datfil temp.dat
    set NUM = `wc temp.dat | awk '{print $1}'`
    echo "      parameter ( spl_num =" $NUM ", dt =" $DT ")" > fort/param1
####Wavelet transform
    gfortran fort/wavelet.f ; ./a.out
    mv wvlet.txt wvlet/$site_yr/$data/wvlet-$datfil
    echo " "
  end
end

rm a.out temp.dat
