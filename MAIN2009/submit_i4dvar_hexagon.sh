#!/bin/csh -f
#
# svn $Id: submit_i4dvar.sh 585 2012-01-03 18:44:28Z arango $
#######################################################################
# Copyright (c) 2002-2012 The ROMS/TOMS Group                         #
#   Licensed under a MIT/X style license                              #
#   See License_ROMS.txt                                              #
################################################## Hernan G. Arango ###
#                                                                     #
#  North Sea 8x8km IS4DVAR assimilation MyOceanII                     #
#                                                                     #
#  This script is used to run ROMS incremental 4DVar algorithm in     #
#  sequential mode through several assimilation cycles. The user      #
#  needs to have the following directory structure:                   #
#                                                                     #
#    $MYROOT/                                                         #
#           /DATA                                                     #
#           /FORWARD                                                  #
#           /IS4DVAR                                                  #
#           /OBS                                                      #
#           /STORAGE                                                  #
#           /MAIN                                                     #
#  and storage directory:                                             #
#                                                                     #
#    $STORAGE                                                         #
#                                                                     #
#  To submit a job in the batch queue.  Use the following command     #
#  in MPI applications to avoid running on the head node NO_LOCAL:    #
#                                                                     #
#      batch now -f submit.sh                                         #
#                                                                     #
#  To check batch use:                                                #
#                                                                     #
#      bbq                                                            #
#                                                                     #
#######################################################################
#
#PBS -N "COPERNICUS"
#PBS -A imr
#PBS -l mppwidth=512
#PBS -l walltime=420:00:00
#PBS -l mppmem=1000mb
#PBS -m abe
#PBS -M trond.kristiansen@imr.no

echo "  "
echo "**************************************************************"
echo "***     ROMS/TOMS Incremental, Strong Constraint 4DVar     ***"
echo "***  Master Execution Script: Sequential State Estimation  ***"
echo "**************************************************************"
echo "***"

#---------------------------------------------------------------------
#  Directories.
#---------------------------------------------------------------------

#  Set ROOT of the directory to run 4DVar.

set MYROOT="/work/users/trondk/NS8km"
set FORCINGROOT="/work/shared/imr/NS8KM/FORCING"
#  Set ROMS/TOMS ROOT directory.

set ROMS_ROOT="/work/users/trondk/src/KATE_ROMS"

#  Set storage directory for some of the relevant output NetCDF files.

set STORAGE="/work/users/trondk/NS8km/STORAGE"

#---------------------------------------------------------------------
#  Application title and IO file prefix.
#---------------------------------------------------------------------

set TITLE="North_Sea_8x8km_IS4DVAR_assimilation_NOWMAPS"

set PREFIX="ns8km"

echo "***  $TITLE"
echo "***"
echo "   "

#---------------------------------------------------------------------
#  Input files.
#---------------------------------------------------------------------

# Set grid NetCDF file.

set GRDname=${FORCINGROOT}/GRID/nordsjoen_8km_grid_hmax20m_v4.nc

# Set Open boundary conditions file, if any.

set BRYname=${FORCINGROOT}/GLORYS2V3/nordsjoen_8km_bry_GLORYS_20091115_to_20131215.nc

# Set Climatology file

set CLMname=${FORCINGROOT}/GLORYS2V3/nordsjoen_8km_clim_GLORYS_20091115_to_20131215.nc

# Set starting sequential assimilation first guess.

set FIRST_GUESS=${FORCINGROOT}/GLORYS2V3/ns8km_ini_22998.nc

# Set background-error covariance standard deviations file.

set STDname=${FORCINGROOT}/OBS/NS8KM_std_2.0.nc

# Set background-error covariance normalization factor file

set NRMname=${FORCINGROOT}/OBS/NS8KM_nrm_i.nc

# Set observations file.

set OBSname=${FORCINGROOT}/OBS/NS8KM_obsSST_2009_to_2012.nc

#---------------------------------------------------------------------
#  Executables and standard input files
#---------------------------------------------------------------------

#  Set ROMS nonlinear and data assimilation executables.

set NL_ROMS="nl_oceanM"
set DA_ROMS="da_oceanM"

#  Set ROMS nonlinear and data assimilation standard input scripts.

set NL_TEMPLATE=nl_ocean_roms37.tmp
set DA_TEMPLATE=da_ocean_roms37.tmp

set NL_STDINP=nl_ocean_${PREFIX}.in
set DA_STDINP=da_ocean_${PREFIX}.in

#  Set ROMS Metatada variables file.

set VARINFO=varinfo.dat
set VARINFO_KATE=varinfo.dat

#  Set 4DVar input script.

set IS4DVAR_TEMPLATE=s4dvar.tmp

set IS4DVAR_PARAM=s4dvar.in

#  Set string manipulations perl script.

set SUBSTITUTE=${ROMS_ROOT}/ROMS/Bin/substitute

#---------------------------------------------------------------------
#  Time window to consider.
#---------------------------------------------------------------------

#  Set  starting and ending year day of the sequential data assimilation.
#  (Reference time: days since 2006-01-01 00:00:00)

set STR_DAY=22998               # July 12, 2006 00:00:00 UTC/GMT
set END_DAY=24106

#  Set data assimilation cycle time window (days).

set DayStep=7

#---------------------------------------------------------------------
#  Set few Parameters.
#---------------------------------------------------------------------

#  Set model parallel partition.

set NtileI=32
set NtileJ=16

#  Set number of parallel nodes to use, NCPUS = NtileI * NtileJ.

set NCPUS=512

#  Set number of outer and inner loops.

set Nouter=1
set Ninner=7

#  Set number of timesteps to write RESTART file.  This is VERY
#  Important since we are using the restart file of the nonlinear
#  model run as the first guess fot the next assimilation cycle.
#  It MUST be equal to NTIMES.

set NRST=6048 #6048 #3024
set DT=100 #100
set DELTA=48 #16 #4 #86400/320 = time steps in 24 hours. Divide that by 10 or so to get deltaDT = 27
set AVGDT=864 #432 #864 #432 #86400/DT

#  Set enviromental variables to avoid running in the head node.

setenv NO_LOCAL 1
setenv EXCLUDE 10

######################################################################
#  Start sequential data assimilation
######################################################################

set cycle=0

set DAY=$STR_DAY

#  Set starting initial conditions file name.

set INIname=${PREFIX}_ini_${DAY}.nc
set ITLname=${PREFIX}_itl_${DAY}.nc

# Delete old log files
/bin/rm ${MYROOT}/STORAGE/*log*
# Delete old error and output files.
/bin/rm -f *.in~ *.tmp~ NS8KM.e* NS8KM.o*

while ($DAY <= $END_DAY)

  @ cycle += 1

  echo ">>> Starting data assimilation cycle: $cycle"
  echo ">>>"

#---------------------------------------------------------------------
# Run 4DVar Algorithm.
#---------------------------------------------------------------------

  cd $MYROOT/MAIN2009

  # Clean directory by removing all existing NetCDF files.

  if (-e $ITLname) then
    /bin/rm -f $MYROOT/MAIN2009/*.nc
  endif
  set ITLname=${PREFIX}_itl_${DAY}.nc

# Set backgound (first guess) state file.

  if ($DAY == $STR_DAY) then
    cp -p $FIRST_GUESS $INIname
    echo "Copied initial file ", $FIRST_GUESS, " to ", $INIname, " day: ", $DAY
  else
    cp -p ${STORAGE}/$INIname .
    echo "Copied initial file ", ${STORAGE}/$INIname, " to ", $INIname, " day: ", $DAY
  endif

# Set tangent linear model initial conditions file.

  cp -p ${MYROOT}/MAIN2009/DATA/${PREFIX}_ini_zero.nc $ITLname

# Get a clean copy of the observation file.  This is really
# important since this file is modified to compute the
# fractional vertical position of the observations when
# they are specified as depth in meter (negative values).

  cp -p $OBSname .

# Modify 4DVar template input script and specify above files.
# FOr NS8KM we s4dvar.tmp as template script

  if (-e $IS4DVAR_PARAM) then
    /bin/rm $IS4DVAR_PARAM
  endif
  cp $IS4DVAR_TEMPLATE $IS4DVAR_PARAM

  $SUBSTITUTE $IS4DVAR_PARAM ocean_std_i.nc $STDname
  $SUBSTITUTE $IS4DVAR_PARAM ocean_nrm_i.nc $NRMname
  $SUBSTITUTE $IS4DVAR_PARAM ocean_obs.nc $OBSname
  $SUBSTITUTE $IS4DVAR_PARAM ocean_mod.nc ${PREFIX}_mod_${DAY}.nc

# Modify 4DVar ROMS standard input script.

  if (-e $DA_STDINP) then
    /bin/rm $DA_STDINP
  endif
  cp $DA_TEMPLATE $DA_STDINP

  $SUBSTITUTE $DA_STDINP MyTITLE $TITLE
  $SUBSTITUTE $DA_STDINP varinfo.dat $VARINFO
  $SUBSTITUTE $DA_STDINP MyNtileI $NtileI
  $SUBSTITUTE $DA_STDINP MyNtileJ $NtileJ
  $SUBSTITUTE $DA_STDINP MyNouter $Nouter
  $SUBSTITUTE $DA_STDINP MyNinner $Ninner
  $SUBSTITUTE $DA_STDINP MyDSTART ${DAY}.0d0
  $SUBSTITUTE $DA_STDINP MyNRST $NRST
  $SUBSTITUTE $DA_STDINP MyDT $DT
  $SUBSTITUTE $DA_STDINP MyAVGDT $AVGDT
  $SUBSTITUTE $DA_STDINP MyDELTA $DELTA

  $SUBSTITUTE $DA_STDINP ocean_grd.nc $GRDname
  $SUBSTITUTE $DA_STDINP ocean_ini.nc $INIname
  $SUBSTITUTE $DA_STDINP ocean_clm.nc $CLMname
  $SUBSTITUTE $DA_STDINP ocean_itl.nc $ITLname
  $SUBSTITUTE $DA_STDINP ocean_bry.nc $BRYname
  $SUBSTITUTE $DA_STDINP ocean_fwd.nc ${PREFIX}_fwd_${DAY}.nc

  $SUBSTITUTE $DA_STDINP ocean_rst.nc ${PREFIX}_rst_${DAY}.nc
  $SUBSTITUTE $DA_STDINP ocean_his.nc ${PREFIX}_his_${DAY}.nc
  $SUBSTITUTE $DA_STDINP ocean_avg.nc ${PREFIX}_avg_${DAY}.nc
  $SUBSTITUTE $DA_STDINP ocean_tlm.nc ${PREFIX}_tlm_${DAY}.nc
  $SUBSTITUTE $DA_STDINP ocean_adj.nc ${PREFIX}_adj_${DAY}.nc
  $SUBSTITUTE $DA_STDINP s4dvar.in ${MYROOT}/MAIN2009/$IS4DVAR_PARAM

# Run incremental 4DVar algorithm.

  echo ">>> Running IS4DVAR algorithm, starting day: $DAY"

  if (-e da_log.$DAY) then
    /bin/rm -f da_log.$DAY
  endif
  aprun -B $DA_ROMS $DA_STDINP > da_log.$DAY

# Move estimated initial conditions, misfit, and log files to storage.

  echo ">>> Done running IS4DVAR, moving initial conditions to storage"


  # Assuming NTIMES=1680 and DT=360----------------------------
  #rm -f $INIname
  #ncks -d ocean_time,126 ${PREFIX}_his_${DAY}.nc $INIname
  mv -f $INIname $STORAGE
  # -----------------------------------------------------------

  mv -f ${PREFIX}_mod_${DAY}.nc $STORAGE
  mv -f da_log.${DAY} $STORAGE

#---------------------------------------------------------------------
# Run Nonlinear model initialized with 4DVAR estimated initial
# conditions for the period of the assimilation time window. It
# will compute the first guess for the next assimilation cycle
#---------------------------------------------------------------------

  cd ${MYROOT}/MAIN2009

# Create ROMS standard input script from template.

  if (-e $NL_STDINP) then
    /bin/rm $NL_STDINP
  endif
  cp $NL_TEMPLATE $NL_STDINP

  cp $STORAGE/$INIname $INIname

  set RSTname=${PREFIX}_rst_${DAY}.nc

  $SUBSTITUTE $NL_STDINP MyTITLE $TITLE
 # $SUBSTITUTE $NL_STDINP varinfo.dat $VARINFO_KATE
  $SUBSTITUTE $NL_STDINP MyNtileI $NtileI
  $SUBSTITUTE $NL_STDINP MyNtileJ $NtileJ
  $SUBSTITUTE $NL_STDINP MyNRST $NRST
  $SUBSTITUTE $NL_STDINP MyDT $DT
  $SUBSTITUTE $NL_STDINP MyDSTART ${DAY}.0d0
  $SUBSTITUTE $NL_STDINP MyAVGDT $AVGDT
  $SUBSTITUTE $NL_STDINP MyDELTA $DELTA

  $SUBSTITUTE $NL_STDINP ocean_grd.nc $GRDname
  $SUBSTITUTE $NL_STDINP ocean_ini.nc $STORAGE/$INIname
  $SUBSTITUTE $NL_STDINP ocean_clm.nc $CLMname
  $SUBSTITUTE $NL_STDINP ocean_bry.nc $BRYname

  $SUBSTITUTE $NL_STDINP ocean_rst.nc $RSTname
  $SUBSTITUTE $NL_STDINP ocean_his.nc ${PREFIX}_his_${DAY}.nc
  $SUBSTITUTE $NL_STDINP ocean_avg.nc ${PREFIX}_avg_${DAY}.nc



# Run nonlinear ROMS.

  echo ">>> Running nonlinear model, starting day: $DAY"

  if (-e nl_log.${DAY}) then
    /bin/rm -f nl_log.${DAY}
  endif
  aprun -B ${NL_ROMS} ${NL_STDINP} > nl_log.${DAY}

# Move current nonlinear history and log files to storage.

   mv -f ${PREFIX}_his_${DAY}.nc ${STORAGE}
   mv -f ${PREFIX}_avg_${DAY}.nc ${STORAGE}
   mv -f nl_log.${DAY} $STORAGE

#---------------------------------------------------------------------
# Advance starting day for next assimilation cycle. Set new initial
# conditions file name.
#---------------------------------------------------------------------

  @ DAY += $DayStep

  set INIname=${PREFIX}_ini_${DAY}.nc

#---------------------------------------------------------------------
# Move next cycle first guess (background state) to storage. It is
# currently stored in the restart file.
#---------------------------------------------------------------------

  cd ${MYROOT}/MAIN2009

  mv -f $RSTname ${STORAGE}/${INIname}
  echo "Moved", $RSTname, " to ", ${STORAGE}/${INIname}

  echo "  "
  echo ">>> Finished data assimilation cycle: $cycle"

end
