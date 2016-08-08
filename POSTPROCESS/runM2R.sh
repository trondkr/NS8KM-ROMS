#!/bin/bash
#
#  Give the job a name
#PBS -N "model2roms_KINO"
#
#  Specify the project the job belongs to
#PBS -A imr
#PBS -q normal
#PBS -l mppwidth=1,walltime=10:00:00
#PBS -l mppmem=1000MB
echo "#PBS -l mppnppn=16"
#
#  Send me an email on  a=abort, b=begin, e=end
#PBS -m abe
#
#  Use this email address (check that it is correct):
#PBS -M trond.kristiansen@imr.no
#
#  Write the standard output of the job to file 'mpijob.out' (optional)
#PBS -o  grd2lonlat.out
#
#  Write the standard error of the job to file 'mpijob.err' (optional)
#PBS -e  grd2lonlat.err
#

#  Make sure I am in the correct directory
cd /work/shared/imr/NS8KM/grid2lonlat
module load python

export MPLCONFIGDIR=${pwd}
export TMP=`pwd`
export PYTHON_EGG_CACHE=/work/shared/imr/NS8KM/grid2lonlat

aprun -B python grid2lonlatZ.py > grid2lonlat.log