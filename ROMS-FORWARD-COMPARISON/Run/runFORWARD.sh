#!/bin/bash
#
#  Give the job a name
#PBS -N "ns8km_forward"
#
#  Specify the project the job belongs to
#PBS -A imr
#PBS -q normal
#PBS -l mppwidth=128,walltime=96:00:00
#PBS -l mppmem=1000MB
#PBS -l mppnppn=32
#
#  Send me an email on  a=abort, b=begin, e=end
#PBS -m abe
#
#  Use this email address (check that it is correct):
#PBS -M trond.kristiansen@imr.no
#
#  Write the standard output of the job to file 'mpijob.out' (optional)
#PBS -o  ns8km_forward.out
#
#  Write the standard error of the job to file 'mpijob.err' (optional)
#PBS -e  ns8km_forward.err
#

#  Make sure I am in the correct directory
cd /work/shared/imr/NS8KM/FORWARD/Run

aprun -B ./oceanM nl_ocean_ns8km.in > ns8km_forward.log