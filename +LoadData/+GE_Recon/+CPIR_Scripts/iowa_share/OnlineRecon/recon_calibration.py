#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
2D gridding reconstruction 
Created: 7/2020
@author: Rolf Schulte
"""

import sys
import sonofrecon_iowa as sor
import subprocess


# set parameters
fidalldir = sor.dir('fidall')
mtx = '0'               # matrix resolution to reconstruct to; 0=acq mtx
mtx_fac = '0'           # factor for matrix resolution to reconstruct to; 0=acq mtx
delay = '0'             # gradient-acquisition delay [us]
comb_time = '0'         # RMS sum-up of time-resolved images

# redirect output to logfile
if len(sys.argv)>6:
  sys.stdout = open(fidalldir+'log/recon_calibration.txt','w')

# print help
if sor.print_help(sys.argv):
  if len(sys.argv)>6:
    sys.stdout.close()
  sys.exit(-1)

# get input parameters
(runnum,exam,series,functionname) = sor.function_input(sys.argv)

# define data directory
datadir = sor.datadir(exam,series)

# get waveform name from symbolic link
wfn = sor.get_wfname(functionname)
if wfn is None:
  print('Error: waveform name not found')
  if len(sys.argv)>6:
    sys.stdout.close()
  sys.exit(-1)

# copy rawdata
if sor.is_ox():
  rawdataname = sor.download_scanarchive(exam,series,runnum,datadir)
else:
  rawdataname = sor.copy_pfile(runnum,datadir)
print(rawdataname)
if rawdataname is None:
  print('Error: no rawdata found; exiting')
  if len(sys.argv)>6:
    sys.stdout.close()
  sys.exit(-1)

# start matlab reconstruction on host
mascri = sor.setup_mcr()        # setup mcr
if mascri=='':
  print('Error: mascri not found')
  if len(sys.argv)>6:
    sys.stdout.close()
  sys.exit(-1)
sor.print_env_vars()

# actual execution of matlab-compiled program
cmd = ['recon_calibration',datadir+rawdataname,'[]',wfn,delay,datadir+rawdataname]
sor.run(cmd,mascri)

# finishing up
print('recon_calibration.py done')
if len(sys.argv)>6:
  sys.stdout.close()

sys.exit(0)

