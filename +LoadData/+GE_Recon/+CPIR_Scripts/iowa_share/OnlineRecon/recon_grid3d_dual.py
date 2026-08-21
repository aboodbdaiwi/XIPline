#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
3D gridding reconstruction 
Created: 7/2020
@author: Rolf Schulte
"""

import sys
import sonofrecon_iowa as sor
import subprocess
import os


# set parameters
fidalldir = sor.dir('fidall')
use_vre = False
mtx = '0'               # matrix resolution to reconstruct to; 0=acq mtx
delay = '0'             # gradient-acquisition delay [us]
lbt = '0'               # Gaussian linebroadening [Hz]
spks = 'false'          # Spike noise detection and removal
do_gascor = 'false' 	# Gas contamination removal flag

# redirect output to logfile
if len(sys.argv)>6:
  sys.stdout = open(fidalldir+'log/recon_grid3d_dual.txt','w')

# print help
if sor.print_help(sys.argv):
  if len(sys.argv)>6:
    sys.stdout.close()
  sys.exit(-1)

# get input parameters
(runnum,exam,series,functionname) = sor.function_input(sys.argv)

# define data directory
if use_vre:
  datadir = sor.datadir_tmp(exam,series)  
  datadir2 = sor.datadir(exam,series)
else:
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
  rawdataname = sor.download_scanarchive(exam,series,runnum,\
    sor.datadir(exam,series))
  if use_vre:
    rawin = sor.vredatadir(exam,series)+rawdataname
  else:
    rawin = datadir+rawdataname
else:
  rawdataname = sor.copy_pfile(runnum,datadir)
  rawin = datadir+rawdataname
print('rawdataname = ')
print(rawdataname)
if rawdataname is None:
  print('Error: no rawdata found; exiting')
  if len(sys.argv)>6:
    sys.stdout.close()
  sys.exit(-1)

# start matlab reconstruction
cmd = ['recon_grid3d_dual',rawin,'[]',wfn,mtx,delay,lbt,spks,datadir+rawdataname,'[]','[]',do_gascor,'3']
# and prep for offline GX processing
cmd_prep = ['prep_dixon_gx',rawin,'[]',wfn,mtx,delay,lbt,spks,datadir+rawdataname,'[]',do_gascor]
if not use_vre:
  # on host
  mascri = sor.setup_mcr()
  if mascri=='':
    print('Error: mascri not found')
    if len(sys.argv)>6:
      sys.stdout.close()
    sys.exit(-1)
  sor.print_env_vars()
  sor.run(cmd,mascri)
  sor.run(cmd_prep,mascri)
else:
  # on vre
  sor.copy_mcr()
  sor.run_on_vre(cmd)
  sor.import_dcm(datadir)              # move dicom to import directory
  sor.move(datadir,datadir2)           # moving files
  datadir = datadir2

# display results on host
cmdstr = 'eog ' + datadir + '*.png &'
print(cmdstr)
os.system(cmdstr)


# finishing up
print('recon_grid3d_dual.py done')
if len(sys.argv)>6:
  sys.stdout.close()

sys.exit(0)

