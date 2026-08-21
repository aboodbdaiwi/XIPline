#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
2D/3D gridding reconstruction for gas exchange
3/2026
@author: Matt Willmering
"""

import sys
import sonofrecon_cchmc as sor
import subprocess
import datetime


# set parameters
fidalldir = sor.dir('fidall')
use_vre = False
dofifo = sor.dofifo()        # use mascri running in background (if use_vre==False)
use_archive = sor.is_ox()    # use ScanArchive (instead of p-file)
mtx = '[]'                   # matrix resolution to reconstruct to
interp_factor = '2'          # factor to interpolate
delay = sor.delay()          # gradient-acquisition delay [us]
lbt = '3'                    # Gaussian linebroadening [Hz] %3 is rough number to approximate offline
mask_fac = '[]'              # Threshold factor for masking
coco = '[]'                  # Coil combination (see mri_coil_combine.m)
plt = '[]'                   # Plotting
grdwrp = sor.grdwrp()        # apply gradwarp distortion correction

# redirect output to logfile
if len(sys.argv)>6:
  now = datetime.datetime.now()
  fnlog = fidalldir + 'log/' + now.strftime("%Y%m%d_%H%M%S") + '_recon_grid3d_gx.txt'
  sys.stdout = open(fnlog,'w')

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
wfn = sor.get_wfname(functionname,sor.datadir(exam,series))
if wfn is None:
  print('Error: waveform name not found')
  if len(sys.argv)>6:
    sys.stdout.close()
  sys.exit(-1)

# copy rawdata
if use_archive:
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

# start matlab reconstruction on host
mascri = sor.setup_mcr()        # setup mcr
if mascri=='':
  print('Error: mascri not found')
  if len(sys.argv)>6:
    sys.stdout.close()
  sys.exit(-1)
sor.print_env_vars()

# start matlab reconstruction
cmd = ['recon_grid3d_gx',rawin,'[]',wfn,mtx,interp_factor,delay,lbt,datadir+rawdataname,mask_fac,coco,plt,grdwrp]
if not use_vre:
  # on host
  if dofifo:
    sor.run_fifo(cmd)
  else:
    sor.run(cmd)
else:
  # on vre
  sor.copy_mcr()
  sor.run_on_vre(cmd)
  sor.import_dcm(datadir)              # move dicom to import directory
  sor.move(datadir,datadir2)           # moving files
  datadir = datadir2

# display results on host
sor.view(datadir)

# finishing up
print('_recon_grid3d_gx.py done')
if len(sys.argv)>6:
  sys.stdout.close()

sys.exit(0)

