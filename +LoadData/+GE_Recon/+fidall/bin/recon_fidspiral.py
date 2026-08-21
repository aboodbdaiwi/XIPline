#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
gridding reconstruction for fidspiral.e
3/2025
@author: Rolf Schulte
"""

import sys
import sonofrecon as sor
import subprocess
import datetime
import shutil


# set parameters
what = sor.get_what()
topdir = sor.dir(what)
use_archive = sor.is_ox()    # use ScanArchive (instead of p-file)
mtx = '[]'                   # matrix resolution to reconstruct to
delay = sor.delay()          # gradient-acquisition delay [us]
lbt = '15'                   # Gaussian linebroadening [Hz]
mask_fac = '[]'              # Threshold factor for masking
coco = '[]'                  # Coil combination (see mri_coil_combine.m)
plt = '0'                    # visible lotting: 0=off(in background), 1=on

# redirect output to logfile
if len(sys.argv)>6:
  now = datetime.datetime.now()
  fnlog = topdir + 'log/' + now.strftime("%Y%m%d_%H%M%S") + '_recon_fidspiral.txt'
  sys.stdout = open(fnlog,'w')

# print help
if sor.print_help(sys.argv):
  if len(sys.argv)>6:
    sys.stdout.close()
  sys.exit(-1)

# get input parameters
(runnum,exam,series,functionname) = sor.function_input(sys.argv)

# define data directory
datadir = sor.datadir(exam,series,what)

# copy rawdata
if use_archive:
  rawdataname = sor.download_scanarchive(exam,series,runnum,datadir)
else:
  rawdataname = sor.copy_pfile(runnum,datadir)
print(rawdataname)
if rawdataname is None:
  print('Error: no rawdata found; exiting')
  if len(sys.argv)>6:
    sys.stdout.close()
  sys.exit(-1)

# copy of product spiral recon dicom
sor.copy_dicom(exam,series,datadir)

# copy trajectory
wfn = datadir + 'fidspiral.grad'
shutil.copy('/export/home/service/log/fidspiral.grad',wfn)

# start matlab reconstruction on host
mascri = sor.setup_mcr([],what)        # setup mcr
if mascri=='':
  print('Error: mascri not found')
  if len(sys.argv)>6:
    sys.stdout.close()
  sys.exit(-1)
sor.print_env_vars()

# actual execution of matlab-compiled program
cmd = ['recon_grid',datadir+rawdataname,'[]',wfn,mtx,delay,lbt,datadir+rawdataname,mask_fac,coco,plt]
sor.run(cmd,mascri)

# finishing up
print('recon_fidspiral.py done')
if len(sys.argv)>6:
  sys.stdout.close()

sys.exit(0)

