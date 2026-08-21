#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Cartesian reconstruction for regular MRI
Created: 3/2025
@author: Rolf Schulte
"""

import sys
import sonofrecon as sor
import subprocess
import datetime


# set parameters
what = sor.get_what()
topdir = sor.dir(what)
use_archive = sor.is_ox()    # use ScanArchive (instead of p-file)
mtx = '-1'                   # matrix resolution to reconstruct to
plt = '0'                    # visible plotting: 0=off(in background), 1=on
coco = '[]'                  # coil combination method
lb = '[]'                    # spatial apodisation
grdwrp = sor.grdwrp()        # apply gradwarp distortion correction

# redirect output to logfile
if len(sys.argv)>6:
  now = datetime.datetime.now()
  fnlog = topdir + 'log/' + now.strftime("%Y%m%d_%H%M%S") + '_recon_mri.txt'
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
print('datadir = ')
print(datadir)

# copy rawdata
if use_archive:
  rawdataname = sor.download_scanarchive(exam,series,runnum,datadir)
else:
  rawdataname = sor.copy_pfile(runnum,datadir)
print('rawdataname = ')
print(rawdataname)
if rawdataname is None:
  print('Error: no rawdata found; exiting')
  if len(sys.argv)>6:
    sys.stdout.close()
  sys.exit(-1)

# start matlab reconstruction on host
mascri = sor.setup_mcr([],what)        # setup mcr
if mascri=='':
  print('Error: mascri not found')
  if len(sys.argv)>6:
    sys.stdout.close()
  sys.exit(-1)
sor.print_env_vars()

# actual execution of matlab-compiled program
cmd = ['recon_mri',datadir+rawdataname,'[]',mtx,datadir+rawdataname,plt,coco,lb,grdwrp]
sor.run(cmd,mascri)

# copy of product spiral recon dicom (executes after matlab window is closed)
sor.copy_dicom(exam,series,datadir)

# finishing up
print('recon_mri.py done')
if len(sys.argv)>6:
  sys.stdout.close()

sys.exit(0)

