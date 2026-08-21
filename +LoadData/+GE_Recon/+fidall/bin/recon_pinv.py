#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
2D/3D pinv reconstruction (+ stack-of-2D)
Created: 7/2025
@author: Rolf Schulte
"""

import sys
import sonofrecon as sor
import subprocess
import datetime
import os


# set parameters
fidalldir = sor.dir('fidall')
use_archive = sor.is_ox()    # use ScanArchive (instead of p-file)
mtx = '[]'                   # matrix resolution to reconstruct to
delay = sor.delay()          # gradient-acquisition delay [us]
lbt = '0'                    # Gaussian linebroadening [Hz]
mask_fac = '[]'              # Threshold factor for masking
coco = '[]'                  # Coil combination (see mri_coil_combine.m)
plt = '[]'                   # Plotting
grdwrp = sor.grdwrp()        # apply gradwarp distortion correction

# redirect output to logfile
if len(sys.argv)>6:
  now = datetime.datetime.now()
  fnlog = fidalldir + 'log/' + now.strftime("%Y%m%d_%H%M%S") + '_recon_pinv.txt'
  sys.stdout = open(fnlog,'w')

# print help
if sor.print_help(sys.argv):
  if len(sys.argv)>6:
    sys.stdout.close()
  sys.exit(-1)

# get input parameters
(runnum,exam,series,functionname) = sor.function_input(sys.argv)

# get waveform name from symbolic link
wfn = sor.get_wfname(functionname)
if wfn is None:
  print('Error: waveform name not found')
  if len(sys.argv)>6:
    sys.stdout.close()
  sys.exit(-1)
else:
  Fn = wfn[0:-4] + '_F.h5'
  print('Searching for')
  print(Fn)
  if os.path.isfile(Fn):
    use_vre=False
    useGPU='0'
    print('Recon matrix found, reconstructing on host.')
  else:
    use_vre = True
    useGPU='1'
    print('Recon matrix not found, reconstructing on VRE.')

  
# define data directory
if use_vre:
  datadir = sor.datadir_tmp(exam,series)  
  datadir2 = sor.datadir(exam,series)
else:
  datadir = sor.datadir(exam,series)



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

# start matlab reconstruction
cmd = ['recon_pinv',rawin,'[]',wfn,mtx,delay,lbt,datadir+rawdataname,mask_fac,coco,plt,grdwrp,useGPU]
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
print('recon_pinv.py done')
if len(sys.argv)>6:
  sys.stdout.close()

sys.exit(0)

