#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Reconstruct zte data
Created: 12/2025
@author: Rolf Schulte
"""

import sys
import sonofrecon as sor
import subprocess
import datetime


# set parameters
fidalldir = sor.dir('fidall')
use_vre = True
dofifo  = False              # use mascri running in background (if use_vre==False)
use_archive = sor.is_ox()    # use ScanArchive (instead of p-file)
mtx = '[]'                   # matrix resolution to reconstruct to
delay = sor.delay()          # gradient-acquisition delay [us]
lbt = '0'                    # Gaussian linebroadening [Hz]
plt = '[]'                   # Plotting
grdwrp = sor.grdwrp()        # apply gradwarp distortion correction


# redirect output to logfile
if len(sys.argv)>6:
  now = datetime.datetime.now()
  fnlog = fidalldir + 'log/' + now.strftime("%Y%m%d_%H%M%S") + '_recon_zte.txt'
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
# recon_zte(d,h,mtx_reco,delay,lb,fname,coil_arr,plt,grdwrp,overrecofa,max_waspi,save_coils)
cmd = ['recon_zte',rawin,'[]',mtx,delay,lbt,datadir+rawdataname,'[]',plt,grdwrp]
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

# copy dicom data
sor.copy_dicom(exam,series,datadir+'dcm/')

# display results on host
sor.view(datadir)

# finishing up
print('recon_zte.py done')
if len(sys.argv)>6:
  sys.stdout.close()

sys.exit(0)

