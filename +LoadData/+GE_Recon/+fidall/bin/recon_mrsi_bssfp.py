#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Reconstruct density-weighted MRSI 
8/2025
@author: Rolf Schulte
"""

import sys
import sonofrecon as sor
import subprocess
import os
import re
import datetime


# set parameters
fidalldir = sor.dir('fidall')
use_archive = sor.is_ox()    # use ScanArchive (instead of p-file)
use_vre = True
# zf = '[-2 -2 -2 -2]'       # zero-fill by factor of two in all dim
zf = '[]'
lbt = '1'                    # Gaussian linebroadening [Hz]
plt_mrsi = False             # plot_mrsi
coco = '[]'                  # coil combination
slowft = '2'                 # 0=slow FT, 1=3D gridding, 2=pinv


# redirect output to logfile
if len(sys.argv)>6:
  now = datetime.datetime.now()
  fnlog = fidalldir + 'log/' + now.strftime("%Y%m%d_%H%M%S") + '_recon_mrsi.txt'
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

# start matlab reconstruction
if not use_vre:
  # on host
  if plt_mrsi:
    plt = '2'
  else:
    plt = '1'
  cmd = ['recon_mrsi',rawin,'[]',wfn,zf,lbt,datadir+rawdataname,plt,coco,slowft,\
    '[]','[]','[]','[]','3','100']
  sor.run(cmd)
else:
  # on vre
  cmd = ['recon_mrsi',rawin,'[]',wfn,zf,lbt,datadir+rawdataname,'0',coco,slowft,\
    '[]','[]','[]','[]','3','100']
  sor.copy_mcr()
  sor.run_on_vre(cmd)  
  sor.import_dcm(datadir)              # move dicom to import directory
  sor.move(datadir,datadir2)           # moving files
  datadir = datadir2

  # display results on host
  sor.view(datadir)

  # plot_mrsi
  if plt_mrsi:
    rdn = rawdataname
    if re.search('\.h5$',rdn):
      rdn = rdn[0:-3]
    if re.search('\.7$',rdn):
      rdn = rdn[0:-2]
    sor.run(['plot_mrsi_file',datadir+rdn+'.mat'])


# finishing up
print('recon_mrsi.py done')
if len(sys.argv)>6:
  sys.stdout.close()

sys.exit(0)

