#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Reconstruct density-weighted MRSI for dissolved-phase xenon
Created: 3/2023
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
use_vre = False
plt_mrsi = 'false'

# redirect output to logfile
if len(sys.argv)>6:
  now = datetime.datetime.now()
  fnlog = fidalldir + 'log/' + now.strftime("%Y%m%d_%H%M%S") + '_recon_mrsi_xenon.txt'
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
wfn = sor.get_wfname(functionname)
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
cmd = ['recon_mrsi_xenon',rawin,'[]',wfn,datadir+rawdataname,plt_mrsi]
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
  cmd.append('3')
  sor.run_on_vre(cmd)  
  sor.import_dcm(datadir)              # move dicom to import directory
  sor.move(datadir,datadir2)           # moving files
  datadir = datadir2

  # display results on host
  sor.view(datadir)
  
  # plot_mrsi
  mascri = sor.setup_mcr()
  rdn = rawdataname
  if re.search('\.h5$',rdn):
    rdn = rdn[0:-3]
  if re.search('\.7$',rdn):
    rdn = rdn[0:-2]
  sor.run(['plot_mrsi_file',datadir+rdn+'.mat'],mascri)


# finishing up
print('recon_mrsi_xenon.py done')
if len(sys.argv)>6:
  sys.stdout.close()

sys.exit(0)

