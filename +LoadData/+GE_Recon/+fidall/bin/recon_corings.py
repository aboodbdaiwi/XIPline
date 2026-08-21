#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Reconstruct concentric rings MRSI 
Created: 7/2025
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
lbt = '1'                    # Gaussian linebroadening [Hz]
zf = '[-2 -2 -2 -2]'         # zero-fill by factor of two in all dim
# zf = '[]'


# redirect output to logfile
if len(sys.argv)>6:
  now = datetime.datetime.now()
  fnlog = fidalldir + 'log/' + now.strftime("%Y%m%d_%H%M%S") + '_recon_corings.txt'
  sys.stdout = open(fnlog,'w')

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
wfn = sor.get_wfname(functionname,datadir)
if wfn is None:
  print('Error: waveform name not found')
  if len(sys.argv)>6:
    sys.stdout.close()
  sys.exit(-1)

# copy rawdata
if use_archive:
  rawdataname = sor.download_scanarchive(exam,series,runnum,\
    sor.datadir(exam,series))
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
cmd = ['recon_corings',rawin,'[]',wfn,zf,lbt,datadir+rawdataname]
# on host
sor.run(cmd)

# finishing up
print('recon_corings.py done')
if len(sys.argv)>6:
  sys.stdout.close()

sys.exit(0)

