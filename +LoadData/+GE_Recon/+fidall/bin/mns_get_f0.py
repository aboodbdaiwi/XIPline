#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Extract MNS frequencies from 1H MRS
Created: 8/2024
@author: Rolf Schulte
"""

import sys
import sonofrecon as sor
import subprocess
import datetime


# set parameters
fidalldir = sor.dir('fidall')
use_archive = sor.is_ox()    # use ScanArchive (instead of p-file)

# redirect output to logfile
if len(sys.argv)>6:
  now = datetime.datetime.now()
  fnlog = fidalldir + 'log/' + now.strftime("%Y%m%d_%H%M%S") + '_mns_get_f0.txt'
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

# start matlab reconstruction on host
mascri = sor.setup_mcr()        # setup mcr
if mascri=='':
  print('Error: mascri not found')
  if len(sys.argv)>6:
    sys.stdout.close()
  sys.exit(-1)
sor.print_env_vars()

# actual execution of matlab-compiled program
cmd = ['mns_get_f0',datadir+rawdataname,'[]',datadir+rawdataname]
sor.run(cmd,mascri)

# finishing up
print('mns_get_f0.py done')
if len(sys.argv)>6:
  sys.stdout.close()

sys.exit(0)

