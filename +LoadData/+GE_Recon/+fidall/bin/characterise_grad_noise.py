#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Characterise gradient noise
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
  fnlog = fidalldir + 'log/' + now.strftime("%Y%m%d_%H%M%S") + '_characterise_grad_noise.txt'
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
cmd = ['characterise_grad_noise',datadir+rawdataname,'[]',datadir+rawdataname]
sor.run(cmd,mascri)

# finishing up
print('characterise_gradient_noise.py done')
if len(sys.argv)>6:
  sys.stdout.close()

sys.exit(0)

