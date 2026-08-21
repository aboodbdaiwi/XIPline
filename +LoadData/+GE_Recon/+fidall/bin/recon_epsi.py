#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Reconstruct fidepsi encoded EPSI data
8/2024
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
lbt = '5'                    # spactral Gaussian linebroadening [Hz]
ws = '0'                     # 0=off;10=automatically determine ecc + water suppression
verb = '2'


# redirect output to logfile
if len(sys.argv)>6:
  now = datetime.datetime.now()
  fnlog = topdir + 'log/' + now.strftime("%Y%m%d_%H%M%S") + '_recon_epsi.txt'
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

# start matlab reconstruction on host
mascri = sor.setup_mcr([],what)
if mascri=='':
  print('Error: mascri not found')
  if len(sys.argv)>6:
    sys.stdout.close()
  sys.exit(-1)
sor.print_env_vars()

# actual execution of matlab-compiled program
cmd = ['recon_csi_all',datadir+rawdataname,'[]','[]',lbt,'[]',ws,datadir+rawdataname,'[]',verb]
sor.run(cmd,mascri)

# finishing up
print('recon_epsi.py done')
if len(sys.argv)>6:
  sys.stdout.close()

sys.exit(0)

