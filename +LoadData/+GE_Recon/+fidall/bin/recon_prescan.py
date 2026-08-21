#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Reconstruct MNS prescan acquired with Bloch-Siegert
Created: 11/2025
@author: Rolf Schulte
"""

import sys
import sonofrecon as sor
import subprocess
import datetime


# set parameters
now = datetime.datetime.now()
fidalldir = sor.dir('fidall')
dofifo = sor.dofifo()        # use mascri running in background
use_archive = sor.is_ox()    # use ScanArchive (instead of p-file)
log_fname = fidalldir + 'blosi.log'

# redirect output to logfile
if len(sys.argv)>6:
  textout = fidalldir + 'log/' + now.strftime("%Y%m%d_%H%M%S") + '_recon_prescan.txt'
  sys.stdout = open(textout,'w')

# print help
if sor.print_help(sys.argv):
  if len(sys.argv)>6:
    sys.stdout.close()
  sys.exit(-1)

# get input parameters
(runnum,exam,series,functionname) = sor.function_input(sys.argv)

# define data directory
datadir = sor.datadir(exam,series)

# copy noise scan
if sor.docopynoise():
  sor.copy_noise(exam,datadir)

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

# actual execution of matlab-compiled program
cmd = ['mns_prescan',datadir+rawdataname,'[]',log_fname,'[]',datadir+rawdataname]
if dofifo:
  sor.run_fifo(cmd)
  # sor.view(datadir+rawdataname)
else:
  sor.run(cmd)

# finishing up
print('recon_prescan.py done')
if len(sys.argv)>6:
  sys.stdout.close()

sys.exit(0)

