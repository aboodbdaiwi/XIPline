#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Reconstruct nutation curve
3/2026
@author: Matt Willmering
"""

import sys
import sonofrecon_cchmc as sor
import subprocess
import datetime

# set parameters
what = sor.get_what()
topdir = sor.dir(what)
dofifo = sor.dofifo()        # use mascri running in background
use_archive = sor.is_ox()    # use ScanArchive (instead of p-file)
delete_pf = False            # delete p-file from mrraw
lb = '[]'
zf = '[]'

# redirect output to logfile
if len(sys.argv)>6:
  now = datetime.datetime.now()
  fnlog = topdir + 'log/' + now.strftime("%Y%m%d_%H%M%S") + '_recon_nutation.txt'
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
mascri = sor.setup_mcr()        # setup mcr
if mascri=='':
  print('Error: mascri not found')
  if len(sys.argv)>6:
    sys.stdout.close()
  sys.exit(-1)
sor.print_env_vars()

# actual execution of matlab-compiled program
cmd = ['recon_nutation',datadir + rawdataname,'[]',zf,lb,'[]','[]']
if dofifo:
  sor.run_fifo(cmd)
  # sor.view(datadir + rawdataname)
else:
  sor.run(cmd)

# clean up mrraw
if delete_pf:
  sor.remove_pfile(runnum)

# finishing up
print('recon_nutation.py done')
if len(sys.argv)>6:
  sys.stdout.close()

sys.exit(0)

