#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Plot spectra to visualise noise
7/2025
@author: Rolf Schulte
"""

import sys
import sonofrecon as sor
import subprocess
import datetime


# set parameters
# plot_spec_noise(d,h,sv,dofreq,thresh,nmax)
fidalldir = sor.dir('fidall')
dofifo = sor.dofifo()        # use mascri running in background
use_archive = sor.is_ox()    # use ScanArchive (instead of p-file)
dofreq = '[]'                # Plot x-axis in absolute frequencies                  (true)
thresh = '[]'                # Factor for printing maxima: thresh*mean(spec(:))        (2)
nmax = '[]'                  # Restrict maximu to nmax                                (20)

# redirect output to logfile
if len(sys.argv)>6:
  now = datetime.datetime.now()
  fnlog = fidalldir + 'log/' + now.strftime("%Y%m%d_%H%M%S") + '_plot_spec_noise.txt'
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

# actual execution of matlab-compiled program
cmd = ['plot_spec_noise',datadir + rawdataname,'[]',datadir + rawdataname,dofreq,thresh,nmax]
if dofifo:
  sor.run_fifo(cmd)
else:
  sor.run(cmd)

# finishing up
print('plot_spec_noise.py done')
if len(sys.argv)>6:
  sys.stdout.close()

sys.exit(0)

