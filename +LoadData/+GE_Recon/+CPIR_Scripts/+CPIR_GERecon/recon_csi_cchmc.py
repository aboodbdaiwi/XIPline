#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Reconstruct Cartesian fidcsi encoded CSI
11/2025
@author: Rolf Schulte
"""

import sys
import sonofrecon_cchmc as sor
import subprocess
import datetime


# set parameters
what = sor.get_what()
topdir = sor.dir(what)
use_archive = sor.is_ox()    # use ScanArchive (instead of p-file)
zf = '[-4, -4, -4, -4]'                    # zero-fill to [#s,#x,#y,#z]
lbt = '15'                    # spectral Gaussian linebroadening [Hz]
plt = '1'                    # plotting: 0=off, 1=normal, 2+=plot_mrsi
coco = '[]'                  # coil combination [do percent_noise4decor]
ws = '0'                     # 0=off;10=automatically determine ecc + water suppression


# redirect output to logfile
if len(sys.argv)>6:
  now = datetime.datetime.now()
  fnlog = topdir + 'log/' + now.strftime("%Y%m%d_%H%M%S") + '_recon_csi_cchmc.txt'
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
print('datadir = ')
print(datadir)

# copy noise scan
if sor.docopynoise():
  sor.copy_noise(exam,datadir)

# copy rawdata
if use_archive:
  rawdataname = sor.download_scanarchive(exam,series,runnum,datadir)
else:
  rawdataname = sor.copy_pfile(runnum,datadir)
print('rawdataname = ')
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
cmd = ['recon_csi_cchmc',datadir+rawdataname,'[]',zf,lbt,datadir+rawdataname,plt,coco,ws]
sor.run(cmd,mascri)

# finishing up
print('recon_csi_cchmc.py done')
if len(sys.argv)>6:
  sys.stdout.close()

sys.exit(0)

