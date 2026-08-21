#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Reconstruct Cartesian CSI SNR Maps
Created: 8/2025
@author: Rolf Schulte
"""

import sys
import sonofrecon as sor
import subprocess
import datetime


# set parameters
fidalldir = sor.dir('fidall')
use_archive = sor.is_ox()    # use ScanArchive (instead of p-file)
zf = '[]'                    # zero-filling

# redirect output to logfile
if len(sys.argv)>6:
  now = datetime.datetime.now()
  fnlog = fidalldir + 'log/' + now.strftime("%Y%m%d_%H%M%S") + '_recon_csi_snr.txt'
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
mascri = sor.setup_mcr()
if mascri=='':
  print('Error: mascri not found')
  if len(sys.argv)>6:
    sys.stdout.close()
  sys.exit(-1)
sor.print_env_vars()

# actual execution of matlab-compiled program
cmd = ['recon_csi_snr',datadir+rawdataname,'[]',zf,datadir+rawdataname]
sor.run(cmd,mascri)

# combine output images + info into single png
if True:
  if use_archive:
    fn = datadir + rawdataname[0:-3]
  else:
    fn = datadir + rawdataname[0:-2]
  cmd = ['montage','-tile','2x2','-mode','concatenate',\
   fn+'_snr_ind3d.png',fn+'_snr_txt.png',fn+'_snr_plot.png',fn+'_signal_noise_plot.png',\
   fn+'_snr_all.png']
  sor.run_os(cmd)

# finishing up
print('recon_csi_snr.py done')
if len(sys.argv)>6:
  sys.stdout.close()

sys.exit(0)

