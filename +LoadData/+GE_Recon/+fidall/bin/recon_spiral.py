#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
2D gridding reconstruction 
Created: 3/2025
@author: Rolf Schulte
"""

import sys
import sonofrecon as sor
import subprocess
import datetime


# set parameters
fidalldir = sor.dir('fidall')
use_archive = sor.is_ox()    # use ScanArchive (instead of p-file)
mtx = '128'                  # matrix resolution to reconstruct to; 0=acq mtx
delay = sor.delay()          # gradient-acquisition delay [us]
lbt = '15'                   # Gaussian linebroadening [Hz]
coco = '[]'                  # coil combination
plt = '[1]'                  # plotting
grdwrp = sor.grdwrp()        # apply gradwarp distortion correction

# redirect output to logfile
if len(sys.argv)>6:
  now = datetime.datetime.now()
  fnlog = fidalldir + 'log/' + now.strftime("%Y%m%d_%H%M%S") + '_recon_spiral.txt'
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
# cmd = ['recon_spiral',datadir+rawdataname,'[]',wfn,mtx,delay,lbt,'false',datadir+rawdataname]
# recon_grid(d,h,wfn,mtx_reco,delay,lb,fname,mask_fac,coco,plt,grdwrp,dns)

cmd = ['recon_grid',datadir+rawdataname,'[]',wfn,mtx,delay,lbt,datadir+rawdataname,'[]',coco,plt,grdwrp]
sor.run(cmd)

# finishing up
print('recon_spiral.py done')
if len(sys.argv)>6:
  sys.stdout.close()

sys.exit(0)

