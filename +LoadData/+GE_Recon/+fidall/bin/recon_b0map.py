#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
gridding/Cartesian reconstruction with subsequent B0map calculation
Created: 8/2025
@author: Rolf Schulte
"""

import sys
import os
import sonofrecon as sor
import subprocess
import datetime


# set parameters
fidalldir = sor.dir('fidall')
dofifo = sor.dofifo()        # use mascri running in background
do_grid = True               # recon_grid.m or recon_cart.m
use_archive = sor.is_ox()    # use ScanArchive (instead of p-file)
start_hos = True             # startup HOS tool
mtx = '[]'                   # matrix resolution to reconstruct to
delay = sor.delay()          # gradient-acquisition delay [us]
lbt = '0'                    # Gaussian linebroadening [Hz]
mask_fac = '[]'              # Threshold factor for masking
coco = '3'                   # Coil combination (see mri_coil_combine.m)
plt1 = '[0 0]'               # plotting in background, no dicom export (recon_cart/recon_grid)
if start_hos:
  # plt2 = '[0 2 0 1]'         # calc_b0map: plt in background, dcm export b0map+bbabs, no dcm df0 norm
  plt2 = '[0 0 0 1 0]'
else:
  plt2 = '[0 2 0 0]'         # (4)=save data in HOS tool format
TE = '[0 2.3d-3]'            # echo-times for df0 scaling [s]
grdwrp = 'false'             # apply gradwarp distortion correction
do_unwrap = 'true'           # phase unwrapping for B0map

# gui_tk = '/usr/g/bin/gui.tk'
gui_tk = fidalldir + 'bin/gui.tk' # HOS tool; create symbolic link in fidall/bin first

# redirect output to logfile
if len(sys.argv)>6:
  now = datetime.datetime.now()
  fnlog = fidalldir + 'log/' + now.strftime("%Y%m%d_%H%M%S") + '_recon_b0map.txt'
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
fname = datadir+rawdataname
print('rawdataname = ')
print(rawdataname)
if rawdataname is None:
  print('Error: no rawdata found; exiting')
  if len(sys.argv)>6:
    sys.stdout.close()
  sys.exit(-1)

# execution of matlab functions
if do_grid:
  cmd1 = ['recon_grid',fname,'[]',wfn,mtx,delay,lbt,fname,mask_fac,coco,plt1,grdwrp]
else:
  cmd1 = ['recon_cart',fname,'[]',wfn,mtx,delay,lbt,fname,'[]',grdwrp,plt1]
cmd2 = ['calc_b0map',fname,'[]',TE,fname,plt2,'[]','[]','[]',do_unwrap]

if dofifo:
  sor.run_fifo(cmd1)
  sor.run_fifo(cmd2)
else:
  sor.run(cmd1)
  sor.run(cmd2)

# display results on host
sor.view(datadir)

# startup HOS tool
if start_hos:
  os.environ['TCL_LIBRARY'] = ''
  os.environ['HO_SHIM_DEBUG'] = '1'
  if fname.endswith('.h5'):
    fname = fname[0:-3]
  if fname.endswith('.7'):
    fname = fname[0:-2]    
  if not sor.hos_reference(fname + '_hos/P99999.7'):
    print('Error: cannot generate/link reference for HOS tool')
    sys.exit(-1)
  cmd = [gui_tk,fname + '_hos/P99999.7']
  sor.run_os(cmd)

# finishing up
print('recon_b0map.py done')
if len(sys.argv)>6:
  sys.stdout.close()

sys.exit(0)

