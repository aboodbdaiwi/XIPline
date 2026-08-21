#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Function for plotting reconstructed fidall MRSI data
Created: 7/2025
@author: Rolf Schulte
Installation as app:
cd fidall/bin
./app_install.py 1 "Plot MRSI" plot_mrsi.py "Plot fidall-reconstructed MRSI data"
"""

import sys
import sonofrecon as sor
import datetime
import os

# set parameters
print('Input arguments ' + str(sys.argv))
fidalldir = sor.dir('fidall')

# print help
if len(sys.argv)<2:
  print('plot_mrsi.py fnin output2stdout')
  print('fnin Reconstructed fidall MRSI mat file: isempty -> gui')
  sys.exit(-1)

# redirect output to logfile
if len(sys.argv)==2:
  now = datetime.datetime.now()
  fnlog = fidalldir + 'log/' + now.strftime("%Y%m%d_%H%M%S") + '_plot_mrsi.txt'
  sys.stdout = open(fnlog,'w')

# get input parameters
datadir = fidalldir+'data/'
fnin = str(sys.argv[1])
if 'SDC_SELECTION_FILE' in fnin:
  print('Called as app: set fnin=[]')
  fnin = '[]'
  val = sor.get_dcm_hdr(['(0x0020,0x0010)','(0x0020,0x0011)'])   # exam,series
  print('Exam='+val[0]+' Series='+val[1])
  if os.path.isdir(datadir+'Exam'+str(val[0])+'/Series'+str(val[1])):
    datadir = datadir+'Exam'+str(val[0])+'/Series'+str(val[1])

  # goto fidall/data directory for good file selection gui default
  print('Changing dir to '+datadir);
  os.chdir(datadir)

# actual execution of matlab-compiled program
cmd = ['plot_mrsi_file',fnin]
sor.run(cmd)

# finishing up
print('plot_mrsi.py done')
if len(sys.argv)==2:
  sys.stdout.close()

sys.exit(0)

