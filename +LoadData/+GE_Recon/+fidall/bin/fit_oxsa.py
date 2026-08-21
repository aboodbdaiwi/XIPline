#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Function for fitting reconstructed fidall MRSI data with OXSA Amares
Created: 7/2025
@author: Rolf Schulte
Installation as app:
cd fidall/bin
./app_install.py 1 "OXSA" fit_oxsa.py "Fit MRS(I) via OXSA Amares"
"""

import sys
import sonofrecon as sor
import datetime
import os

# set parameters
print('Input arguments ' + str(sys.argv))
fidalldir = sor.dir('fidall')

lb = '0'                               # Applied Gaussian linebroadening -> deapodisation [Hz]
fnout = 'true'                         # Save to file (if fname given; or true->fname from fnin)
plt = '3'                              # Plotting (0=off,1=final,2=all,3=random15%)
pk_fname = '[]'                        # OXSA prior knowledge file name
offset = '0'                           # Frequency offset                      [ppm]
time_delay = '0'                       # Acquisition delay (linear phase)      [s]
figoff = '[]'                          # Figure ID (+dcm series) offset

# print help
if len(sys.argv)<2:
  print('fit_oxsa.py fnin output2stdout')
  print('fnin Reconstructed fidall MRSI mat file: isempty -> gui')
  sys.exit(-1)

# redirect output to logfile
if len(sys.argv)==2:
  now = datetime.datetime.now()
  fnlog = fidalldir + 'log/' + now.strftime("%Y%m%d_%H%M%S") + '_fit_oxsa.txt'
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

# start matlab reconstruction on host
mascri = sor.setup_mcr()        # setup mcr
if mascri=='':
  print('Error: mascri not found')
  if len(sys.argv)==2:
    sys.stdout.close()
  sys.exit(-1)
sor.print_env_vars()

# goto fidall/data directory for good file selection gui default
print('Changing dir to '+datadir);
os.chdir(datadir)

# actual execution of matlab-compiled program
cmd = ['fit_oxsa',fnin,'[]',lb,fnout,plt,pk_fname,'[]',offset,time_delay,figoff]
sor.run(cmd,mascri)

# finishing up
print('fit_oxsa.py done')
if len(sys.argv)==2:
  sys.stdout.close()

sys.exit(0)

