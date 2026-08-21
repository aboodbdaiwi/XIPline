#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Function to read 1H centre frequency from dicom image and convert into various MNS frequencies
Created: 10/2025
@author: Rolf Schulte
Installation as app:
cd fidall/bin
./app_install.py 1 "MNS Freq" mnsfreq_app.py "Determine MNS Frequencies from 1H Dicom"
"""

import sys
import os
import subprocess
import re
import sonofrecon as sor
import datetime


# set parameters
fidall_log = '/export/home/sdc/fidall/log/'
if os.path.isdir(fidall_log):
  fnlog = fidall_log
else:
  dname = os.path.dirname(os.path.realpath(__file__))+'/'
  print(dname)
  if not os.path.isdir(dname + 'log'):
    os.mkdir(dname + 'log')
  fnlog = dname + 'log/'
now = datetime.datetime.now()
mflog = fnlog + 'mnsfreqs.log'
fnlog = fnlog + now.strftime("%Y%m%d_%H%M%S") + '_mnsfreq_app.log'

# open logfile
sys.stdout = open(fnlog,"wt")
sys.stderr = sys.stdout
print('Input arguments ' + str(sys.argv))

# function help
if len(sys.argv)<2:
  sor.xmessage('Error: at least 1 input arguments required')
  sys.stdout.close()
  sys.exit(-1)

# get info from selected dicom file
val = sor.get_dcm_hdr(['(0x0018,0x0084)','(0x0018,0x0085)', \
                       '(0x0043,0x1002)','(0x0043,0x1003)','(0x0043,0x1004)'])   # cfreq,nucleus,x+y+z shim
print(val)
if len(val)<5:
  sor.xmessage('Error in get_dcm_hdr; exiting')
  sys.stdout.close()
  sys.exit(-1)

# parameters
cfreq = float(val[0])*1000000
nucleus = val[1]
xshim = val[2]
yshim = val[3]
zshim = val[4]

# checks
if nucleus != '1H':
  sor.xmessage('Error: nucleus (=\''+nucleus+'\) != \'1H\'')
  sys.stdout.close()
  sys.exit(-1)
if cfreq<1:
  sor.xmessage('Error: cfreq<1; aborting')
  sys.stdout.close()
  sys.exit(-1)

# mns frequency ratios
mns_tab = [\
    ['1H ','H2O',       1],\
    ['2H ','D2O',       0.153506071852142],\
    ['13C','lactate',   1/3.976224812],\
    ['13C','pyruvate',  1/3.976271443],\
    ['13C','bicarb',    1/3.976311538],\
    ['13C','urea',      1/3.976302499380294],\
    ['13C','oil',       0.25145],\
    ['13C','peg',       0.251465],\
    ['23Na','',         0.264517727834842],\
    ['31P','Ph buff',   0.404806204],\
    ['31P','PCr',       0.404803954],\
    ['129Xe','gas bag', 0.2766032749],\
    ['129Xe','gas lung',0.2766029834],\
    ['129Xe','blood',0.2766632951],\
    ['129Xe','tissue',0.2766577631]]
 
print(mns_tab)
print(len(mns_tab))

# convert to MNS freq
txt = ''
for ll in range(len(mns_tab)):
  f0 = int(round(cfreq*mns_tab[ll][2]))
  line = mns_tab[ll][0] + ' ' + mns_tab[ll][1]
  line = line.ljust(14)
  if not f0>99999999:
    line = line + ' '
  txt = txt + line + str(f0) + ' [Hz]'
  if ll+1<len(mns_tab):
    txt = txt + '\n'

# shim info
txt = txt + '\n----------------------------\n'
txt = txt + 'Shim: x=' + xshim + ',y=' + yshim + ',z=' + zshim


# write results into mnsfreqs.log
f = open(mflog,'a')
f.write('----- ' + now.strftime("%Y%m%d_%H%M%S") + ' -----\n')
f.write(txt + '\n')
f.close()

# display results in message box
sor.xmessage(txt)

# close logfile
sys.stdout.close()

# main function exit code
sys.exit(0)


