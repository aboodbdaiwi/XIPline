#!/usr/bin/python2
# -*- coding: utf-8 -*-
"""
Interface function to run reconstructions on vre
Created: 9/2025
@author: Rolf Schulte
"""

import sys
import sonofrecon as sor
import subprocess
import datetime
import time
import re
import os


# set parameters
max_mascri = 1          # maximum number of matlab_scripter processes allowed (0=no limit)
fidalldir = sor.dir('fidall')

# always redirect output to logfile
now = datetime.datetime.now()
fnlog = fidalldir + 'log/' + now.strftime("%Y%m%d_%H%M%S") + '_recon_vre.txt'
sys.stdout = open(fnlog,'w')

# print function help
print(sys.argv)
sys.stdout.flush()

# redirect display to host for robustness
if True:
  os.environ['DISPLAY'] = '10.0.1.1:0.0'

# setup mcr
mascri = sor.setup_mcr()
if mascri=='':
  print('Errror: mascri not found')
  sys.stdout.close()
  sys.exit(-1)
sor.print_env_vars()

# limit #processes running
if max_mascri>0:
  while True:
    # procs = subprocess.check_output(['ps', '-a', '-c', '-x', '-ocomm=']).splitlines()
    # count = procs.count('matlab_scripter')
    output = subprocess.Popen(['ps','-e'],stdout=subprocess.PIPE).communicate()[0]
    output = output.decode('ISO-8859-1')
    count = len(re.findall('matlab_scripter',output))
    print('#matlab_scripter processes running = ' + str(count))
    sys.stdout.flush()
    if count<max_mascri:
      break
    else:
      time.sleep(10)

# actual execution of matlab compiled program
cmdstr = sys.argv[1:]
cmdstr.insert(0,mascri)
print('cmdstr = ')
print(cmdstr)
sys.stdout.flush()
mascriproc = subprocess.Popen(cmdstr,shell=False,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
(stdoutdata, stderrdata) = mascriproc.communicate()
stdoutdata = stdoutdata.decode('ISO-8859-1')
stderrdata = stderrdata.decode('ISO-8859-1')
print('stdoutdata =\n' + str(stdoutdata))
print('stderrdata =\n' + str(stderrdata))

# finishing up
print('recon_vre2.py done')
sys.stdout.close()

sys.exit(0)

