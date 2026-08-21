#!/usr/bin/python
# 6/2021 Rolf Schulte

import sys
import os
import subprocess


# function help
if len(sys.argv)<4:
  print('Error: at least three input argument required')
  print('Add/remove button to GE Apps section')
  print('Usage:\n\tapp_install.py install appname fname description')
  print('\tinstall=0    deinstall app')
  print('\tinstall=1    install app')
  print('\tappname      Name of app button')
  print('\tfname        Function to call')
  print('\tdescription  Button description (optional)')
  sys.exit(-1)

# input parameters
install = int(sys.argv[1])
appname = sys.argv[2]
fname = sys.argv[3]
if len(sys.argv)>4:
  desc = sys.argv[4]
else:
  desc = 'Run ' + fname

# command
if install==0:   # deinstall app
  print('Deinstalling app ' + appname)  
  cmd = ['addApp2System','-m','FLOAT','-r',appname]
else:         # install app
  print('Installing app ' + appname)
  if not os.path.isfile(fname):
    print('Error: fname=\''+fname+'\' not found')
    sys.exit(-1)
  fpname = os.path.abspath(fname)
  cmd = ['addApp2System','-d',desc,'-n','MULT_INST_ALLOWED','-m','FLOAT','-t','n', \
    '-s','SELECTED','-v',appname,fpname]

# execute
print(cmd)
subprocess.call(cmd)

# main function exit code
sys.exit(0)

