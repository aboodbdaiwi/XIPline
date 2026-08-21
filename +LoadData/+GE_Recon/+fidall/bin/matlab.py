#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Interface function to start matlab_scripter on host
Created: 7/2025
@author: Rolf Schulte
"""

import sys
import sonofrecon as sor
import os

# set parameters
what = sor.get_what()

# print info if argument == 'help'
print(sys.argv)
if len(sys.argv)>1:
  if sys.argv[1]=='help':
    print('matlab.py <MATLAB FUNCTION> input')
    print('  Run various matlab functions')
    print('  For included matlab functions and more info look into')
    print('    ' + sor.dir('fidall') + 'matlab/')
    print('  Example call:')
    print('      matlab.py \'\"plot(pinv(randn(1,100)));\"\'\n')
    sys.exit(0)

if len(sys.argv)>2:
  print('Error: only one argument allowed; place in brackets')
  sys.exit(-1)

# setup mcr
mascri = sor.setup_mcr([],what)
if mascri=='':
  print('Error: mascri not found')
  sys.exit(-1)

# actual execution of matlab compiled program
if len(sys.argv)>1:
  print(sys.argv)
  arg_str = str(sys.argv[1])
  cmdstr = mascri + ' ' + arg_str
else:
  cmdstr = mascri
print(cmdstr)
os.system(cmdstr)

# finishing up
print('matlab.py done')

sys.exit(0)

