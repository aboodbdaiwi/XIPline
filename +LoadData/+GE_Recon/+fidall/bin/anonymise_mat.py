#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Anonymise Matlab mat files
remove patient data from h (header) structure 
Created: 7/2025
@author: Rolf Schulte
"""

import sys
import sonofrecon as sor
import subprocess


# set parameters
fidalldir = sor.dir('fidall')
dofifo = sor.dofifo()        # use mascri running in background

# print help
if len(sys.argv)<2:
  print('anonymise_mat.py fname')
  print('fname of Matlab mat file containing h structure')
  sys.exit(-1)

# get input parameters
fname = str(sys.argv[1])

# actual execution of matlab-compiled program
cmd = ['anonymise_mat',fname]
if dofifo:
  sor.run_fifo(cmd)
else:
  sor.run(cmd)

sys.exit(0)

