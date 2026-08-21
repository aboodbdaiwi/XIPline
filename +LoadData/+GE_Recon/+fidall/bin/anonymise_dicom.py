#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Anonymise dicom files
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
  print('anonymise_dicom.py fnames')
  print('fnames Dicom files to anonymise')
  sys.exit(-1)

# get input parameters
if len(sys.argv)==2:
  fnames = str(sys.argv[1])
else:
  # convert list into semi-colon-separated string
  fnames = ';'.join(sys.argv[1:-1])

print(fnames)
# actual execution of matlab-compiled program
cmd = ['anonymise_dicom',fnames]
if dofifo:
  sor.run_fifo(cmd)
else:
  sor.run(cmd)

sys.exit(0)

