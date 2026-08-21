#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Copy dicom files to fidall data
12/2023
@author: Rolf Schulte
"""

import sys
import sonofrecon as sor
import subprocess
import datetime
import shutil
import os


# set parameters
what = sor.get_what()
topdir = sor.dir(what)

# redirect output to logfile
if len(sys.argv)>6:
  now = datetime.datetime.now()
  fnlog = topdir + 'log/' + now.strftime("%Y%m%d_%H%M%S") + '_archive_dicom.txt'
  sys.stdout = open(fnlog,'w')

# print help
if sor.print_help(sys.argv):
  if len(sys.argv)>6:
    sys.stdout.close()
  sys.exit(-1)

# get input parameters
(runnum,exam,series,functionname) = sor.function_input(sys.argv)

# define data directory
datadir = sor.datadir(exam,series,what)

# copy dicom files
sor.copy_dicom(exam,series,datadir)

# finishing up
print('archive_dicom.py done')
if len(sys.argv)>6:
  sys.stdout.close()

sys.exit(0)

