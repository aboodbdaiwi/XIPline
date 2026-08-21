#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Copy raw data and dicom files to fidall data
6/2025
@author: Rolf Schulte
"""

import sys
import sonofrecon as sor
import subprocess
import datetime
import shutil
import os


# set parameters
get_dicom = True                  # copy dicom files
get_pfile = False                 # copy p-file (requires autolock=2)
get_archive = True                # download ScanArchive from vre

what = sor.get_what()
topdir = sor.dir(what)

# redirect output to logfile
if len(sys.argv)>6:
  now = datetime.datetime.now()
  fnlog = topdir + 'log/' + now.strftime("%Y%m%d_%H%M%S") + '_archive_all.txt'
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
if get_dicom:
  print('Copying dicom files')
  sor.copy_dicom(exam,series,datadir)

# copy p-file
if get_pfile:
  print('Copying p-file')
  rawdataname = sor.copy_pfile(runnum,datadir)
  print(rawdataname)
  if rawdataname is None:
    print('Warning: no rawdata found')

# download ScanArchive
if get_archive and sor.is_ox():
  print('Downloading ScanArchive')
  rawdataname = sor.download_scanarchive(exam,series,runnum,datadir)
  print(rawdataname)
  if rawdataname is None:
    print('Warning: no rawdata found')

# finishing up
print('archive_all.py done')
if len(sys.argv)>6:
  sys.stdout.close()

sys.exit(0)

