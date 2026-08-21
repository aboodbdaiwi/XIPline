#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Download p-files and/or ScanArchive files to fidall data
Created: 10/2025
@author: Rolf Schulte
"""

import sys
import sonofrecon as sor
import datetime


# set parameters
get_pfile = False                 # copy p-file
get_image_archive = True          # download ScanArchive from vre
fidalldir = sor.dir('fidall')

# redirect output to logfile
if len(sys.argv)>6:
  now = datetime.datetime.now()
  fnlog = fidalldir + 'log/' + now.strftime("%Y%m%d_%H%M%S") + '_archive_pfile.txt'
  sys.stdout = open(fnlog,'w')

# print help
if sor.print_help(sys.argv):
  if len(sys.argv)>6:
    sys.stdout.close()
  sys.exit(-1)

# get input parameters
(runnum,exam,series,functionname) = sor.function_input(sys.argv)

# define data directory
datadir = sor.datadir(exam,series)

# copy rawdata
if get_pfile:
  print('Copying p-file')
  rawdataname = sor.copy_pfile(runnum,datadir)
  print(rawdataname)
  if rawdataname is None:
    print('Warning: no rawdata found')

# download ScanArchive
if get_image_archive and sor.is_ox():
  print('Downloading ScanArchive')
  rawdataname = sor.download_scanarchive(exam,series,runnum,datadir)
  print(rawdataname)
  if rawdataname is None:
    print('Warning: no rawdata found')

# finishing up
print('archive_pfile.py done')
if len(sys.argv)>6:
  sys.stdout.close()

sys.exit(0)

