#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Download ScanArchive file from vre
Created: 7/2020
@author: Rolf Schulte
"""

import sys
import sonofrecon as sor


# start main function
print(str(sys.argv))

# function help
if len(sys.argv)<3:
  print('Insufficient input arguments: Usage')
  print('\tget_h5_scan_archive_file.py <exam#> <series#> <runnum> <datadir> <fpos>')
  print('\tDownload ScanArchive from vre from directory')
  print('\t/data/arc/Closed/Exam<exam#>/Series<series#>/')
  print('\t<runnum>      check header for run number (optional)')
  print('\t<datadir>     target directory (default=current)')
  print('\t<fpos>==last  download last (matching) file (default)')
  print('\t<fpos>==all   download all (matching) files')
  print('\t<fpos>==first download first (matching) files')
  sys.exit(-1)

# input and default parameters
exam   = str(sys.argv[1])
series = str(sys.argv[2])
if len(sys.argv)>3:
  runnum = str(sys.argv[3])
else:
  runnum = '0'
if len(sys.argv)>4:
  datadir = sys.argv[4]
else:
  datadir = '.'
if len(sys.argv)>5:
  fpos = sys.argv[5]
else:
  fpos = 'last'

# start download
if sor.download_scanarchive(exam,series,runnum,datadir,fpos) is None:
  print('Error: downloading ScanArchive failed')
  sys.exit(-1)

# finishing up
print('Download finshed without error')

sys.exit(0)                  # main function exit code


