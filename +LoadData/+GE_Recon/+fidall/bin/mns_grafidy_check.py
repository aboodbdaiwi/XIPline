#!/usr/bin/python

# 11/2019 Rolf Schulte
import sys
import os
import shutil
# import hashlib
# hashlib not working because of stupic FIPS mode on scanner


# parameters
fnames = ["grafidy.cal","b0_ecc_coeff.dat","b0_ecc_coeff_ice.dat"]
md5fname = 'mns_grafidy_md5sum.txt'
caldir = '/usr/g/caldir/'
# caldir = './'


# function help
if len(sys.argv)<2:
  print('Insufficient input arguments: Usage')
  print('\tmns_grafidy_check.py anyinput')
  print('Function for checking md5sum of calibration files')
  print('\t'+' '.join(fnames))
  print('\tagainst values stored in '+md5fname+' ; path='+caldir)
  sys.exit(-1)


# checking if files exist
if not os.path.isfile(caldir+md5fname):
  print('Error: file not found: '+caldir+md5fname)
  sys.exit(-1)

for fname in fnames:
  if not os.path.isfile(caldir+fname):
    print('Error: file not found: '+caldir+fname)
    sys.exit(-1)


# checking md5sum
print('Checking md5sum of '+' '.join(fnames)+' files')
print('md5checksum file existing: '+caldir+md5fname)
fin = open(caldir+md5fname, 'r')
md5file = fin.read()
fin.close()

if 1:
  os.chdir(caldir)
  md5cmd = 'md5sum '+' '.join(fnames)
  print('Executing system command:\n\t'+md5cmd)
  p = os.popen(md5cmd)
  md5check = p.read()
else:
  md5check = ''
  for fname in fnames:
    fin = open(caldir+fname, 'r')
    m = hashlib.md5()
    data = fin.read()
    m.update(data)
    fin.close()
    md5 = m.hexdigest()
    md5check = md5check+md5+'  '+fname+'\n'

if md5file!=md5check:
  print('Attention: md5sum not matching')
  print('File=\n'+md5file)
  print('Check=\n'+md5check)
  print('You must restore proton calibrations manually from backup')
  sys.exit(-1)
else:
  print('Checksum test successful')
  print('-> original proton calibration files correct')


# main function exit code
sys.exit(0)
