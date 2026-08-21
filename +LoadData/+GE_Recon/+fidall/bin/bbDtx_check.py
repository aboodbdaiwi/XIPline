#!/usr/bin/python
# 4/2023 Rolf Schulte

import sys
import os
#import subprocess
#import shutil
import glob
import re
#if os.path.isdir(npp):
#  sys.path.insert(0,npp)
try:
  import numpy as np
except:
  sys.path.append('/export/home/sdc/vxtl/python/')
  llp = os.environ.get('LD_LIBRARY_PATH','')
  if not 'vxtl' in llp:
    print('Add /export/home/sdc/vxtl/lib to LD_LIBRARY_PATH')
  import numpy as np
  
debug = 0

# ---------------------------------------------
# sub-function to read bbDtx calibration file
def read_bbDtx(fname,debug=0):
  data = []
  ln = -100000               # line number: "Coarse Attenuator" -> ln=0
  lnuc = 0                   # nuclei counter
  fid = open(fname,"rt")
  for line in fid:
    # print(line)
    # look for nucleus
    if 'Coarse Attenuator' in line:
      nuc = int(re.findall('(\d+)',line)[0])
      if debug>1:
        print(line)
        print(nuc)
      if nuc not in [2,3,7,13,19,23,31,129]:
        print('Warning: nuc(='+str(nuc)+') not in list')
      ln = 0                 # start line counting

    # extract signal level at step=0 (0dB)
    if ln==4:
      if debug>1:
        print(line)
      val = re.findall('([\d\.]+)',line)
      if debug>1:
        print(val)
      if val[0] != '0':
        print('Warning: val[0](='+val[0]+') != 0')
      data.insert(lnuc,[nuc,float(val[1])])
      # print(data)
      lnuc+=1
    ln+=1
  fid.close()
  if debug>1:
    print(data)
  return(data)


# ---------------------------------------------
# sub-function to read (old) dtx2 calibration file
def read_dtx2(fname,debug=0):
  data = []
  ln = -100000               # line number: "Original Calibration" -> ln=0
  x = re.search('dtx2_cal_',fname)
  if not x:
    print('Error: nuc not found')
  nuc = int(re.findall('(\d+)',fname[x.end():x.end()+3])[0])
  if debug>1:
    print(nuc)
  fid = open(fname,"rt")
  for line in fid:
    # print(line)
    if 'Original Calibration' in line:
      ln = 0                 # start line counting

    # extract signal level at step=0 (0dB)
    if ln==3:
      if debug>1:
        print(line)
      val = re.findall('(\S+)',line)
      if debug>1:
        print(val)
        print(len(val))
      if len(val) != 6:
        print('Warning: len(val)(='+str(len(val))+') != 6:')
      if val[0] != '00':
        print('Warning: val[0](='+val[0]+') != 00')
      data = [[nuc,float(val[4])]]
    ln+=1
  fid.close()
  if debug>1:
    print(data)
  return(data)


# ---------------------------------------------
# executed code
# function help
if len(sys.argv)<2:
  print('Error: at least one input argument required')
  print('Load, quantify and check bbDtx calibration files')
  print('Usage:\n\tbbDtx_check.py dir')
  print('\tdir    Directory name')
  print('\t       default=/export/home/service/log')
  sys.exit(-1)


# input parameters
datadir = str(sys.argv[1])
if not os.path.isdir(datadir):
  print('Error: directory (=\''+datadir+'\') not found')
  sys.exit(-1)


# look for files
new = True
fnames = glob.glob(datadir + '/bbDtx_1_cal*')
# print(len(fnames))
if len(fnames)==0:
  fnames = glob.glob(datadir + '/dtx2_cal*')
  if len(fnames)==0:
    print('Error: no bbDtx_1_cal or dtx2_cal files found')
    sys.exit(-1)
  else:
    new = False

# read data
nucs = []
data = []
for fn in fnames:
  if debug>0:
    print(fn)
  if new:
    d1 = read_bbDtx(fn)
  else:
    d1 = read_dtx2(fn)
  if debug>0:
    print(d1)
  for d2 in d1:
    if d2[0] not in nucs:
      nucs.append(d2[0])
      data.append([])
    inuc = nucs.index(d2[0])
    data[inuc].append(d2[1])

# sort data
ind = [nucs.index(x) for x in sorted(nucs)]
tmp = []
for l in ind:
  tmp.append(data[l])
data = tmp
nucs.sort()
nn = len(nucs)

# print data
print('---------------------------------------------')
print('nucs = '+str(nucs))
for lnuc in range(nn):
  print('% nuc = ' + str(nucs[lnuc]))
  print('data(:,' + str(lnuc+1) + ') = ...')
  print(data[lnuc])

# print info
i31 = nucs.index(31)
m31 = np.mean(data[i31])
print('---------------------------------------------')
print('nucleus\tsignal\t\tmin(sig)\tmax(sig)\tsig/sig_31P\tmax/min')
print('\t[SigLvl]\t[SigLvl]\t[SigLvl]\t[dB]\t\t[dB]')
for lnuc in range(nn):
  md = np.mean(data[lnuc])
  mind = np.min(data[lnuc])
  maxd = np.max(data[lnuc])
  print('%d \t%.0f+-%.2f\t%.0f\t\t%.0f\t\t%.2f\t\t%.3f' % \
    (nucs[lnuc], md, np.std(data[lnuc],ddof=1), mind, maxd, \
    10*np.log10(md/m31),10*np.log10(maxd/mind)))


# main function exit code
sys.exit(0)

