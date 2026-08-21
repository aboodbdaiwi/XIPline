#!/usr/bin/python
# 11/2025 Rolf Schulte

import sys
import os
import subprocess
import shutil
import re
import time

stars = '*******************************************************'
dcm_input_list = True

# function help
if len(sys.argv)<2:
  print('Error: at least one input argument required')
  print('Create anonymised copy <dir>anon of directory <dir>')
  print('Usage:\n\tanon_dir.py dir')
  print('\tdir    Directory name')
  sys.exit(-1)


# input parameters
dir = str(sys.argv[1])


# check if input is exam number
isInt = True
try:
   # converting to integer
   int(dir)
except ValueError:
   isInt = False
if isInt:
   dir = 'Exam'+str(dir)

# convert relative to absolute paths
print(dir)
dir = os.path.abspath(dir)
print(dir)

# check if directories exist
print(stars)
print('Anonymising directory \''+dir+'\'')
print(stars)
if not os.path.isdir(dir):
  print('Error: directory not found')
  sys.exit(-1)
if os.path.isdir(dir+'anon'):
  print('Error: directory existing: \''+dir+'anon\'')
  sys.exit(-1)


# copy everything
shutil.copytree(dir,dir+'anon')
time.sleep(0.1)

if dcm_input_list:
  fndcm = []

# anonymise files
err = False
fnnas = []         # file names not anonymised
for path, subdirs, files in os.walk(dir+'anon'):
  print(path)
  print(subdirs)
  print(files)
  for name in files:
    cmd = ''
    fname = os.path.join(path, name)
    
    # p-files
    if re.search('P\d{5}.7',name):
      cmd = ['pfile_anon',fname]
      fnanon = fname + '.anon'
    
    # ScanArchives
    if name.startswith('ScanArchive') and name.endswith('.h5'):
      fnanon = fname + '.anon'
      cmd = ['ArchiveTool','--anon-header','--anon-mode','Partial',\
        '--input-file',fname,'--output-file',fnanon]
        
    # mat files
    if name.endswith('.mat'):
      if name.startswith('ScanArchive') or re.search('P\d{5}',name):
        cmd = ['anonymise_mat.py',fname]
        fnanon = fname[0:-4] + '_anon.mat'
      else:
        print(stars)
        print(name + ' not a fidall data file')
        
    # dicom files
    if re.search('MRDC.\d+$',name) or name.endswith('.dcm'):
      if dcm_input_list:
        fndcm.append(fname)
      else:
        cmd = ['anonymise_dicom.py',fname]
    
    # header files
    if name.startswith('ScanArchive') and name.endswith('.hdr'):
      print(stars)
      print('removing ' + fname)
      os.remove(fname)

    # execute cmd
    if cmd:
      print(stars)
      print(cmd)
      # subprocess.call(cmd)
      proc = subprocess.Popen(cmd,shell=False,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
      (stdoutdata, stderrdata) = proc.communicate()
      stdoutdata = stdoutdata.decode('ISO-8859-1')
      stderrdata = stderrdata.decode('ISO-8859-1')
      print('stdoutdata =\n' + str(stdoutdata))
      sys.stdout.flush()
      if str(stderrdata):
        print('Error: stderrdata')
        print('stderrdata =\n' + str(stderrdata))
        err = True      
      
      if os.path.isfile(fnanon):
        os.remove(fname)
        shutil.move(fnanon,fname)
      else:
        print('Warning: File not found: \'' + fnanon + '\'')
        err = True
        fnnas.append(fname)

# anonymise dicom files
if dcm_input_list:
  cmd = fndcm
  cmd.insert(0,'anonymise_dicom.py')
  # cmd = ['anonymise_dicom.py',fndcm[0]]
  print(stars)
  print(cmd)
  proc = subprocess.Popen(cmd,shell=False,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
  (stdoutdata, stderrdata) = proc.communicate()
  stdoutdata = stdoutdata.decode('ISO-8859-1')
  stderrdata = stderrdata.decode('ISO-8859-1')
  print('stdoutdata =\n' + str(stdoutdata))
  sys.stdout.flush()
  if str(stderrdata):
    print('Error: stderrdata')
    print('stderrdata =\n' + str(stderrdata))
    err = True

# print info
print(stars)
if err:
  print('Error: anonymising data failed.')
  print('Not anonymised files:')
  print(fnnas)
  sys.exit(-1)
else:
  print('Anonymising data OK.')

# main function exit code
sys.exit(0)

