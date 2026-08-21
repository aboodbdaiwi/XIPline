#!/usr/bin/python3

import sys
import os
import subprocess
import re
import datetime

# ****************************
# set parameters
VERBOSITY   = '4'
topdir = '/export/home/sdc/mrshim/'
MODEL_FILE  = topdir + 'sphmodel_ge.csv'
LIC         = topdir + 'dev.lic'
do_interpol = False
interp_fac = 2

if True:
  now = datetime.datetime.now()
  tstr = now.strftime("%Y%m%d_%H%M%S") + '_'
else:
  tstr = ''
ZIP_DATA    = 'b0map.zip'
TUNEUP_FILE = 'tuneup.csv'
MASK_FILE   = tstr + 'GRx.txt'
SOLN_FILE   = tstr + 'soln.json'
SPH_FILE    = tstr + 'sph.npy'


# ****************************
# print help
if len(sys.argv)<3:
  print('mrshim.py <exam> <series>')
  print('Calculate local shim coil setting')
  print('Prescribe MRS voxel and input exam and series number from B0map')
  sys.exit(-1)

# input paraters
print('Input arguments ' + str(sys.argv))
exam   = str(sys.argv[1])
series = str(sys.argv[2])


# ****************************
# create datadir
if not os.path.isdir(topdir + 'data'):
  os.mkdir(topdir + 'data')
datadir = topdir + 'data/Exam' + exam
if not os.path.isdir(datadir):
  os.mkdir(datadir)
datadir = topdir + 'data/Exam' + exam + '/Series' + series + '/'
if not os.path.isdir(datadir):
  os.mkdir(datadir)
print('datadir = ' + datadir)


# ****************************
# get dicom directory name via pathExtract
cmd = ['pathExtract',exam,series]
print(cmd)
sp = subprocess.Popen(cmd,shell=False,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
dcm = sp.stdout.readlines()
if sys.version_info[0]>2:    # python3 requires conversion to str
  dcm = [x.decode('ISO-8859-1') for x in dcm]
dcm.pop(0)                   # remove first line
DCM_FOLDER = str(re.sub('/i\d+.MRDC.\d+\n','/',dcm[0]))
print('DCM_FOLDER = ' + DCM_FOLDER)
if not DCM_FOLDER:
  print('Error: DCM_FOLDER not found in regular expression')
  sys.exit(-1)
if not os.path.isdir(DCM_FOLDER):
  print('Error: directory not found: DCM_FOLDER = ' + DCM_FOLDER)
  sys.exit(-1)


# ****************************
# create mask file
if True:
  cmd = ['printSHM']
  print(cmd)
  proc = subprocess.Popen(cmd,shell=False,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
  (stdoutdata, stderrdata) = proc.communicate()
  stdoutdata = stdoutdata.decode('ISO-8859-1')
  stderrdata = stderrdata.decode('ISO-8859-1')
  # print('stdoutdata =\n' + str(stdoutdata))
  # print('stderrdata =\n' + str(stderrdata))
  # sys.stdout.flush()
  if not str(stdoutdata):
    print('Error: stdoutdata')
    sys.exit(-1)
  if str(stderrdata):
    print('Error: stderrdata')
    print('stderrdata =\n' + str(stderrdata))
    sys.stdout.flush()
    sys.exit(-1)

  # check voxelExists
  for line in stdoutdata.splitlines():
    if re.search('voxelExists',line):
      print(line)
      voex = re.findall('\d+',line)[0]
      if voex=='0':
        print('Error: voxel not existing')
        sys.exit(-1)
      else:
        if voex!='1':
          print('Error: voxelExists ~= 0 or 1')
          sys.exit(-1)

  fid = open(datadir + MASK_FILE,'wt')
  if fid.closed:
    print('Error: cannot open ' + MASK_FILE)
    sys.exit(-1)
  fid.write(stdoutdata)
  fid.close()
        

# ****************************
# Convert dicom files to MRShim .zip format
cmd = [topdir+'bin/import_b0fields',DCM_FOLDER,'--output:'+datadir+ZIP_DATA,'--overwrite','--license:'+LIC,'-v:'+VERBOSITY]
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
  sys.stdout.flush()
  sys.exit(-1)
 
 
# ****************************
# Convert dicom files to MRShim .zip format
if do_interpol:
  cmd = [topdir+'bin/gen_interpolated',datadir+ZIP_DATA,'--model:resample','--resolution:'+str(interp_fac),'--output:'+datadir+ZIP_DATA,'--license:'+LIC,'--v:4']
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
    sys.stdout.flush()
    sys.exit(-1)


# ****************************
# Generate the mask and update the .zip file
#   @maskfile: File of the exported protocol parameters that have a selected shim volume
# cmd = [topdir+'gen_mask',datadir+ZIP_DATA,'--maskfile:'+datadir+MASK_FILE,'--license:'+LIC,'-v:'+VERBOSITY]
cmd = [topdir+'bin/gen_mask',datadir+ZIP_DATA,'--maskfile:'+datadir+MASK_FILE,'--license:'+LIC,'-v:'+VERBOSITY,\
  '--pos_offset:[-0.004,-0.004,-0.004]']
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
  sys.stdout.flush()
  sys.exit(-1)


# ****************************
# Generate the spherical harmonic fields for the co-ordinates of the input file
#   @prepend_f0: Ensures that the frequency offset field is included
#   @model: Model file of the spherical harmonic fields
#   @output: Output file (as .npy)
cmd = [topdir+'bin/gen_sphharm',datadir+ZIP_DATA,'--prepend_f0','--model:'+MODEL_FILE,\
  '--output:'+datadir+SPH_FILE,'--license:'+LIC,'-v:'+VERBOSITY]
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
  sys.stdout.flush()
  sys.exit(-1)


# ****************************
# Calculate the shim values
#   @sphharm: Degree of spherical harmonics to use in calculation
#   @sphfile: Spherical harmonic fields (output file from `gen_sphharm`)
#   @tuneup_file: File of tune-up values and sensitivities
#   @output: File to save the shim value solution in json format
cmd = [topdir+'bin/calc_shims',datadir+ZIP_DATA,'--sphharm:1','--sphfile:'+datadir+SPH_FILE,\
  '--tuneup_file:'+datadir+TUNEUP_FILE,'--output:'+datadir+SOLN_FILE,'--license:'+LIC,'-v:'+VERBOSITY]
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
  sys.stdout.flush()
  sys.exit(-1)

# extract linear shim values
shim_term = ['f0','Gx','Gy','Gz']
for line in stdoutdata.splitlines():
  if re.search('Sph vals',line):
    print(line)
    shims = re.findall('[0-9ed\-\.]+',line)
    # print(shims)
    if len(shims)!=4:
      print('Error: len(shims)!=4')
      sys.exit(-1)
    for ll in range(4):
      print(shim_term[ll] + ' = ' + str(round(float(shims[ll])*10e6)))
    

# ****************************
sys.exit(0)


