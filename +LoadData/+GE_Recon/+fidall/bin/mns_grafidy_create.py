#!/usr/bin/python

# 2/2020 Rolf Schulte
import sys
import os
import shutil
import re
# import hashlib
# hashlib not working because of stupic FIPS mode on scanner
import inspect


# function help
if len(sys.argv)<2:
  print('Insufficient input arguments: Usage')
  print('\tmns_grafidy_create.py nucnum')
  print('\tnucnum  Nucleus number: 0,2,3,13,23,31, or 129')
  print('\t        0 -> set B0 coefficients to zero')
  print('\texecute inside caldir')
  print('\tcd /export/home/signa/caldir')
  print('  Service key required for grafidy')
  sys.exit(-1)


# parameters
nucnum = int(sys.argv[1])
ratios = {0:0.0,2:0.1535,3:0.7618,13:0.2515,23:0.2645,31:0.4048,129:0.2786}
fnames = ["grafidy.cal","b0_ecc_coeff.dat","b0_ecc_coeff_ice.dat"]
nuclei = {0:'zero',2:'deuterium',3:'helium',13:'carbon',23:'sodium',31:'phosphorus',129:'xenon'}
md5fname = 'mns_grafidy_md5sum.txt'


# misc checks
nucname = nuclei.get(nucnum,'None')
if nucname=='None':
  print('Error: nucleus',nucnum,'not supported')
  sys.exit(-1)
ratio = ratios.get(nucnum,'None')
if ratio=='None':
  print('Error: nucleus',nucnum,'not supported')
  sys.exit(-1)
print('nucnum=',nucnum,'nucname=',nucname,'ratio=',ratio)

if os.getcwd()!='/export/home/signa/caldir':
  print('Warning: cwd (=', os.getcwd(), ') != /export/home/signa/caldir')
  tmp = raw_input('press return to continue')
  
# rename original proton files: grafidy.cal, b0_ecc_coeff.dat and b0_ecc_coeff_ice.dat
for fname in fnames:
  print('Checking',fname)
  if os.path.isfile('./'+fname):
    if not os.path.isfile('./'+fname+'.proton'):
      shutil.copyfile(fname,fname+'.proton')
  else:
    print('Error: file',fname,'not found')
    sys.exit(-1)

# create file with md5sum of proton files
if not os.path.isfile('./'+md5fname):
  print('Calculating md5sum of '+'.proton '.join(fnames)+'and store in '+md5fname)
  if 1:
    md5cmd = 'md5sum '+' '.join(fnames)+' > '+md5fname
    print('Executing system command:\n\t'+md5cmd)
    if os.system(md5cmd) != 0:
      print('Warning: executing md5sum failed')
  else:
    fout = open(md5fname,'w')
    for fname in fnames:
      fin = open(fname+'.proton', 'r')
      m = hashlib.md5()
      data = fin.read()
      m.update(data)
      fin.close()
      md5 = m.hexdigest()
      txt = md5+'  '+fname
      print(txt)
      fout.write(txt+'\n')
    fout.close()

# modify grafidy.cal for each nucleus
# look for lines with b0 -> scale amplitude (last entry) by gyromagnetic ratio
b0 = re.compile('b0')
print('Creating file= '+fnames[0]+'.'+nucname)
fin  = open(fnames[0]+'.proton',"r")
fout = open(fnames[0]+'.'+nucname,"w")
for line in fin:
  m = b0.search(line)
  if m:
    x = line.split()
    xx = float(x[-1])*ratio
    # print(x[-1],' ',xx)
    x[-1] = str(xx)
    newline = '  '+'    '.join(x)+'\n'
    fout.write(newline)
  else:
    fout.write(line)
fin.close()
fout.close()


# rename original grafidy.cal and create symbolic link to grafidy.cal.nucleus
for fname in fnames:
  os.rename(fname,fname+'.orig')
os.symlink(fnames[0]+'.'+nucname,fnames[0])

# instruct user to start and exit grafidy
print('Open and exit grafidy tool:')
print('-> Service Desktop Manager -> Service Browser -> Calibration')
print('   -> Calibration Tools -> Grafidy3')
print('Do NOT abort via ctrl-c (leaves system calibration in wrong state)')
tmp = raw_input('when done with grafidy, press enter to continue:')

# renaming new files
for l in range(1,3):
  fname = fnames[l]
  print('Renaming '+fname+' to '+fname+'.'+nucname)
  if os.path.isfile(fname):
    os.rename(fname,fname+'.'+nucname)
  else:
    print('Attention: file '+fname+' not found')
    print('-> translating '+fnames[0]+' into '+fname+' NOT successful')
    print('-> rerun mns_grafidy_create.py and grafidy tool')

# restoring original proton files
print('Restoring original proton calibration files')
os.unlink(fnames[0])
for fname in fnames:
  os.rename(fname+'.orig',fname)

# checking md5sum
def foo():
  pass

if os.system('mns_grafidy_check.py 1') != 0:
  chkcmd = os.path.abspath(inspect.getfile(foo))
  chkcmd = chkcmd[0:-9]+'check.py'
  print('mns_grafidy_check.py not in path; trying')
  print(chkcmd)
  if os.system(chkcmd+' 1') != 0:
    print('Attention: could not start mns_grafidy_check.py')
    print('md5sum not verified')


# main function exit code
sys.exit(0)

