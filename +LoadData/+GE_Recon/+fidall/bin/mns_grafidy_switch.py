#!/usr/bin/python

# 12/2019 Rolf Schulte
import sys
import os
import shutil

print('Input arguments', str(sys.argv))

# function help
if len(sys.argv)<2:
  print('Insufficient input arguments: Usage')
  print('\tmns_grafidy_switch.py nucleus')
  print('\tnucleus=0,1,2,3,13,23,31 or 129')
  print('\tCreate files b0_ecc_coeff.dat.nucleus b0_ecc_coeff_ice.dat.nucleus first')
  print('\t -> mns_grafidy_create.py')
  print('\tRequires proactive option key -> check via sok_test')
  sys.exit(-1)

# parameters
nucleus = int(sys.argv[1])
atp_script = '/usr/g/bin/update_ecc.atp'
atp_cmd = 'atp '+atp_script
caldir = '/usr/g/caldir/'
nuclei = {0:'zero',1:'proton',2:'deuterium',3:'helium',13:'carbon',23:'sodium',31:'phosphorus',129:'xenon'}
nucname = nuclei.get(nucleus)
fnames = ["grafidy.cal","b0_ecc_coeff.dat","b0_ecc_coeff_ice.dat"]


# check if md5sum file is existing
fnmd5 = caldir+'mns_grafidy_md5sum.txt'
print('Checking if md5checksum file existing: '+fnmd5)
if not os.path.isfile(fnmd5):
  print('Error: file not found: '+fnmd5)
  sys.exit(-1)

# check if proton files existing
for l in range(1,3):
  fn1h = caldir+fnames[l]+'.proton'
  print('Checking if file existing: '+fn1h)
  if not os.path.isfile(fn1h):
    print('Error: file ('+fn1h+') not found')
    print('  run mns_grafidy_create.py first')
    sys.exit(-1)

# check if nucleus files existing
if nucleus>1:
  for l in range(1,3):
    fnnuc  = caldir+fnames[l]+'.'+nucname
    print('Checking if file existing: '+fnnuc)
    if not os.path.isfile(fnnuc):
      print('Error: file ('+fnnuc+') not found')
      print('  run mns_grafidy_create.py first')
      sys.exit(-1)

# copying files
for l in range(1,3):
  fnnuc  = caldir+fnames[l]+'.'+nucname
  fnorig = caldir+fnames[l]
  print('Removing original file: '+fnorig)
  os.remove(fnorig)
  print('Copying nucleus file '+fnnuc+' to '+fnorig)
  shutil.copyfile(fnnuc,fnorig)

# create update_ecc.atp file if not existing
if os.path.isfile(atp_script):
  print('atp script found: '+atp_script)
else:
  print('Creating atp script: '+atp_script)
  fout = open(atp_script,'w')
  fout.write('DOWNLOAD_ECC;\n')
  fout.close()

# execute atp command  
if os.system(atp_cmd) != 0:
  print('Warning: running atp failed; starting tps_reset')
  if os.system('tps_reset') != 0:
    print('Warning: tps_reset failed')
    print('Changing grafidy to '+nucname+' failed')
  else:
    print('Changing grafidy to '+nucname+' successful (via tps_reset)')
else:
  print('Changing grafidy to '+nucname+' successful (via atp)')
  

# checking md5sum
def foo():
  pass

if nucleus==1:
  ossysval = os.system('mns_grafidy_check.py 1')
  if ossysval!=0:
    if ossysval==32512:
      print('Warning: mns_grafidy_check.py command not found')
      print('Add to binary path; eg via "source ~/fidall/cshrc"')
      print('Trying "/export/home/sdc/fidall/mns_grafidy_check.py 1"')
      if os.system('/export/home/sdc/fidall/bin/mns_grafidy_check.py 1') != 0:
        print('Warning: mns_grafidy_check.py failed')
    else:
      print('Warning: mns_grafidy_check.py failed')
else:
  print('Do NOT forget to switch back to protons!!!')


# main function exit code
sys.exit(0)

