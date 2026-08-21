#!/usr/bin/python3
# 7/2025 Rolf Schulte

import sys
import os
import glob
import sonofrecon as sor
import datetime
import pathlib

# set parameters
fidalldir = sor.dir('fidall')
# fidalldir = '/home/schulter/projects/mnsrp/'
dirs = [fidalldir + 'waveforms/*', fidalldir + 'parameters/*', '/usr/g/bin/recon2*', fidalldir + 'bin/*']
stars = '***************************************************************************\n'
md5chk = True
cmds = ['getver','checkHardware','gemshid']

now = datetime.datetime.now()
fninfo = fidalldir + 'log/' + now.strftime("%Y%m%d_%H%M%S") + '_mnsrp_info.txt'
fid = open(fninfo, 'w')

# print function help
if True:
  print('Log mnsrp info to ' + fninfo)
  print('Files, symbolic links, md5sum in:\n\t' + '\n\t'.join(dirs))
  print('System commands:\n\t' + '\n\t'.join(cmds))

# get info from directories (symbolic links + md5sum)
for dname in dirs:
  print(stars + dname + '\n' + stars)
  fid.write(stars + dname + '\n' + stars)
  for fname in sorted(glob.glob(dname)):
    if md5chk:
      if os.path.isdir(fname):
        md5txt = '/'
      else:
        if pathlib.Path(fname).is_fifo():
          md5txt = '|'
        else:
          md5cmd = 'md5sum ' + fname
          p = os.popen(md5cmd)
          md5txt = p.read()
          md5txt = ' ' + md5txt[0:32]
    else:
      md5txt = ''

    if os.path.islink(fname):
      txt = os.path.basename(fname) + ' -> ' + os.readlink(fname) + md5txt
    else:
      txt =  os.path.basename(fname) + md5txt
    print(txt)
    fid.write(txt+'\n')
    
# get systeminfo
for cmd in cmds:
  print(stars + cmd + '\n' + stars)
  fid.write(stars + cmd + '\n' + stars)
  p = os.popen(cmd)
  txt = p.read()
  print(txt)
  fid.write(txt)

# print file dir+name
txt = stars + 'Info written to\n'+ fninfo + '\n'
print(txt)
fid.write(txt)

# close info file again
fid.close()

# main function exit code
sys.exit(0)

