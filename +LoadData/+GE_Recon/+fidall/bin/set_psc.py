#!/usr/bin/python
# 3/2020 Rolf Schulte

import sys
import os

# start main function
print('--- set_psc.py ---')
print(str(sys.argv))


# function help
if len(sys.argv)<5:
  print('Setting prescan values via script')
  print('Requires proactive option key (check with sok_test) or service key')
  print('Insufficient input arguments: Usage')
  print('\tset_psc.py <R1> <R2> <TG> <freq>')
  print('\t<R1>     Analogue receiver gain: 1-13')
  print('\t<R2>     Digital receiver gain: 1-30 (w/ EDR)')
  print('\t<TG>     Transmit gain: 0-200 [dB/10]')
  print('\t<freq>   Centre frequency [Hz]')
  print('\tinput=-1 no update')
  sys.exit(-1)

# input and default parameters
r1   = int(sys.argv[1])
r2   = int(sys.argv[2])
tg   = int(sys.argv[3])
freq = int(sys.argv[4])

# check ranges
if r1<-1:
  print('Error: R1 (='+str(r1)+') <- 1')
  sys.exit(-1)
if r1==0:
  print('Error: R1 (='+str(r1)+') == 0')
  sys.exit(-1)
if r1>13:
  print('Error: R1 (='+str(r1)+') > 13')
  sys.exit(-1)
if r2<-1:
  print('Error: R2 (='+str(r2)+') <- 1')
  sys.exit(-1)
if r2==0:
  print('Error: R2 (='+str(r2)+') == 0')
  sys.exit(-1)
if r2>30:
  print('Error: R2 (='+str(r2)+') > 30')
  sys.exit(-1)
if tg<-1:
  print('Error: TG (='+str(tg)+')< -1')
  sys.exit(-1)
if tg>200:
  print('Error: TG (='+str(tg)+') > 200')
  sys.exit(-1)
if freq<-1:
  print('Error: freq (='+str(freq)+') <- 1')
  sys.exit(-1)
if freq>128000000:
  print('Error: frequency (='+str(freq)+'[Hz]) too high')
  sys.exit(-1)

# write atp script
atp_cmd = 'PSC_UPDATE_VAL \"'+str(r1)+' '+str(r2)+' '+str(tg)+' '+str(freq)+'\";'
fname = 'psc_tmp134.atp'
fid = open(fname,"w")
fid.write(atp_cmd)
fid.close()

# execute atp script
cmd = 'atp '+fname
print(cmd)
if os.system(cmd) != 0:
  print('Warning: running atp failed')


sys.exit(0)                  # main function exit code

