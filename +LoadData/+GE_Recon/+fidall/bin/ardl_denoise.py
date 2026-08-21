#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Interface function for ARDL Denoising of fidall data
Created: 2/2025
@author: Rolf Schulte
"""

import sys
import sonofrecon as sor

# set parameters
fidalldir = sor.dir('fidall')
noise_df = '0.75'                      # denoising level:  0=off, 1=max
ring_df  = '1'                         # de-ringing level: 0=off, 1=max
inf_cmd  = '[]'                        # command for python inference
ardl_dir = fidalldir + 'ardl/'         # directory with REALInference.py + model
model = 'Standard_1562_July22_2019/'   # tensorflow model name
# py_cmd   = 'ssh vre /usr/bin/sudo /export/home/signa/vre/reset/vre/vreSudoWrapper.sh ' + \
#   '/usr/bin/docker run -v ' + ardl_dir + ':/ardl mnsrp_ardl:01 python3'
py_cmd = ''                            # command to start python
plt = 'true'                           # plotting

# print help
if len(sys.argv)<2:
  print('ardl_denoise.py fname')
  sys.exit(-1)

# get input parameters
fname = str(sys.argv[1])

# start matlab reconstruction on host
mascri = sor.setup_mcr()        # setup mcr
if mascri=='':
  print('Error: mascri not found')
  sys.exit(-1)
sor.print_env_vars()

# actual execution of matlab-compiled program
cmd = ['ardl_denoise',fname,'[]',noise_df,ring_df,py_cmd,inf_cmd,ardl_dir,model,plt]
sor.run(cmd,mascri)

sys.exit(0)

