#!/bin/tcsh

set matlabroot=/b/apps/linux/matlab/matlab-2022a/

echo "Starting mcc R2022a 9.12"
echo "matlabroot="
echo "$matlabroot"

if ( -f matlab_scripter.ctf ) then
  echo "matlab_scripter.ctf existing; moving to old"
  rm old/matlab_scripter
  rm old/matlab_scripter.ctf
  mv matlab_scripter matlab_scripter.ctf old/
else
  echo "matlab_scripter.ctf not existing"
endif


nohup $matlabroot/inst/bin/mcc -o matlab_scripter -Cmv matlab_scripter.m \
-I ../matlab \
-I ../matlab_private \
-I ../read_mr \
-I ../Ox3.0.3 \
-I ../OXSA/main \
-I ../OXSA/utils \
-a ../read_mr/hdr_sizes.txt

