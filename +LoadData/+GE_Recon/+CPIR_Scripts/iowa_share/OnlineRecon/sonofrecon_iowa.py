# -*- coding: utf-8 -*-
"""
Python functions for son-of-host reconstruction
11/2025
@author: Rolf Schulte

"""

import os
import sys
import shutil
import re
import subprocess
import time


# ******************************************
# default settings
# adjustment for recon delay (default='0')
def delay():
  delay = '0'
  return delay

# print more extensive debugging info (default=False)
def print_debug():
  return False

# start mascri instance or run in background (default=True)
def dofifo():
  return True

# display images via sor.view with eog/display/gimp (default=True)
def doview():
  return True
  
# apply gradwarp distortion correction (default='true')
def grdwrp():
  grdwrp = 'true'
  # grdwrp = 'false'
  return grdwrp

# copy noise scan from /usr/g/recon_calib/ExamData to datadir
def docopynoise():
  return True
  

# ******************************************
# print function help: inputs from son-of-host recon
def print_help(argv):
  print(str(argv))
  print(' '.join(str(x) for x in argv))
  if len(argv)<6:
    print('Error: at least five input arguments required')
    print('<pfile_no> <other> <scan_no> <exam_no> <series_no> <last_image>')
    print('\t<pfile_no>  Number of p-file (w/o P and .7)')
    print('\t<other>     Input ignored (use 0)')
    print('\t<scan_no>   Input ignored (use 0)')
    print('\t<exam_no>   Exam number')
    print('\t<series_no> Series number')
    print('\t<last image> Input irrelevant; if given, print output to text log file, otherwise to screen')
    sys.stdout.flush()
    return True
  else:
    return False


# ******************************************
# define directories
def dir(what):
  if what=='fidall':
    dname = '/export/home/sdc/fidall/'
  else:
    if what=='mrraw':
      dname = '/usr/g/mrraw/'
    else:
      if what=='meti':
        dname = '/export/home/sdc/meti/'
      else:
        print('Warning: dir(\''+what+'\') not found')
        dname = ''
  return dname


# ******************************************
# ScanArchives directory on vre
def vredatadir(exam,series):
  datadir = '/data/arc/Closed/Exam' + exam + '/Series' + series + '/'
  return datadir


# ******************************************
# data directory
def datadir(exam='000',series='000',what='fidall'):
  datadir = dir(what) + 'data/Exam' + exam + '/'
  if not os.path.isdir(datadir):
    os.mkdir(datadir)
  datadir = datadir + 'Series' + series + '/'
  if not os.path.isdir(datadir):
    os.mkdir(datadir)
  if os.environ.get('MASCRILOGDIR') is None:
    os.environ['MASCRILOGDIR'] = datadir
  return datadir


# ******************************************
# temporary data directory (mounted on vre as well)
def datadir_tmp(exam='000',series='000',what='fidall'):
  datadir = dir(what) + 'ex' + exam + '_se' + series + '/'
  if os.path.isdir(datadir):
    print('Warning: datadir (\'' + datadir + '\') already existing')
  else:
    os.mkdir(datadir)
  os.environ['MASCRILOGDIR'] = datadir
  return datadir


# ******************************************
# determine if fidall or meti
def get_what():
  what = 'fidall'       # default is fidall
  cwdir = os.path.dirname(os.path.abspath(__file__))
  if re.search('meti',cwdir) is not None:
    what = 'meti'
  print(what)
  
  return what


# ******************************************
# function input
def function_input(argv):
  runnum = str(argv[1])
  exam = str(argv[4])
  series = str(argv[5])
  functionname = str(argv[0])
  return (runnum,exam,series,functionname)


# ******************************************
# check if scanner has Orchestra (>=DV26)
def is_ox():
  prgname = 'ArchiveTool'
  import distutils.spawn
  return distutils.spawn.find_executable(prgname) is not None


# ******************************************
# copy p-file
def copy_pfile(runnum,datadir):
  import datetime

  num_str = runnum
  pfname = 'P' + num_str.zfill(5) + '.7'
  if not os.path.isfile(dir('mrraw')+pfname):
    xmessage('Error: pfile (\'' + pfname + '\') not found')
    return None

  # add date and time string to front of p-file name
  now = datetime.datetime.now()
  pftname = now.strftime("%Y%m%d_%H%M%S_") + pfname

  # copy p-file
  shutil.copyfile(dir('mrraw')+pfname,datadir+pftname)
  print(pfname)

  sys.stdout.flush()
  return pftname


# ******************************************
# delete p-file from mrraw dir
# mrraw is not cleaned up; plus autolock>0 required in fidall to prevent Ox error
def remove_pfile(runnum):
  num_str = runnum
  pfname = 'P' + num_str.zfill(5) + '.7'
  if not os.path.isfile(dir('mrraw')+pfname):
    # xmessage('Error: pfile (\'' + pfname + '\') not found')
    print('Warning: deleting p-file: \'' + pfname + '\' not found')
    return None

  # delete p-file
  print('Deleting p-file ' + dir('mrraw')+pfname)
  os.remove(dir('mrraw') + pfname)

  sys.stdout.flush()
  return pfname


# ******************************************
# copy dicom files
def copy_dicom(exam,series,datadir):
  # check/create data (sub)dir
  if not os.path.isdir(datadir):
    print('Datadir (='+datadir+') not existing: creating')
    os.mkdir(datadir)
  
  # get dicom file names via pathExtract
  cmd = ['pathExtract',exam,series]
  print(cmd)
  sp = subprocess.Popen(cmd,shell=False,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
  dcm = sp.stdout.readlines()
  if sys.version_info[0]>2:  # python3 requires conversion to str
    dcm = [x.decode('ISO-8859-1') for x in dcm]
  dcm.pop(0)              # remove first line

  # copy files
  print(dcm)
  if not dcm:
    print('Warning: no dicom files found')
  else:
    print('Copying ' + str(len(dcm)) + ' dicom files')
    for fname in dcm:
      fname = fname.rstrip()
      # print(fname)
      if os.path.isfile(fname):
        shutil.copy2(fname,datadir)
      else:
        print('Warning: file not found \'' + fname + '\'')


# ******************************************
# extract waveformname from symbolic link
def get_wfname(functionname,copy2dir=''):

  m = re.search('\d+$',functionname)
  rhrecon = m.group(0)
  print('rhrecon = ' + rhrecon)
  if int(rhrecon)<5000:
    traj = '%02d' % (int(rhrecon)-2600)
  else:
    traj = '%03d' % (int(rhrecon)-5000)
  ak_name = dir('fidall') + 'waveforms/'+'ak_grad'+traj+'.wav'
  print('Looking for symbolic link from: '+ak_name)

  link_name = os.readlink(ak_name)          # where symbolic link points to
  if link_name is None:
    xmessage('Error: symbolic link not found')
    return None

  link_name = link_name[0:-4] + '.mat'      # replace suffix
  # if tmp starts with '/' it is an absolute path
  if link_name[0]=='/':
    wfn = link_name
  else:
    wfn = dir('fidall') + 'waveforms/' + link_name
  print('wfn = ' + wfn)
  if not os.path.isfile(wfn):
    xmessage('Warning wfn (=' + wfn + ') not found')

  # copy trajectory mat file if copy2dir given
  dst = copy2dir + os.path.split(wfn)[1]
  if os.path.isdir(copy2dir) and os.path.isfile(wfn) and not os.path.isfile(dst):
    print('Copying waveform to ' + copy2dir)
    shutil.copy2(wfn,copy2dir)
  
  sys.stdout.flush()
  return wfn


# ******************************************
# setup variables for matlab runtime environment
def setup_mcr(versions = [],what='fidall'):
  if not versions:
    versions = [242,912,98,95,81]
  if os.path.isfile('startup'):
    print('Warning: startup file existing; changing directory to\n\t'+dir(what))
    sys.stdout.flush()
    os.chdir(dir(what))    # go to different directory 
    # (if started in /usr/g/bin matlab_scripter conflicts with startup file)
  
  mcr_root = ''    # MCR root directory
  mascri = ''      # matlab_scripter bin
  hostname = os.uname()[1]
  hostname = re.findall('\w+',hostname)[0]
  print('hostname: '+hostname)

  hostname = os.uname()[1]
  hostname 

  # list of list of potential MCR directories
  dnames = []
  for i in range(1000):
    dnames.append(0)
  dnames[81] = ['/usr/local/MATLAB/MATLAB_Compiler_Runtime/v81/',\
                '/export/home/signa/research/MATLAB/v81/']
  dnames[95] = ['/opt/mathworks_matlab_runtime_r2018b/root/v95/',\
                '/export/home/signa/research/MATLAB/v95/']
  dnames[98] = ['/opt/mathworks_matlab_runtime_r2020a/root/v98/',\
		'/export/home/signa/research/MATLAB/v98/']
  dnames[912] = ['/opt/mathworks_matlab_runtime_r2022a/root/v912/',\
		'/export/home/signa/research/MATLAB/v912/']
  dnames[242] = ['/opt/mathworks_matlab_runtime_r2024b/root/R2024b/',\
                 '/export/home/signa/research/MATLAB/v242/']

  # try is mascri and mcr exist
  for ver in versions:
    print('Trying mcr version = ' + str(ver))
    mascri1 = dir(what) + 'bin/mascri' + str(ver) + '_iowa/matlab_scripter'
    print(mascri1)
    if os.path.isfile(mascri1):
      print('Found mascri = ' + mascri1)
      for dn in dnames[ver]:
        if os.path.isdir(dn):
          print('Found dn = ' + dn)
          mcr_root = dn
          mascri = mascri1
          break
    if not mascri=='':
      break

  # print error if nothing found
  if mascri=='':
    xmessage('Error: no mascri found')
  if mcr_root=='':
    xmessage('Error: no mcr_root found')

  # variables
  mcr_arch = 'glnxa64'
  mcr_proc = 'amd64'
  mcr_java = mcr_root + 'sys/java/jre/' + mcr_arch + '/jre/lib/' + mcr_proc + '/'

  # environment variables
  os.environ['PATH'] = dir(what) + 'bin:' +  os.environ['PATH']

  LD_LIBRARY_PATH = mcr_root + 'runtime/' + mcr_arch + ':' + \
    mcr_root + 'bin/' + mcr_arch + ':' + mcr_root + 'sys/os/' + mcr_arch + ':' + \
    mcr_java + 'native_threads:' + mcr_java + 'server:' + mcr_java + 'client:' + \
    mcr_java[0:-1]
  if hostname=='vre':   # bug fix for missing libraries on MR29.0_R01
    LD_LIBRARY_PATH = LD_LIBRARY_PATH+':/export/home/signa/research/MATLAB/lib64'
  if os.environ.get('LD_LIBRARY_PATH') is None:
    os.environ['LD_LIBRARY_PATH'] = LD_LIBRARY_PATH
  else:
    os.environ['LD_LIBRARY_PATH'] = LD_LIBRARY_PATH + ':' + os.environ.get('LD_LIBRARY_PATH','')

  # print(os.environ.get('LD_LIBRARY_PATH'))

  if os.environ.get('XAPPLRESDIR') is None:
    os.environ['XAPPLRESDIR'] = mcr_root + 'X11/app-defaults'
  else:
    os.environ['XAPPLRESDIR'] = mcr_root + 'X11/app-defaults:' + os.environ.get('XAPPLRESDIR','')
  if os.environ.get('DCMDICT') is None:
    os.environ['DCMDICT'] = '/export/home/sdc/app-defaults/dicom/gems-dicom-dict.txt';
  if os.environ.get('homedir') is None:
    os.environ['homedir'] = dir(what)

  # workaround for inverted MNS Rx freq bug
  # os.environ['mnsrpinvmns'] = '-1'
  
  # set default data reading to Ox
  # os.environ['mnsrpuseox'] = '1'
  
 
  sys.stdout.flush()
  return mascri


# ******************************************
# copy matlab runtime libraries to directory mounted on vre
def copy_mcr(versions = [242,912,98,95,81]):
  mascri = ''      # matlab_scripter bin

  # list of list of potential MCR directories
  dnames = []
  for i in range(1000):
    dnames.append(0)
  dnames[81] = '/usr/local/MATLAB/MATLAB_Compiler_Runtime/v81/'
  dnames[95] = '/opt/mathworks_matlab_runtime_r2018b/root/v95/'
  dnames[98] = '/opt/mathworks_matlab_runtime_r2020a/root/v98/'
  dnames[912] = '/opt/mathworks_matlab_runtime_r2022a/root/v912/'
  dnames[242] = '/opt/mathworks_matlab_runtime_r2024b/root/R2024b/'
  mcr_copy = '/export/home/signa/research/MATLAB/'
  if not os.path.isdir(mcr_copy):
    print('Creating directory mcr_copy = \n\t'+mcr_copy)
    os.mkdir(mcr_copy)

  # check if mascri and mcr exist
  for ver in versions:
    print('Trying mcr version = ' + str(ver))
    mascri = dir('fidall') + 'bin/mascri' + str(ver) + '_iowa/matlab_scripter'
    print(mascri)
    if os.path.isfile(mascri):
      print('Found mascri = ' + mascri)
      if os.path.isdir(dnames[ver]):
        print('Found mcr = ' + dnames[ver])
        dncpy = mcr_copy+'v'+str(ver)
        if not os.path.isdir(dncpy):
          print('Copying mcr from\n\t'+dnames[ver]+'\n\tto\n\t'+dncpy)
          shutil.copytree(dnames[ver],dncpy,symlinks=True)
          time.sleep(5)      # dir not immediately found on vre
        else:
          print('Found dncpy = '+dncpy)
        if ver==95:
          dncpy = mcr_copy+'lib64'
          if not os.path.isdir(dncpy):
            print('Copying /usr/lib64 to '+'\n\tto\n\t'+dncpy)
            os.mkdir(dncpy)
            for ff in os.listdir('/usr/lib64'):
              f = '/usr/lib64/'+ff
              if not os.path.isdir(f):
                if os.path.islink(f):
                  linkto = os.readlink(f)
                  os.symlink(linkto,dncpy+'/'+ff)
                else:
                  try:
                    shutil.copy(f,dncpy)
                  except Exception:
                    pass

            # shutil.copytree('/usr/lib64',dncpy,symlinks=True)
            time.sleep(5)    # dir not immediately found on vre
      else:
        print('Warning: did not find mcr = '+dnames[ver])
      
  # print error if nothing found
  if mascri=='':
    xmessage('Error: no mascri found')
 

# ******************************************
# print environment variables
def print_env_vars():
  # print environment variables
  print('Environment variables')
  print('PATH=\n' + os.environ.get('PATH',''))
  print('LD_LIBRARY_PATH=\n' + os.environ.get('LD_LIBRARY_PATH',''))
  print('XAPPLRESDIR=\n' + os.environ.get('XAPPLRESDIR',''))
  print('DCMDICT=\n' + os.environ.get('DCMDICT',''))
  print('homedir=\n' + os.environ.get('homedir',''))
  print('IMPORT_IMAGE_DIR=\n' + os.environ.get('IMPORT_IMAGE_DIR',''))
  print('DISPLAY=\n' + os.environ.get('DISPLAY',''))
  sys.stdout.flush()


# ******************************************
# download ScanArchive from vre from directory:
# /data/arc/Closed/Exam<exam#>/Series<series#>
# runnum        check header for run number (optional)
# datadir       target directory (default=current)
# fpos=='last'  download last (matching) file (default)
# fpos=='all'   download all (matching) files
# fpos=='first' download first (matching) file
def download_scanarchive(exam,series,runnum='0',datadir='.',fpos='last'):
  import time

  # remote directory name
  rmtdir = vredatadir(exam,series)
  fnames = []
  contn = False

  print('Remote directory = ' + rmtdir)
  print('runnum = '+ runnum)

  # check if datadir is a directory
  if not os.path.isdir(datadir):
    xmessage('Error: datadir (=\'' + datadir + '\') not a directory')
    return None

  # check if remote directory exists
  # keep looking for 180secs
  for t in range(180):
    print(t)
    if not ls_rmtdir(rmtdir) == []:
      break
    time.sleep(1)
  else:
    xmessage('Error: remote directory not existing')
    return None

  # get filenames for downloading
  time.sleep(1)
  if runnum=='0':
    # time.sleep(1)
    fnames = ls_rmtdir(rmtdir)
    for fname in fnames:
      print('fname = '+fname)
  else:
    # look for ScanArchives with matching runnum
    # if nothing found: continue looking for 60secs
    for t in range(60):
      print(t)
      fns = ls_rmtdir(rmtdir)
      for fname in reversed(fns):
        print('Checking '+fname)
        chk = checkVREScanArchiveFileTag(rmtdir+fname,'rdb_hdr_image.rawrunnum',runnum)
        if chk:
          print(' -> True')
          fnames.append(fname)
          contn = True
          break 
        else:
          print(' -> False')
      if contn:
        break;
      time.sleep(1)
    else:
      xmessage('Error: no ScanArchive with runnum=' + runnum + ' found')
      return None

  # select file positions: last (default), all, first
  if fpos=='last':
    del fnames[0:-1]
  else:
    if fpos=='first':
      del fnames[1:]
    else:
      if fpos!='all':
        print('Warning: unkown fpos; downloading all found files')
   
  if len(fnames)==0:
    xmessage('Error: no files for download')
    return None

  # actual download
  print('Downloading '+str(len(fnames))+' files: '+' '.join(fnames))
  for fname in fnames:
    print('\t'+fname)
    rcp(rmtdir+fname,datadir)

  sys.stdout.flush()
  time.sleep(0.2)
  return fname


# ******************************************
# check ScanArchive tags
def checkVREScanArchiveFileTag(arcFile,tagName,tagValue):
  cmd = '/usr/g/bin/ArchiveTool --get '+tagName+' --input-file '+arcFile+' | tail -1'
  print(cmd)
  try:
    rsh = subprocess.Popen(['rsh', 'vre', cmd],shell=False,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
  except:
    rsh = subprocess.Popen(['ssh', '-oStrictHostKeyChecking=no', 'vre', cmd],shell=False,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
  result = rsh.stdout.readlines()
  print(result)
  if result == []:
    print('ERROR: tag %s not found in %s' % tagName, arcFile)
    return None
  else:
    string = str(result[0]).strip()
    pattern = '^\s*\S+\s*\(\d+\)\s*=\s*('+str(tagValue)+').*'
    ans = re.search(pattern,string)
    if ans:
      return True
    else:
      return False


# ******************************************
# list remote directory
def ls_rmtdir(rmtdir):
  cmd = 'ls '+rmtdir
  print(cmd)
  try:
    rsh = subprocess.Popen(['rsh', 'vre', cmd],shell=False,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
  except:
    rsh = subprocess.Popen(['ssh', '-oStrictHostKeyChecking=no', 'vre', cmd],shell=False,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
  result = rsh.stdout.readlines()
  if sys.version_info[0]>2:  # python3 requires conversion to str
    result = [x.decode('ISO-8859-1') for x in result]
  print(result)
  result = [x.rstrip("\n") for x in result]
  return result


# ******************************************
# remote copy
def rcp(fname,datadir):
  print('rcp/scp vre:'+fname+' '+datadir)
  try:
    rsh = subprocess.Popen(['rcp', 'vre:'+fname, datadir],shell=False,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
  except:
    # rsh = subprocess.Popen(['scp', 'vre:'+fname, datadir],shell=False,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    rsh = subprocess.Popen(['scp', '-oStrictHostKeyChecking=no', 'vre:'+fname, datadir],shell=False,stdout=subprocess.PIPE,stderr=subprocess.PIPE)


# ******************************************
# move dicom files to import directory: 
# all files/directories in datadir with suffix .dcm or _dcm
def import_dcm(datadir,copyfiles=False):
  import_image_dir = os.environ.get('IMPORT_IMAGE_DIR','/export/home1/sdc_image_pool/import')
  if not datadir.endswith('/'):
    datadir = datadir+'/'
  print('Importing dicom images from datadir \''+datadir+'\'')
  if copyfiles:
    print('by copying files ending with dcm or sdcopen to import dir \''+import_image_dir+'\'')
  else:
    print('by moving files ending with dcm or sdcopen to import dir \''+import_image_dir+'\'')
  sys.stdout.flush()

  if not os.path.isdir(import_image_dir):
    print('Warning: import_image_dir (\'' + import_image_dir + '\') not found')
    sys.stdout.flush()
    return False
  
  if not os.path.isdir(datadir):
    print('Warning: datadir (\'' + datadir + '\') not found')
    sys.stdout.flush()
    return False

  files = os.listdir(datadir)
  for f in files:
    if f.endswith('dcm'):
      if os.path.exists(datadir+f):
        print('\t'+f)
        if copyfiles:
          if os.path.isfile(datadir+f):
            shutil.copy(datadir+f,import_image_dir+'/'+f+'.sdcopen')
          else:
            shutil.copytree(datadir+f,import_image_dir+'/'+f+'.sdcopen')
        else:
          shutil.move(datadir+f,import_image_dir+'/'+f+'.sdcopen')
      else:
        print('Warning: not a file \''+f+'\'')
      sys.stdout.flush()
    if f.endswith('sdcopen'):
      if os.path.exists(datadir+f):
        print('\t'+f)
        if copyfiles:
          if os.path.isfile(datadir+f):
            shutil.copy(datadir+f,import_image_dir+'/'+f)
          else:
            shutil.copytree(datadir+f,import_image_dir+'/'+f)
        else:
          shutil.move(datadir+f,import_image_dir+'/'+f)
      else:
        print('Warning: not a file \''+f+'\'')
      sys.stdout.flush()
  return True


# ******************************************
# move files of src directory to dst and remove src
def move(src,dst,deldir=True):
  if not os.path.isdir(src):
    print('Warning: directory src not found\n\t'+src)
    sys.stdout.flush()
    return False
  if not os.path.isdir(dst):
    print('Warning: directory dst not found\n\t'+dst)
    sys.stdout.flush()
    return False

  print('Moving files from\n\t'+src+'\n\tto\n\t'+dst)
  sys.stdout.flush()
  files = os.listdir(src)
  for f in files:
    if os.path.isfile(dst+f):
      print('Warning: file existing; overwriting\n\t'+dst+f)
      sys.stdout.flush()
      os.remove(dst+f)
    shutil.move(src+f,dst)

  if deldir:
    os.rmdir(src)

  return True


# ******************************************
# copy receiver noise scans from
# dname = '/usr/g/recon_calib/ExamData/' + str(exam) +'/'
# fpos=='last'  download last (matching) file (default)
# fpos=='all'   download all (matching) files
# fpos=='first' download first (matching) file
# returns fnames
def copy_noise(exam,datadir='',fpos='last'):
  # directory with noise data
  dname = '/usr/g/recon_calib/ExamData/' + str(exam) +'/'
  if not os.path.isdir(dname):
    print('Error: noise directory not found: ' + dname)
    return None

  # get noise files
  fnames = [os.path.join(dname,fname) for fname in os.listdir(dname)]
  fnames.sort(key=os.path.getmtime)
  if not fnames:
    print('Error: no noise files found in: ' + dname)
    return None

  # select file positions: last (default), all, first
  if fpos=='last':
    del fnames[0:-1]
  else:
    if fpos=='first':
      del fnames[1:]
    else:
      if fpos!='all':
        print('Warning: unkown fpos; downloading all found files')

  # actual coyping
  print('Copying ' + str(len(fnames)) + ' noise files')
  for fname in fnames:
    print('\t' + fname)
    shutil.copyfile(fname,datadir+os.path.basename(fname))

  sys.stdout.flush()
  return fnames


# ******************************************
# starting subprocess + printing output (normal Linux command)
def run_os(cmd,verb=True):
  if verb:
    print('cmd = ')
    print(cmd)
  sys.stdout.flush()
  mascriproc = subprocess.Popen(cmd,shell=False,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
  (stdoutdata, stderrdata) = mascriproc.communicate()
  stdoutdata = stdoutdata.decode('ISO-8859-1')
  stderrdata = stderrdata.decode('ISO-8859-1')
  if verb:
    print('stdoutdata =\n' + str(stdoutdata))
    print('********************************')
    print('stderrdata =\n' + str(stderrdata))
    print('********************************')
    sys.stdout.flush()


# ******************************************
# starting subprocess + printing output
def run(cmd,mascri=''):
  # setting up mmatlab runtime environment
  if mascri=='':
    mascri = setup_mcr([],get_what())
    if mascri=='':
      print('Error: mascri not found')
      if len(sys.argv)>6:
        sys.stdout.close()
      sys.exit(-1)
    if print_debug():
      print_env_vars()

  cmd.insert(0,mascri)
  print('cmd = ')
  print(cmd)
  sys.stdout.flush()
  mascriproc = subprocess.Popen(cmd,shell=False,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
  (stdoutdata, stderrdata) = mascriproc.communicate()
  stdoutdata = stdoutdata.decode('ISO-8859-1')
  stderrdata = stderrdata.decode('ISO-8859-1')
  print('stdoutdata =\n' + str(stdoutdata))
  print('********************************')
  print('stderrdata =\n' + str(stderrdata))
  print('********************************')
  sys.stdout.flush()


# ******************************************
# start matlab reconstruction via recon_vre.py on vre
def run_on_vre(cmd):
  if True:
    # enable xwindows for robustness, in case anything is plotted on vre by accident
    run_os(['xhost','+','vre'])
  cmd.insert(0,dir('fidall')+'bin/recon_vre.py')
  cmd.insert(0,'vre')
  try:
    subprocess.Popen(['rsh'],shell=False,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    cmd.insert(0,'rsh')
  except:
    cmd.insert(0,'ssh')
    cmd.insert(1,'-oStrictHostKeyChecking=no')
  print('cmd = ')
  print(cmd)
  sys.stdout.flush()

  rsh = subprocess.Popen(cmd,shell=False,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
  (stdoutdata, stderrdata) = rsh.communicate()
  stdoutdata = stdoutdata.decode('ISO-8859-1')
  stderrdata = stderrdata.decode('ISO-8859-1')
  print('stdoutdata =\n' + str(stdoutdata))
  print('********************************')
  print('stderrdata =\n' + str(stderrdata))
  print('********************************')
  sys.stdout.flush()


# ******************************************
# run matlab_scripter in background and feed cmd via fifo pipe in there
def run_fifo(cmd_inp):
  import pathlib
  import datetime
  import time
  
  # os.environ['DISPLAY'] = ':0'
  what = get_what()
  topdir = dir(what)
  
  # check and create fifo pipe
  fnfifo = topdir + 'bin/fifo_host'
  if not pathlib.Path(fnfifo).is_fifo():
    print('fnfifo (=\'' + fnfifo + '\') not existing; creating')
    cmd = 'mkfifo ' + fnfifo
    print(cmd)
    sys.stdout.flush()
    os.system(cmd)
    time.sleep(1)
    if not pathlib.Path(fnfifo).is_fifo():
      print('Error: cannot create ' + fnfifo)
      return False
  
  # check if matlab_scripter process running
  spawnmascri = True
  cmd = ['ps','aux']
  print(cmd)
  sys.stdout.flush()
  sp = subprocess.Popen(cmd,shell=False,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
  lines = sp.stdout.readlines()
  if sys.version_info[0]>2:  # python3 requires conversion to str
    lines = [x.decode('ISO-8859-1') for x in lines]
  if lines == []:
    print('ERROR: lines empty')
    return False
  else:
    for line in lines:
      # print(line)
      if re.search('matlab_scripter',line):
        if re.search('fifopipe',line):
          if re.search(what,line):
            print('\'matlab_scripter fifopipe\' process running for ' + what)
            spawnmascri = False
  
  # start matlab_scripter process in background
  now = datetime.datetime.now()
  if spawnmascri:
    print('matlab_scripter fifo_pipe not running: starting in background')
    fnlog = topdir + 'log/' + now.strftime("%Y%m%d_%H%M%S") + '_mascri.txt'
    print('logfile = ' + fnlog)
    fidlog = open(fnlog,'wt')
    if fidlog.closed:
      xmessage('Error opening file ' + fnlog)
    # cmd = 'nohup ' + topdir + 'bin/matlab.py \"fifopipe ' + fnfifo + '\" &'
    cmd = [topdir + 'bin/matlab.py','fifopipe ' + fnfifo]
    print(cmd)
    sys.stdout.flush()
    fidlog.write('Starting matlab_scripter fifopipe ' + fnfifo + '\n')
    subprocess.Popen(cmd,shell=False,stdout=fidlog,stderr=fidlog)
    # os.system(cmd)
  
  # start actual process
  fnstatus = topdir + 'log/status/' + now.strftime("%Y%m%d_%H%M%S") + '_status.txt'
  cmd_inp.append(fnstatus)
  print('cmd_inp= ')
  print(cmd_inp)
  sys.stdout.flush()
  
  # write into fifo pipe to start recon in background
  fidfifo = open(fnfifo,'wt')
  fidfifo.write(str(cmd_inp))
  fidfifo.close()
  
  print('looking for status file ' + fnstatus)
  sys.stdout.flush()
  
  for ll in range(300):
    if os.path.isfile(fnstatus):
      print('status file found; continuing (ll= ' + str(ll) + ' [s])')
      fidstatus = open(fnstatus,'rt')
      if fidstatus.closed:
        xmessage('Error opening file ' + fnstatus)
        return False
      line = fidstatus.readline()
      fidstatus.close()
      print('status = ' + line)
      if int(line)==0:
        xmessage('Error:\n' + '\n'.join(cmd_inp))
        return False
      else:
        if int(line)==1:
          print('Recon successful')
        else:
          xmessage('Error: run_fifo:\nunknown status ' + line)
      break
    else:
      time.sleep(1)
    if (ll==299):
      print('Warning: reconstruction time >5min (continuing)')
  
  sys.stdout.flush()
  
  return True
  

# ******************************************
# get runnum,exam,series info via selected dicom file
def get_dcm_hdr(codes):
  dcm_parser = '/export/home/sdc/bin/dicmCompParser'
  fnselect = '/usr/tmp/sdc_selection'

  # check that program and file exist
  if not os.path.isfile(dcm_parser):
    xmessage('Error: dcm_parser not found (\''+dcm_parser+'\')')
    return(['-1'])
  if not os.path.isfile(fnselect):
    xmessage('Error: fnselect not found (\''+fnselect+'\')')
    return(['-1'])

  # get selected dicom file
  print('Extracting selected dicom file from')
  print('\t'+fnselect)
  fid = open(fnselect,'rt')
  line = fid.readline()
  # print(line)
  fid.close()
  fndcm = re.findall('^\S+',line)[0]
  print('selected file\n\t'+fndcm)

  # extract parameters from fndcm
  cmd = [dcm_parser,'-i',fndcm]
  print(cmd)
  sp = subprocess.Popen(cmd,shell=False,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
  lines = sp.stdout.readlines()
  if sys.version_info[0]>2:  # python3 requires conversion to str
    lines = [x.decode('ISO-8859-1') for x in lines]
  val = []
  if lines == []:
    print('ERROR: lines empty')
    return(['-1'])
  else:
    print(codes)
    for code in codes:
      for line in lines:
        if re.search(code,line):
          # print(line)
          tmp = re.findall('\[[\w\.\-]+\]$',line)[0]
          # print(tmp)
          val.append(tmp[1:-1])
  if len(codes) != len(val):
    print('Warning: len(codes) != len(val)')

  return(val)


# ******************************************
# display xmessage
def xmessage(msg,prnt=True):
  if prnt:
    print(msg)
  sys.stdout.flush()
  cmd = 'echo \"' + msg + '\" | xmessage -file - &'
  os.system(cmd)


# ******************************************
# view images in directory
def view(datadir):
  import glob
  if re.search('\.h5$',datadir) is not None:
    datadir = datadir[0:-3]
  if re.search('\.7$',datadir) is not None:
    datadir = datadir[0:-2]
      
  if doview():
    if glob.glob(datadir + '*.png'):
      if os.path.isfile('/usr/bin/eog'):
        cmd = '/usr/bin/eog'
      else:
        if os.path.isfile('/usr/bin/display'):
          cmd = '/usr/bin/display -resize 50%'
        else:
          cmd = '/usr/bin/gimp'
      cmd = cmd + ' ' + datadir + '*.png &' 
      print(cmd)
      os.system(cmd)
    else:
      xmessage('No png files found in ' + datadir)
  else:
    print('doview=False: skipping image viewing')


# ******************************************
# generate/link reference for HOS tool
def hos_reference(fname):
  # fname=P99999.7 file+dir with B0map info in HOS tool format
  # dir with HOS reference files
  dnref = '/w/config/ho_shim/reference/'

  # check/complete params fname
  if not fname.endswith('.params'):
    fname = fname + '.params'
  print('Reading pars from ' + fname)
  fid = open(fname,'rt')
  if fid.closed:
    print('Error opening file ' + fname)
    xmessage('Error opening file ' + fname)
    return False
  
  # read P99999.7.params file
  lines = fid.readlines()
  if len(lines)<6:
    print('#lines (=' + str(len(lines)) + ') in param file < 6')
    fid.close()
    return False
  mtx = [int(i) for i in lines[0].split(' ')]
  slthick = float(lines[4])
  opfov = float(lines[5])
  fov = [opfov,opfov,slthick*mtx[2]/10]
  fid.close()
  
  # check/generate HOS reference file
  # fovstr = '{:s}'.format('_'.join(['{:.0f}'.format(x) for x in fov]))
  # mtxstr = '{:s}'.format('_'.join(['{:d}'.format(x) for x in mtx]))
  fovstr = '_'.join(['%.0f' % x for x in fov])
  mtxstr = '_'.join(['%d' % x for x in mtx])
  refname = 'HO_MAP_GRD_X.FOV_' + fovstr + '.RES_' + mtxstr
  print('HOS reference file: ' + dnref + refname)
  if os.path.isfile(dnref + refname):
    print('HOS reference file existing')
  else:
    print('HOS reference file not existing -> generating')
    if not os.path.isfile(dnref + 'hos.cal'):
      print('Error: hos.cal missing -> lvshim system calibration required')
      return False
    cmd = ['/usr/g/tools/lvshim/runRefmaps',str(fov[0]),str(fov[1]),str(slthick/10),\
      str(mtx[0]),str(mtx[1]),str(mtx[2])]
    run_os(cmd)
    if not os.path.isfile(dnref + refname):
      print('Error generating HOS reference')
      print('File not generated/found: ' + dnref + refname)
      return False

  # check/set symbolic links
  fnlinks = ['reference.field','reference.mask','reference.params']
  suf = ['','.mask','.params']
  if os.path.isfile(dnref + fnlinks[0]):
    fnlinto = os.readlink(dnref + fnlinks[0])
  else:
    fnlinto = ''
  print(fnlinks[0])
  print(fnlinto)

  if dnref + refname == fnlinto:
    print('Symbolic links already set correctly')
  else:
    print('Setting symbolic links ' + dnref + fnlinks[0] + ' + .mask + .params')
    for ll in range(3):
      if os.path.islink(dnref + fnlinks[ll]):
        os.unlink(dnref + fnlinks[ll])
      os.symlink(refname + suf[ll], dnref + fnlinks[ll])

  return True

