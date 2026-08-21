#! /usr/bin/env python

# -*- coding: utf-8 -*-
"""
Created on Thu Aug  6 17:00:17 2015

@author: Josh Kaggie


To use

import GEwav as ge

ge.list_directory('/usr/g/mrraw/')

or from command line

python GEwav.py

or if chmod u+x GEwav.py

./GEwav.py

"""

import numpy as np
path = '/home/jk636/c/data/'
filen = 'ak_grad20.wav'
import glob



    
try:
    mlab.addpath('/home/jk636/c/code/matlab')
except:
    pass


deffile = 'P10752.7'


import os



#################
import glob
import numpy as np

Stanford0 = np.dtype([('desc', np.uint8, 256), ('N.gpts', '>u2', 1), ('N.groups', '>u2', 1)]) #ugh Stanford, why not a normal order?



fdl_scaling = {
    'te': 1.,
    'tr': 1.,        
    'ph': 10000.,                
    'th': 10000.,                        
    'fl': 1.,
    'du': 1.,                        
    'sc': np.inf
}

fdlheader_dtype = np.dtype([
            ('npars',np.int32,1),
            ('scaling',np.double,1),
            ('spare',np.int32,1),
            #('list',np.int32),
            ])




def read_par(parfile='simulations'):
    import scipy.io as sio
    par = sio.loadmat(parfile)
    grad_delay = par['par']['grad'][0][0][0][0][0][0][0]
    method_path = par['par']['path'][0][0]['method'][0][0][0]
    B1slice_path = par['par']['path'][0][0]['B1slice'][0][0][0]
    dictionary_path = par['par']['path'][0][0]['dictionary'][0][0][0]
    par['par']['recon'][0]
    par0['bp'][0][0][0].dtype
    T1 = par0['bp'][0][0][0]['T1'][0][0]
    T2 = par0['bp'][0][0][0]['T1'][0][0]
    frequency = par0['bp'][0][0][0]['T1'][0][0]
    mrftype = par0['bp'][0][0][0]['T1'][0][0]
    method_name = par0['bp'][0][0][0]['method_name'][0][0]
    

def rescale_array(inarray, vmin=None, vmax=None,center_zero = False,offset = 0):
    if vmin == None:
        vmin = np.min(inarray)
    if vmax == None:
        vmax = np.max(inarray)
    inarray = inarray-np.min(inarray)
    inarray = inarray/np.max(inarray)
    if center_zero:
        inarray = inarray-0.5
        inarray = inarray*2.
    inarray = (vmax-vmin)*inarray+vmin+offset
    return inarray


def read_MRF_pisasim(seqname = 'Gris_979', path = 'simulations/'):    
    full_dictionary = np.fromfile(path+seqname+'_real.dat',dtype=np.single) + 1j*np.fromfile(path+seqname+'_imag.dat',dtype=np.single)
    norm_dictionary = np.fromfile(path+seqname+'_norm.dat',dtype=np.single)
    svd_dictionary = np.fromfile(path+seqname+'_sVd_compression.dat',dtype=np.single)
    return full_dictionary, norm_dictionary, svd_dictionary
    
    

def read_fdl(fname='fisp_mrf_flip.fdl',ptype='tr',x=np.zeros(10)):
    fid = open(fname, 'rb')
    fdlheader = np.fromfile(fid,dtype=fdlheader_dtype,count=1)
    listval = np.fromfile(fid,dtype=np.int32)    
    listval = listval/fdlheader['scaling']
    return fdlheader, listval


def write_fdl(fname='vap_tr90b.fdl',ptype='tr',x=np.zeros(10)):
    fid = open(fname, 'wb')
        
    nparams = len(x)
    scaling = fdl_scaling[ptype]
    spare = 0
    
    fid.write(np.array(nparams).astype(np.int32).tobytes())
    fid.write(np.array(scaling).astype(np.double).tobytes())
    fid.write(np.array(spare).astype(np.int32).tobytes())
    
    fid.write(x.astype(np.int32).tobytes())
    return 0


def write_mod(fname='flip_rf1.mod',x=np.zeros(10)):
    #function write_mod(fname,x)
    #%WRITE_MOD  Write vector x into binary file readable for fidall
    #%  write_mod(fname,x)
    #%   fname  Output filename (-> "/export/home/sdc/fidall/flip_rf1.mod")
    #%       x  binary vector 
    #%          (0..90): flip angles [deg]
    #% 07/2012 Rolf Schulte
    #if (nargin<1), help(mfilename); return; end;
    if type(fname) != type(''):
        print ('Filename not string')
        print fname
        return -1
    if np.any(x<0):
        print 'x<0'
        return -1
    if np.any(x>90):
        print 'x>90'
        return -1
    x = np.array(x).astype(np.float32)
    fid = open(fname,'wb')
    lenx = np.array(len(x))
    fid.write(lenx.tobytes())
    fid.write(x.tobytes())
    fid.close()


def write_modi(fname='vartr<CV20>.modi',x=np.zeros(10)):
    if type(fname) != type(''):
        print ('Filename not string')
        print fname
        return -1        
    x = np.array(x).astype(np.int32)
    fid = open(fname,'wb')
    fid.write(x.tobytes())
    fid.close()


def read_modi(fname='varTR1.modi',dtype=np.int32):
    fid = open(fname, 'rb')
    return np.fromfile(fname,dtype=dtype)
    
def read_mod(fname='varFA1.mod',dtype=np.float32):
    fid = open(fname, 'rb')
    n = np.fromfile(fname,dtype=dtype,count=1)
    thearray = np.fromfile(fname,dtype=dtype)    
    return thearray #, n


def check_wav(wav):
    print np.max(np.max(np.diff(wav,axis=0),axis=0),axis=0)
    print np.max(np.max(np.cumsum(wav,axis=0),axis=0),axis=0)
    print np.max(np.max(wav,axis=0),axis=0)
    print np.sum(np.sum(wav**2,axis=0),axis=0)
    
def read_mat_wav(filename=None,fixN=None):
    import scipy.io as sio
    data=sio.loadmat(filename)
    trajx  = data['k'].real
    trajy  = data['k'].imag
    grad = data['k']
    ffs = data['out']
    kpts = ffs['N'][0][0]['kpts']
    wave = ffs['wave'][0][0]
    #ff['wave'][0][0].shape
    #ff['N']
    if fixN:
        wave = wave[:fixN]
    return wave



def read_ak_wav(fname):
    '''
    COMPLETE
    read_ak_wav(fname)
    
    #reads in a filename with the Stanford waveform type
    
    Original matlab code written by Rolf Schulte.  Python code implemented by Josh Kaggie Aug 2015.
     
    %READ_AK_WAV  Reads waveforms stored to external file using Stanford format
    % [grad,bw,fov,desc,N,params] = read_ak_wav(fname)
    % 
    %      fname  File name of output file                       [string]
    %       grad  Gradient waveforms with 4us time resolution    [T/m]
    %             [#pts/interleave,#interleaves,#groups]
    %             with: #groups = 2 for 2d-imaging, =3 for 3d-imaging
    %         bw  Full bandwidth (opuser0 = 2d3*oprbw)           [Hz]
    %        fov  Field-of-view (old=nucleus; new=1H equivalent) [m]
    %       desc  Description string (256 chars)
    %          N  N.gpts   # input gradient pts/interleave
    %             N.gpts   # input gradient pts/interleave
    %             N.kpts   # readout pts
    %             N.groups # groups
    %             N.intl   # interleaves
    %             N.params # parameters
    %     params  Header file parameters (scanner units)
    %             grad_type fov N.intl gmax N.gpts gdt N.kpts kdt 0 0 0
    %
    % 06/2013     Rolf Schulte
    %  Josh Kaggie 2015-2017, python version.  '''
    
    Stanford0 = np.dtype([('desc', np.uint8, 256), ('N.gpts', '>u2', 1), ('N.groups', '>u2', 1)]) #ugh Stanford, why not a normal order?    
    #Stanford0 = np.dtype([('desc', np.uint8, 256), ('N.gpts', np.uint16, 1), ('N.groups', np.uint16, 1)]) #ugh Stanford, why not a normal order?    
    f = open(fname, "rb")
    f.seek(0)
    data1 = np.fromfile(f, dtype=Stanford0,count=1)
    Ngpts = data1['N.gpts']
    Ngroups = data1['N.groups']
    desc = data1['desc'].tostring()
    Nintl = np.fromfile(f, dtype = '>u2', count = Ngroups)#N.intl = fread(fid, N.groups, 'uint16');
    #Nintl = np.fromfile(f, dtype = np.uint16, count = Ngroups)#N.intl = fread(fid, N.groups, 'uint16');
    intl = Nintl[0]
    #print 'josh'
    print Ngpts, intl, Ngroups
    #print str(desc)
    print
    #print 'what'
    Nparams = np.fromfile(f, dtype = '>u2', count = 1)#N.params = fread(fid, 1, 'uint16');
    #Nparams = np.fromfile(f, dtype = np.uint16, count = 1)#N.params = fread(fid, 1, 'uint16');
    params = np.fromfile(f, dtype = '>f8', count = Nparams) #float64
    #params = np.fromfile(f, dtype = np.float64, count = Nparams) #float64
    #wav = np.fromfile(f,'>i2') #int8
    wav = np.fromfile(f,'>i2') #int8
    print Ngroups, intl, Ngpts, np.shape(wav)
    try:
        grad = np.transpose(np.reshape(wav,[int(Ngroups), int(intl), int(Ngpts)  ]));
    except:
        print('Error occurred in transposing wavfile.')
        print(primes(Ngpts), primes(intl), primes(Ngroups))
        return wav
    if np.any(np.abs(wav)) > 32769:
        print 'Waveform not scaled to 32767'    
    if np.any(grad[-1]) != 1:
        print 'Stop bit not set to 1'
    grad[-1] = 0
    grad = grad/32767.*params[3]/100.
    bw = 1e6/params[7]
    fov = params[1]/100.
    #[grad,bw,fov,desc,N,params,wave]
    N = {'gpts': Ngpts, 'groups':Ngroups, 'intl':Nintl}
    #grad = np.transpose(grad)
    print np.shape(wav)
    print 'The shape of the gradient is ', np.shape(grad)
    print 'The maximum of the gradient is ', np.max(grad)
    return {'grad':grad, 'bw':bw, 'fov':fov, 'desc':desc,'N':N, 'params':params, 'wav':wav}
   

def write_ak_wav(fname, grad,bw=62500.,fov=0.02,desc='',n_kpts = [],params=None, rfs_fastfid = -10, ak_dynwf_mode = 0):
    '''
    COMPLETE
    write_ak_wav(fname, grad,bw,fov,desc,n_kpts, rfs_fastfid, ak_dynwf_mode)  
    %WRITE_AK_WAV  Saves waveforms to external file using Stanford format
    % out = write_ak_wav(fname,grad,bw,fov,desc,n_kpts,rfs_fastfid,ak_dynwf_mode)
    % INPUT (SI Units)
    %        fname  File name of output file                      [string]
    %         grad  Gradient waveforms with 4us time resolution   [T/m]
    %               [#pts/interleave,#interleaves,#groups]
    %               with: #groups = 2 for 2d-imaging, =3 for 3d-imaging
    %               if #groups=1 and complex grad -> 2D
    %           bw  Full sampling bandwidth; scalar               [Hz]
    %          fov  Field-of-view relative to 1H; scalar          [m]
    %          des  Description string (opt; up to 254 chars)
    %       n_kpts  Number of acq points (opuser1) (opt)
    %               default = ceil(N.gpts*gdt/kdt)
    %  rfs_fastfid  Reduce TR of fidall (opt; default=-10): 
    %               -10=fidall default,0=orig,1+2=reduced delay,3=no crusher
    %ak_dynwf_mode  Mode for dynamic waveforms (scanmode=2) (opt) 
    %               0=read from file in rsp (default); 1=read into AGP memory
    %
    % OUTPUT (GE Scanner Units):
    % out               Output structure
    % out.descr         Description
    % out.N.gpts        # input gradient pts/interleave.  Number of gradient points in a single gradient, per interleave.
    % out.N.kpts        # readout pts
    % out.N.groups      # groups (=2; real and imaginary)
    % out.N.intl        # interleaves
    % out.N.params      # parameters (=14)
    % out.parms=params  Header file parameters
    % out.wave=wave     Output waveform
    %
    % Difference to previous versions: save indices to worst (3D) trajectories 
    % with >1 interleaves into header for heating calculations (pgen on host)
    %
    % (C)   2007        Adam@mrsrl.stanford.edu 
    % (M)   07/2014     Rolf Schulte
    Python implementation, Aug 2015, Josh Kaggie
    write_ak_wav('/home/jk636/c/data/testwav.wav', hi['grad'], hi['bw'], hi['fov'], hi['desc'])
    '''        
    #FIXED PARAM
    gdt = 4. # gradiest sampling time in us  
    ###### Josh:  UGgggggghhhh....   It will be nice to rename all of this.
    grad_type = 0  
    if len(np.shape(grad)) == 2:
        if np.iscomplex(grad):
            print 'Data is complex!'
            tmp = np.zeros(np.hstack([np.shape(grad), 2]))
            tmp[:,:,0] = np.real(grad)
            tmp[:,:,1] = np.imag(grad)
            grad = tmp
        else:
            print np.shape(grad), ' is the shape of the gradient!'
    elif len(np.shape(grad)) != 3:
        print 'warning: gradient shape error'
    if np.any(grad[0])>1e-10:
        print 'warning: first element not zero'
    if np.any(grad[-1])>1e-10:
        print 'warning: last element not zero'
    ## derived params
    grad = np.transpose(grad)
    grad = grad*1.0e2 # T/m to G/cm    
    gmax = np.max(np.abs(grad))
    smax = np.max(np.abs(np.diff(grad)))/gdt*1.0e3
    fov = fov*1.0e2 # m to cm
    kdt = 1.0e6 / bw # sampling time in us
    if np.abs(np.ceil(kdt)-kdt)> 1.0e-10:
        print 'kdt is not a multiple of 1us'
    Ngroups, Nintl, Ngpts = np.shape(grad)
    if len(np.shape(n_kpts)):
        n_kpts = np.ceil(Ngpts*gdt/kdt)        
    if Ngpts > 32766:
        print 'N.gpts exceeds maximum of fidall'        
    if Nintl > 16382:
        print 'N.intl exceeds maximum of fidall'
    print 'gmax = ' + str((gmax*10)) + ' mT/m (or 10 G/cm)'
    print 'smax = ' + str((smax*10)) + ' mT/m/s (or 10 G/cm/ms)'        
    desc = desc.ljust(253)[:253] + '\n\f\n'   #pad to 253 characters, then trim to 253 characters, then add 3 characters        
    #determine worst (sum-of-squares) sub-waveform for pulsegen on host    
    '''
    wf4pg = np.zeros((1,3))    
    if Nintl>1:
        tmp = np.sum(grad**2,axis=0)
        for l1 in np.arange(Ngroups):
            sos_val = 0
            for l2 in np.arange(Nintl):
                if tmp[0,l2,l1] > sos_val*1.001:
                    sos_val = tmp[0,l2,l1]
                    wf4pg[0,l1] = l2-1'''
    #if Nintl.any()>1
    print 'something something worst trajectories'

    print gmax
    wav = 2.*( (2.**14-1.)/gmax*grad)
    wav[:,:,-1] = 1.
    wav = wav.astype(np.int16)
    #return wav,grad
    
    print np.shape(wav), wav.max(),wav.min(), ' wavvals'
    #Stanford0 = np.dtype([('desc', np.uint8, 256), ('N.gpts', '>u2', 1), ('N.groups', '>u2', 1)]) #ugh Stanford, why not a normal order?    
    f = open(fname, "wb")
    f.seek(0)
    #desc_int = np.array(map(ord,str(desc))).astype(np.uint8)
    f.write(desc)
    #np.dtype([('desc', np.uint8, 256), ('N.gpts', '>u2', 1), ('N.groups', '>u2', 1)])
    print Ngpts, Ngroups, Nintl, Nparams, params
    f.write(np.array(Ngpts).astype('>u2').tostring())
    f.write(np.array(Ngroups).astype('>u2').tostring())
    #f.write(desc)
    Nintlgroups = np.outer(Nintl,np.ones(Ngroups))
    f.write(np.array(Nintlgroups).astype('>u2').tostring())
    f.write(np.array(Nparams).astype('>u2').tostring())
    f.write(np.array(params).astype('>f8').tostring())
    print 'wavshape is this ', wav.shape
    f.write(np.array(wav).flatten().astype('>i2').tostring())
    return wav




def read_rho(filename = '/home/jk636/c/Downloads/bir4-12-3-20.rho', hasHeader = True,unwrap = True):
    ''' reads in pulse shape from GE waveform'''
    f = open(filename, 'rb')
    f.seek(0)
    h = None
    h1 = None
    if hasHeader:
        h1 = np.fromfile(f, dtype=np.uint8, count = 4 )
        h = np.fromfile(f, dtype = np.int32, count = 3 )
        header = np.fromfile(f, dtype = np.uint8, count = 64).tostring()
        f.seek(64,0) # skip to 64 bytes after beginning of file
        pulse = np.fromfile(f,dtype = '>i2')
    else:
        pulse = np.fromfile(f,dtype = '>i2')
    #DataByteCount = h(1) - 64 - 4 - 4;    
    
    if unwrap:
        pulse = np.unwrap(pulse/pulse.max()*np.pi)
    return pulse, h, h1, header
    
 
    

class read_MR_rawdata():
    def __init__(self,filename= '/home/jk636/fidalltest/P10752.7', save_baseline=False, phaselist = False, echolist = False, slicelist = False, rcvrlist = False):
        ''' hello?'''
        #from mlab.releases import latest_release as matlab
        ec=0;
        non_standard_pfile=0;
        #% Determine if the input parameter IN_NAME points to a directory of a file
        f = open(filename)
        
        all_data = readMRheader(filename, header_only = True)
        headers = all_data.datas['rdb_hdr_short']
        self.da_xres = headers['da_xres']
        self.da_yres = headers['da_yres']
        self.nslices = headers['nslices']
        self.nechoes = headers['nechoes']
        self.point_size = headers['point_size']
        self.coils = headers['dab(1)']-headers['dab(0)']+1
        
        self.elements = self.da_xres*2
        self.frame_sz = 2* self.da_xres*self.point_size
        self.echo_sz = self.frame_sz*self.da_yres
        self.slice_sz = self.echo_sz*self.nechoes
        self.mslice_sz = self.slice_sz*self.coils
        self.data_elements = self.da_xres*self.da_yres
        
        self.all_data = all_data
        self.raw_data = np.fromstring(all_data.datas['remainder'],np.int32)
        self.baseline_offset = self.da_xres*2*self.point_size
        self.data_elements = self.da_xres*(self.da_yres-1)*2
        #ff=a[0x7001,0x1008]
        
def read_MR_dicomrawdata(filename='/home/jk636/tmp/P70144.dcm'):
    import dicom
    filename = deffile
    data = dicom.read_file(filename)
    rdb_raw_data = data[0x7001,0x1005]
    #(0019, 10ab)  #979, userdata4
    #0019,10a8 $1024, userdata1
    spokes= data[0x0019,0x10a8]
    points= data[0x0019,0x10ab]
    raw=np.fromstring(rdb_raw_data.value,np.int32)
    


class readMRheader():
    def __init__(self,filename= '/home/jk636/fidalltest/P10752.7',headerID = 'all', 
                 datatype = 'raw',header_only = True):
        header = {}    
        self._DV24_sizes = np.array([157276,4096,16384,16384,98304,2052,2052,2048,1500,1960,2560,2448,7488])
        _DV24_sizes = self._DV24_sizes
        #1, rdb_hdr, 2 per_pass, 3 unlock_raw, 4 data_acq_tab, 5 nex_tab, 6 nex_abort_tab, 7 tool, 
        #8 prescan, 9 exam, 10 series, 11 image, 12 grad_data
        actual_offsets = np.cumsum(self._DV24_sizes)    
        #[ formatID, endianID, header_list, header_lengths, rdbm_rev ] = get_file_format( filename,datatype );
        letsread = open(filename,'rb')
        datas = {}
        letsread.seek(0)
        datas['rdb_hdr_short']=np.fromfile(letsread,_DV24_rdbhdr,1)
        letsread.seek(self._DV24_sizes[0],1)
        letsread.seek(0)
        self.run_char = datas['rdb_hdr_short']['run_char'].tostring().replace('\x00', ' ' )
        self.scan_date = datas['rdb_hdr_short']['scan_date'].tostring().replace('\x00', ' ' )
        self.scan_time = datas['rdb_hdr_short']['scan_time'].tostring().replace('\x00', ' ' )
        itemX = 1
        rdb_hdr = letsread.read(_DV24_sizes[itemX])
        self.rdb_header = datas['rdb_hdr_short']

        itemX = itemX + 1
        datas['per_pass'] = letsread.read(_DV24_sizes[itemX])
        self.per_pass = datas['per_pass']
        itemX = itemX + 1
        datas['unlock_raw'] = letsread.read(_DV24_sizes[itemX])
        datas['data_acq_tab'] = letsread.read(_DV24_sizes[itemX+1])
        datas['data_acq_tab_short']=np.fromstring(datas['data_acq_tab'][:_DV24_data_acq_tab.itemsize],_DV24_data_acq_tab)
        self.data_acq_tab = datas['data_acq_tab_short']
        #_DV24_data_acq_tab
        datas['nex_tab'] = letsread.read(_DV24_sizes[itemX+2])
        datas['nex_abort_tab'] = letsread.read(_DV24_sizes[itemX+3])
        datas['tool'] = letsread.read(_DV24_sizes[itemX+4])
        datas['prescan'] = letsread.read(_DV24_sizes[itemX+5])
        datas['exam'] = letsread.read(_DV24_sizes[itemX+6])
        print len(datas['exam']),_DV24_exam_header.itemsize
        datas['scan_date'] = datas['rdb_hdr_short']['scan_date']
        
        #
        '''
        off_data 	[157276]
        off_per_pass 	[4096]
        off_unlock_raw 	[20480]
        off_data_acq_tab 	[36864]
        off_nex_tab 	[135168]
        off_nex_abort_tab 	[137220]
        off_tool 	[139272]
        off_exam 	[142820]
        off_series 	[144780]
        off_image 	[147340]
        off_ps 	[141320]'''
        datas['exam_short']=np.fromstring(datas['exam'][:_DV24_exam_header.itemsize],_DV24_exam_header)
        self.exam = datas['exam_short']
        self.ex_desc = datas['exam_short']['ex_desc'].tostring().replace('\x00','_')
        self.scan_time = self.rdb_header['scan_time'].tostring().replace('\x00','_')
        self.scan_date = self.rdb_header['scan_date'].tostring().replace('\x00','_')
        self.coils = self.rdb_header['dab(1)']-self.rdb_header['dab(0)']+1
        

        #_DV24_exam_header 
        datas['series'] = letsread.read(_DV24_sizes[9])
        datas['series_short']=np.fromstring(datas['series'][:_DV24_series_header.itemsize],_DV24_series_header)
        self.series = datas['series_short']
        self.scan_name = self.series['c_key'].tostring().replace('\x00','_')
        datas['image'] = letsread.read(_DV24_sizes[10])
        #_DV24_image_header, GEwav._DV24_image_header.itemsize, gg['image'][:788]
        datas['image_short']=np.fromstring(datas['image'][:_DV24_image_header.itemsize],_DV24_image_header)
        self.image = datas['image_short']
        if header_only:
            return
        datas['grad_data'] = letsread.read(_DV24_sizes[11])
        #datas['grad_data_short']=np.fromstring(datas['grad_data'][:_DV24_gradient_header.itemsize],_DV24_gradient_header)
        #self.gradient_data = datas['grad_data_short']
        datas['remainder'] = letsread.read()
        datas['basics'] = {}
        datas['basics']['n_coils'] = datas['rdb_hdr_short']['dab(1)']-datas['rdb_hdr_short']['dab(0)']+1
        datas['basics']['x_n'] = datas['rdb_hdr_short']['da_xres']
        datas['basics']['y_n'] = datas['rdb_hdr_short']['da_yres']-1
        datas['basics']['n_slices'] = datas['rdb_hdr_short']['nslices']
        datas['basics']['tot_n'] = tot_n = datas['basics']['n_coils']*datas['basics']['n_coils']*datas['basics']['n_coils']*2
        self.basics = datas['basics']
        self.n_x = datas['basics']['x_n']
        self.n_y = datas['basics']['y_n']
        self.n_coils = datas['basics']['n_coils']
        self.n_slices = datas['basics']['n_slices']
        self.n_tot = self.n_x*self.n_y*self.n_coils*self.n_slices*2*4
        #print 
        #letsread.seek(-1*int(self.n_tot))
        #self.tester=np.fromfile(letsread,np.int32)
        import os
        self.filesize = os.fstat(letsread.fileno()).st_size
        self.letsread = letsread
        letsread.seek(int(self.filesize-self.n_tot),0)
        self.act_data = np.fromfile(letsread,dtype=np.int32)
        self.act_data = self.act_data[::2]+1j*self.act_data[1::2]
        #int(len(self.act_data)/self.n_y/self.n_x/self.n_slices/self.n_coils),
        ###THIS SHOULD BE ITT??!!!
        #self.act_data = np.reshape(self.act_data, [ int(self.n_y),int(self.n_coils), int(self.n_x), int(self.n_slices)])
        self.act_data = np.reshape(self.act_data, [ int(self.n_x),int(self.n_coils), int(self.n_y), int(self.n_slices)])
        if False:
            try:
                #datas['remainder'] = np.fromstring(remainder,'i4')
                datas['remainder2'] = np.fromstring(datas['remainder'],np.int32)
            except:
                pass
            try:
                datas['basics']['orig_n'] = len(datas['remainder'])
                datas['basics']['new_n'] = np.floor(datas['basics']['orig_n']/datas['basics']['tot_n'])        
            except:
                pass
        self.datas = datas
        letsread.close()
    def return_MRF_params(self):
        self.n_inters = self.rdb_header['user4']
        self.n_ = self.rdb_header['user4']
    def print_usercvs(self):
        for xi in xrange(48):
            print 'user'+str(xi), '\t',
            print self.rdb_header['user'+str(xi)]
    def print_rdb_names(self):
        for name in self.rdb_header.dtype.names:
            print name, '\t', self.rdb_header[name]
    def print_iamge_names(self):
        for name in self.image.dtype.names:
            print name, '\t', self.image[name]
    def print_exam_names(self):
        for name in self.exam.dtype.names:
            print name, '\t', self.exam[name]
    def print_series_names(self):
        for name in self.series.dtype.names:
            print name, '\t', self.series[name]
    def print_all_headers(self):
        self.print_exam_names()
        self.print_rdb_names()
        self.print_series_names()
    
    def print_scanparams(self):
        #print self.rdb_header['user'+str(xi)]        
        print self.rdb_header['da_yres']
        print self.rdb_header['da_xres']
        print self.rdb_header['nslices']
    def enum(self,inputval):
        enum_dtype(inputval)
    
    
def rdb_hdr_to_stuff(rdb_hdr):
    Stanford0 = np.dtype([
                ('desc', np.uint8, 256),    
                ('desc1', np.uint8, 256),    
                ])
    
 
def enum_dtype(inputval):
    import numpy as np
    for name in inputval.dtype.names:
        print name,
        print inputval[name]

    
 
def get_file_format(filename='/home/jk636/c/data/P30208.7', datatype='raw'):
    formatID = ''
    endianID = ''
    header_list = []
    header_lengths = []
    if datatype=='image':
        f = open(filename,'rb')
        test = np.fromfile(f,dtype='>u8',count=3352).tostring()
        test_110 = test.find('DICM')
        testo = test
        test = test.find('IMGF')-1
        f.close()
        if test == 0 :
            formatID = 'ximg'
            endianID = 'ieee-be'
            header_list = [12, 11, 8, 9, 10]
            header_lengths = get_header_lengths(1)
        elif test ==3228:
            formatID = 'listsel'
            endianID = 'ieee-be'
            header_list = [11,8,9,10,12]
            header_lengths = get_header_lengths(1)
        elif test_110<3552:
            formatID = 'DICM'
            endianID = ' '
            header_list = ' ' 
            header_lengths = ' ' 
        else:
            formatID = 'ximg'
            endianID = 'ieee-le'
            header_list =  [12,11,8,9,10]
            header_lengths = get_header_lengths(1)
        rdbm_rev = 'image'
        return testo
            
            
    print formatID, endianID, header_list, header_lengths


    
    

def make_recon_directory(template_recon_path='/home/jk636/d/MRF/mrfv2/mrfv2/reconserver50',
                         template_out_path='/home/jk636/d/MRF/mrfv2/mrfv2/tmp50',recon_out_path='/home/jk636/d/MRF/mrfv2/mrfv2/reconserver54',
                         tmp_out_path='/home/jk636/d/MRF/mrfv2/mrfv2/tmp54'):
    import os

    try:
        os.stat(recon_out_path)
    except:
        os.mkdir(recon_out_path)       
    
    try:
        os.stat(tmp_out_path)
    except:
        os.mkdir(tmp_out_path)       
    import glob
    recon_files = glob.glob(template_recon_path)
    for filename in recon_files:
        os.symlink(filename, filename.replace(template_recon_path,recon_out_path))

    tmp_files = glob.glob(template_tmp_path)
    for filename in recon_files:
        os.symlink(filename, filename.replace(template_recon_path,recon_out_path))    


def fix_path(path=''):
    if path[-1] != '/':
        path = path+'/'
    return path

#header = GEwav.readMRheader('/home/jk636/c/data/JamieConesKneesMRF/P47104.7',header_only=True)
def list_directory(path='/usr/g/mrraw',pattern='P*7',printlist=True,returnlist=False, new_method = True):
    import glob
    path = fix_path(path)    
    files = glob.glob(path+pattern)
    #files.sort(key=os.path.getmtime)
    returnvalues = {}
    for filename in files:
        header = readMRheader(filename,header_only=True) 
        if new_method:
            printlist = False
            print filename[-8:]
            print header.ex_desc
            print header.scan_name, header.scan_date, header.scan_time
            print header.rdb_header['da_yres'],
            print header.rdb_header['da_xres'],
            print header.rdb_header['nslices'],
            print header.rdb_header['broad_band_select'],
            print header.rdb_header['ps_aps_or_mps'],
            print header.rdb_header['ps_mps_freq'],
            print header.rdb_header['ps_aps_freq']
            print '#################################\n\n'
        if printlist:
            print filename[-8:]
            print header.rdb_header['scan_time'].tostring().replace('\x00',''),
            print header.rdb_header['scan_date'].tostring().replace('\x00',''),
            print header.exam['ex_desc'].tostring().replace('\x00','')[:20],
            #print header.series['se_desc'].tostring().replace('\x00','')[-20:],
            print header.series['c_key'].tostring().replace('\x00',''),
            print header.rdb_header['da_yres'],
            print header.rdb_header['da_xres'],
            print header.rdb_header['nslices'],
            print header.rdb_header['broad_band_select'],
            print header.rdb_header['ps_aps_or_mps'],
            print header.rdb_header['ps_mps_freq'],
            print header.rdb_header['ps_aps_freq']
        if returnlist:
            returnvalues[filename] = {}
            returnvalues[filename]['filename'] = filename[-8:]
            returnvalues[filename]['scan_time'] = header.rdb_header['scan_time'].tostring().replace('\x00',''),
            returnvalues[filename]['scan_date'] =  header.rdb_header['scan_date'].tostring().replace('\x00',''),
            returnvalues[filename]['ex_desc'] =  header.exam['ex_desc'].tostring().replace('\x00',''),
            #print header.series['se_desc'].tostring().replace('\x00','')[-20:],
            returnvalues[filename]['c_key'] =  header.series['c_key'].tostring().replace('\x00',''),
            returnvalues[filename]['da_xres'] =  header.rdb_header['da_xres'],
            returnvalues[filename]['da_yres'] =  header.rdb_header['da_yres'],
            returnvalues[filename]['user1'] =  header.rdb_header['user1'],
            returnvalues[filename]['user4'] =  header.rdb_header['user4'],
            returnvalues[filename]['user7'] =  header.rdb_header['user7'],
            returnvalues[filename]['user2'] =  header.rdb_header['user2'],
            returnvalues[filename]['user3'] =  header.rdb_header['user3'],
            returnvalues[filename]['user5'] =  header.rdb_header['user5'],
            returnvalues[filename]['user6'] =  header.rdb_header['user6'],
            returnvalues[filename]['user8'] =  header.rdb_header['user8'],
            returnvalues[filename]['user10'] =  header.rdb_header['user10'],
            returnvalues[filename]['user9'] =  header.rdb_header['user9'],

            returnvalues[filename]['nslices'] =  header.rdb_header['nslices'],
            returnvalues[filename]['broad_band_select'] =  header.rdb_header['broad_band_select'],
            returnvalues[filename]['ps_aps_or_mps'] =  header.rdb_header['ps_aps_or_mps'],
            returnvalues[filename]['ps_mps_freq'] =  header.rdb_header['ps_mps_freq'],
            returnvalues[filename]['ps_aps_freq'] =  header.rdb_header['ps_aps_freq']
    if returnlist:
        return returnvalues
    
def find_all_pfiles(start_path = '/home/jk636/d',outputfile = 'Poutput.txt'):
    outputfile = open(outputfile,'w')
    list_all  = {}
    for path, subs, files in os.walk(start_path):
        list_all[path] = []
        print_path = True    
        for filename in files:
            if 'P' in filename and '.7' in filename:
                if print_path:
                    outputfile.write(path+'\n\n')
                    print_path = False
                print filename
                try:
                    header = readMRheader(path+'/'+filename,header_only=True) 
                    description =  header.ex_desc.strip().replace('_','')#re.sub(r'[^\x00-\x7f]',r'', header.ex_desc.strip().replace('_','')) 
                    outputfile.write( filename + ' ' + header.scan_name +' '+ description[-20:] +' '+ header.scan_date+' '+ header.scan_time + ' ' + str(header.rdb_header['ps_mps_freq']) + '\n\n')
                except:
                    pass
                list_all[path].append(filename)
                    
    

def list_sort(path='/home/jk636/d/Projects/MRF Repeat Brain/data/1.5TCam'):
    lists = list_directory(path,returnlist=True)
    
    for key in np.sort(lists.keys()):
        val = lists[key]
        print key[-8:],val['ex_desc'][0][-10:], val['ps_mps_freq'][0][0]


def recursive_search(directory, pattern):
    import os, fnmatch
    for root, dirs, files in os.walk(directory):
        for basename in files:
            if fnmatch.fnmatch(basename, pattern):
                filename = os.path.join(root, basename)
                yield filename


def  get_header_lengths( format_num ):
    #% Define the lengths for the individual headers present in each file format;
    #% zero lengths mean the header is not in the file.
    #%       XIMG   LIST
    
    a  = [    [ 0   ,   0],#   %  1 rdb_rec
        [0 ,     0 ],#  %  2 per_pass
        [0  ,    0 ],#  %  3 unlock_raw
        [0   ,   0 ],#  %  4 data_acq_tab
        [0    ,  0 ],#  %  5 nex_tab
        [0     , 0 ],#  %  6 nex_abort_tab
        [0     , 0 ],#  %  7 tool
        [1040  , 1040],#   %  8 exam
        [1028  , 1028],#   %  9 series
        [1044  , 1044],#   % 10 image
        [116   , 116 ],#  % 11 suite
        [124   , 124 ]]#  ] ;   % 12 pixel
    a = np.array(a)
    header_lengths = a[ int(format_num-1),:]
    return header_lengths 
            

def copy_header(file_header='/media/jk636/e/P35328.7', file_raw = '/media/jk636/e/P23552.7', file_out = '/media/jk636/e/P35329.7')   :
    f1 = open(file_header, 'rb')
    f2 = open(file_raw, 'rb')
    f3 = open(file_out, 'wb')
    
    import GEwav
    GEwav._DV24_sizes[0]
    
    f1.seek(0)
    f3.write(f1.read(GEwav._DV24_sizes[0]))
    f2.seek(GEwav._DV24_sizes[0])
    f3.write(f2.read())
    f3.close()
    f2.close()
    f1.close()


_DV24_coil = np.dtype([
    ('cttEntry(id).logicalCoilName', '>i1',128),
    ('cttEntry(id).clinicalCoilName', '>i1',32),
    ('cttEntry(id).configUID', np.uint32),
    ('cttEntry(id).coilConnector', np.int32),
    ('cttEntry(id).isActiveConfig', np.uint32),
    ('cttEntry(id).channelTranslationMap', np.int16,32)
    ])

## converted from read_rdb_hdr
_DV24_rdbhdr = np.dtype( [
  ('rdbm_rev',   np.float32),
  ('run_int',   '<i4'),
  ('scan_seq',   '<i2'),
  ('run_char', '>i1'  , 6), 
  ('scan_date', '<i1'  , 10), 
  ('scan_time', '<i1'  , 8), 
  ('logo', '<i1'  , 10), 
  ('file_contents',   '<i2'),
  ('lock_mode',   '<i2'),
  ('dacq_ctrl',   '<i2'),
  ('recon_ctrl',   '<i2'),
  ('exec_ctrl',   '<u2'),
  ('scan_type',   '<i2'),
  ('data_collect_type',   '<i2'),
  ('data_format',   '<i2'),
  ('recon',   '<i2'),
  ('datacq',   '<i2'),
  ('npasses',   '<i2'),
  ('npomp',   '<i2'),
  ('nslices',   '<u2'),
  ('nechoes',   '<i2'),
  ('navs',   '<i2'),
  ('nframes',   '<i2'),
  ('baseline_views',   '<i2'),
  ('hnover',   '<i2'),
  ('frame_size',   '<u2'),
  ('point_size',   '<i2'),
  ('vquant',   '<i2'),
  ('cheart',   '<i2'),
  ('ctr',   np.float32),
  ('ctrr',   np.float32),
  ('initpass',   '<i2'),
  ('incrpass',   '<i2'),
  ('method_ctrl',   '<i2'),
  ('da_xres',   '<u2'),
  ('da_yres',   '<i2'),
  ('rc_xres',   '<i2'),
  ('rc_yres',   '<i2'),
  ('im_size',   '<i2'),
  ('rc_zres',   '<i4'),
  ('raw_pass_size_deprecated',   '<i4'),
  ('sspsave_deprecated',   '<i4'),
  ('udasave_deprecated',   '<i4'),
  ('fermi_radius',   np.float32),
  ('fermi_width',   np.float32),
  ('fermi_ecc',   np.float32),
  ('clip_min',   np.float32),
  ('clip_max',   np.float32),
  ('default_offset',   np.float32),
  ('xoff',   np.float32),
  ('yoff',   np.float32),
  ('nwin',   np.float32),
  ('ntran',   np.float32),
  ('scalei',   np.float32),
  ('scaleq',   np.float32),
  ('rotation',   '<i2'),
  ('transpose',   '<i2'),
  ('kissoff_views',   '<i2'),
  ('slblank',   '<i2'),
  ('gradcoil',   '<i2'),
  ('ddaover',   '<i2'),
  ('sarr',   '<i2'),
  ('fd_tr',   '<i2'),
  ('fd_te',   '<i2'),
  ('fd_ctrl',   '<i2'),
  ('algor_num',   '<i2'),
  ('fd_df_dec',   '<i2'),
  ('dab(0)',   '<i2'),
  ('dab(1)',   '<i2'),  
  ('dab(2)',   '<i2'),
  ('dab(3)',   '<i2'),
  ('dab(4)',   '<i2'),
  ('dab(5)',   '<i2'),
  ('dab(6)',   '<i2'),
  ('dab(7)',   '<i2'),  
  ('user0',   np.float32),
  ('user1',   np.float32),
  ('user2',   np.float32),
  ('user3',   np.float32),
  ('user4',   np.float32),
  ('user5',   np.float32),
  ('user6',   np.float32),
  ('user7',   np.float32),
  ('user8',   np.float32),
  ('user9',   np.float32),
  ('user10',   np.float32),
  ('user11',   np.float32),
  ('user12',   np.float32),
  ('user13',   np.float32),
  ('user14',   np.float32),
  ('user15',   np.float32),
  ('user16',   np.float32),
  ('user17',   np.float32),
  ('user18',   np.float32),
  ('user19',   np.float32),
  ('v_type' , np.int32),
  ('v_coefxa' , np.float32),
  ('v_coefxb' , np.float32),
  ('v_coefxc' , np.float32),
  ('v_coefxd' , np.float32),
  ('v_coefya' , np.float32),
  ('v_coefyb' , np.float32),
  ('v_coefyc' , np.float32),
  ('v_coefyd' , np.float32),
  ('v_coefza' , np.float32),
  ('v_coefzb' , np.float32),
  ('v_coefzc' , np.float32),
  ('v_coefzd' , np.float32),
  ('vm_coef1' , np.float32),
  ('vm_coef2' , np.float32),
  ('vm_coef3' , np.float32),
  ('vm_coef4' , np.float32),
  ('v_venc' , np.float32),
  ('spectral_width' , np.float32),
  ('csi_dims' , np.int16),
  ('xcsi' , np.int16),
  ('ycsi' , np.int16),
  ('zcsi' , np.int16),
  ('roilenx' , np.float32),
  ('roileny' , np.float32),
  ('roilenz' , np.float32),
  ('roilocx' , np.float32),
  ('roilocy' , np.float32),
  ('roilocz' , np.float32),
  ('numdwell' , np.float32),
  ('ps_command' , np.int32),
  ('ps_mps_r1' , np.int32),
  ('ps_mps_r2' , np.int32),
  ('ps_mps_tg' , np.int32),
  ('ps_mps_freq' , np.uint32),
  ('ps_aps_r1' , np.int32),
  ('ps_aps_r2' , np.int32),
  ('ps_aps_tg' , np.int32),
  ('ps_aps_freq' , np.uint32),
  ('ps_scalei' , np.float32),
  ('ps_scaleq' , np.float32),
  ('ps_snr_warning' , np.int32),
  ('ps_aps_or_mps' , np.int32),
  ('ps_mps_bitmap' , np.int32),
  ('ps_powerspec', '<i1'  , 256), 
  ('ps_filler1' , np.int32),
  ('ps_filler2' , np.int32),
  ('obsolete1' , np.float32,16),  ## BEGIN AUTOCONVERT
  ('obsolete2' , np.float32,16),
  ('halfecho', np.int16),
  ('im_size_y', np.int16),
  ('data_collect_type1', np.int32),
  ('freq_scale', np.float32),
  ('phase_scale', np.float32),
  ('ovl', np.int16),
  ('pclin', np.int16),
  ('pclinnpts', np.int16),
  ('pclinorder', np.int16),
  ('pclinavg', np.int16),
  ('pccon', np.int16),
  ('pcconnpts', np.int16),
  ('pcconorder', np.int16),
  ('pcextcorr', np.int16),
  ('pcgraph', np.int16),
  ('pcileave', np.int16),
  ('hdbestky', np.int16),
  ('pcctrl', np.int16),
  ('pcthrespts', np.int16),
  ('pcdiscbeg', np.int16),
  ('pcdiscmid', np.int16),
  ('pcdiscend', np.int16),
  ('pcthrespct', np.int16),
  ('pcspacial', np.int16),
  ('pctemporal', np.int16),
  ('pcspare', np.int16),
  ('ileaves', np.int16),
  ('kydir', np.int16),
  ('alt', np.int16),
  ('reps', np.int16),
  ('ref', np.int16),
  ('pcconnorm', np.float32),
  ('pcconfitwt', np.float32),
  ('pclinnorm', np.float32),
  ('pclinfitwt', np.float32),
  ('pcbestky', np.float32),
  ('vrgf', np.int32),
  ('vrgfxres', np.int32),
  ('bp_corr', np.int32),
  ('recv_freq_s', np.float32),
  ('recv_freq_e', np.float32),
  ('hniter', np.int32),
  ('fast_rec', np.int32),
  ('refframes', np.int32),
  ('refframep', np.int32),
  ('scnframe', np.int32),
  ('pasframe', np.int32),
  ('user_usage_tag', np.uint32),
  ('user_fill_mapMSW', np.uint32),
  ('user_fill_mapLSW', np.uint32),
  ('user20', np.float32),
  ('user21', np.float32),
  ('user22', np.float32),
  ('user23', np.float32),
  ('user24', np.float32),
  ('user25', np.float32),
  ('user26', np.float32),
  ('user27', np.float32),
  ('user28', np.float32),
  ('user29', np.float32),
  ('user30', np.float32),
  ('user31', np.float32),
  ('user32', np.float32),
  ('user33', np.float32),
  ('user34', np.float32),
  ('user35', np.float32),
  ('user36', np.float32),
  ('user37', np.float32),
  ('user38', np.float32),
  ('user39', np.float32),
  ('user40', np.float32),
  ('user41', np.float32),
  ('user42', np.float32),
  ('user43', np.float32),
  ('user44', np.float32),
  ('user45', np.float32),
  ('user46', np.float32),
  ('user47', np.float32),
  ('user48', np.float32),
  ('pcfitorig', np.int16),
  ('pcshotfirst', np.int16),
  ('pcshotlast', np.int16),
  ('pcmultegrp', np.int16),
  ('pclinfix', np.int16),
  ('pcconfix', np.int16),
  ('pclinslope', np.float32),
  ('pcconslope', np.float32),
  ('pccoil', np.int16),
  ('vvsmode', np.int16),
  ('vvsaimgs', np.int16),
  ('vvstr', np.int16),
  ('vvsgender', np.int16),
  ('zip_factor', np.int16),
  ('maxcoef1a', np.float32),
  ('maxcoef1b', np.float32),
  ('maxcoef1c', np.float32),
  ('maxcoef1d', np.float32),
  ('maxcoef2a', np.float32),
  ('maxcoef2b', np.float32),
  ('maxcoef2c', np.float32),
  ('maxcoef2d', np.float32),
  ('maxcoef3a', np.float32),
  ('maxcoef3b', np.float32),
  ('maxcoef3c', np.float32),
  ('maxcoef3d', np.float32),
  ('ut_ctrl', np.int32),
  ('dp_type', np.int16),
  ('arw', np.int16),
  ('vps', np.int16),
  ('mcReconEnable', np.int16),
  ('fov', np.float32),
  ('te', np.int32),
  ('te2', np.int32),
  ('dfmrbw', np.float32),
  ('dfmctrl', np.int32),
  ('raw_nex', np.int32),
  ('navs_per_pass', np.int32),
  ('dfmxres', np.int32),
  ('dfmptsize', np.int32),
  ('navs_per_view', np.int32),
  ('dfmdebug', np.int32),
  ('dfmthreshold', np.float32),
  ('grid_control', np.int16),
  ('b0map', np.int16),
  ('grid_tediff', np.int16),
  ('grid_motion_comp', np.int16),
  ('grid_radius_a', np.float32),
  ('grid_radius_b', np.float32),
  ('grid_max_gradient', np.float32),
  ('grid_max_slew', np.float32),
  ('grid_scan_fov', np.float32),
  ('grid_a2d_time', np.float32),
  ('grid_density_factor', np.float32),
  ('grid_display_fov', np.float32),
  ('fatwater', np.int16),
  ('fiestamlf', np.int16),
  ('app', np.int16),
  ('rhncoilsel', np.int16),
  ('rhncoillimit', np.int16),
  ('app_option', np.int16),
  ('grad_mode', np.int16),
  ('pfile_passes', np.int16),
  ('asset', np.int32),
  ('asset_calthresh', np.int32),
  ('asset_R', np.float32),
  ('coilConfigUID', np.uint32),
  ('asset_phases', np.int32),
  ('scancent', np.float32),
  ('position', np.int32),
  ('entry', np.int32),
  ('lmhor', np.float32),
  ('last_slice_num', np.int32),
  ('asset_slice_R', np.float32),
  ('asset_slabwrap', np.float32),
  ('dwnav_coeff', np.float32),
  ('dwnav_cor', np.int16),
  ('dwnav_view', np.int16),
  ('dwnav_corecho', np.int16),
  ('dwnav_sview', np.int16),
  ('dwnav_eview', np.int16),
  ('dwnav_sshot', np.int16),
  ('dwnav_eshot', np.int16),
  ('a3dwin_type', np.int16),
  ('a3dwin_apod', np.float32),
  ('a3dwin_q', np.float32),
  ('ime_scic_enable', np.int16),
  ('clariview_type', np.int16),
  ('ime_scic_edge', np.float32),
  ('ime_scic_smooth', np.float32),
  ('ime_scic_focus', np.float32),
  ('clariview_edge', np.float32),
  ('clariview_smooth', np.float32),
  ('clariview_focus', np.float32),
  ('scic_reduction', np.float32),
  ('scic_gauss', np.float32),
  ('scic_threshold', np.float32),
  ('ectricks_no_regions', np.int32),
  ('ectricks_input_regions', np.int32),
  ('psc_reuse', np.int16),
  ('left_blank', np.int16),
  ('right_blank', np.int16),
  ('acquire_type', np.int16),
  ('retro_control', np.int16),
  ('etl', np.int16),
  ('pcref_start', np.int16),
  ('pcref_stop', np.int16),
  ('ref_skip', np.int16),
  ('extra_frames_top', np.int16),
  ('extra_frames_bot', np.int16),
  ('multiphase_type', np.int16),
  ('nphases', np.int16),
  ('pure', np.int16),
  ('pure_scale', np.float32),
  ('off_data', np.int32),
  ('off_per_pass', np.int32),
  ('off_unlock_raw', np.int32),
  ('off_data_acq_tab', np.int32),
  ('off_nex_tab', np.int32),
  ('off_nex_abort_tab', np.int32),
  ('off_tool', np.int32),
  ('off_exam', np.int32),
  ('off_series', np.int32),
  ('off_image', np.int32),
  ('off_ps', np.int32),
  ('off_spare_b', np.int32),
  ('new_wnd_level_flag', np.int32),
  ('wnd_image_hist_area', np.int32),
  ('wnd_high_hist', np.float32),
  ('wnd_lower_hist', np.float32),
  ('pure_filter', np.int16),
  ('cfg_pure_filter', np.int16),
  ('cfg_pure_fit_order', np.int16),
  ('cfg_pure_kernelsize_z', np.int16),
  ('cfg_pure_kernelsize_xy', np.int16),
  ('cfg_pure_weight_radius', np.int16),
  ('cfg_pure_intensity_scale', np.int16),
  ('cfg_pure_noise_threshold', np.int16),
  ('wienera', np.float32),
  ('wienerb', np.float32),
  ('wienert2', np.float32),
  ('wieneresp', np.float32),
  ('wiener', np.int16),
  ('flipfilter', np.int16),
  ('dbgrecon', np.int16),
  ('ech2skip', np.int16),
  ('tricks_type', np.int32),
  ('lcfiesta_phase', np.float32),
  ('lcfiesta', np.int16),
  ('herawflt', np.int16),
  ('herawflt_befnwin', np.int16),
  ('herawflt_befntran', np.int16),
  ('herawflt_befamp', np.float32),
  ('herawflt_hpfamp', np.float32),
  ('heover', np.int16),
  ('pure_correction_threshold', np.int16),
  ('swiftenable', np.int32),
  ('numslabs', np.int16),
  ('numCoilConfigs', np.uint16),
  ('ps_autoshim_status', np.int32),
  ('dynaplan_numphases', np.int32),
  ('medal_cfg', np.int16),
  ('medal_nstack', np.int16),
  ('medal_echo_order', np.int16),
  ('medal_kernel_up', np.int16),
  ('medal_kernel_down', np.int16),
  ('medal_kernel_smooth', np.int16),
  ('medal_start', np.int16),
  ('medal_end', np.int16),
  ('rcideal', np.uint32),
  ('rcdixproc', np.uint32),
  ('df', np.float32),
  ('bw', np.float32),
  ('te1_deprecated', np.float32),
  ('esp_deprecated', np.float32),
  ('feextra', np.int32),
  ('raw_pass_size', np.uint64),
  ('sspsave', np.uint64),
  ('udasave', np.uint64),
  ('vibrant', np.int16),
  ('asset_torso', np.int16),
  ('asset_alt_cal', np.int32),
  ('kacq_uid', np.int32),
  ('coils',_DV24_coil,4), ### MAY BE WRONG!
  ('psc_ta', np.int32),
  ('disk_acq_ctrl', np.int32),
  ('asset_localTx', np.int32),
  ('a3dscale', np.float32),
  ('broad_band_select', np.int32),
  ('scanner_mode', np.int16),
  ('numbvals', np.int16),
  ('numdifdirs', np.int16),
  ('difnext2', np.int16),
  ('difnextab(id)', np.int16,100),
  ('channel_combine_method', np.int16),
  ('nexForUnacquiredEncodes', np.int16),
  ('a2dscale', np.float32),
  ('dd_mode', np.int16),
  ('dd_q_ta_offset', np.int16),
  ('dd_q_phase_delay', np.float32),
  ('dacq_ctrl_chksum', np.uint32),
  ('patient_checksum', np.uint32),
  ('rcidealctrl', np.uint32),
  ('echotimes(id)', np.float32,64),
  ('asl_perf_weighted_scale', np.int16),
  ('echo_pc_extra_frames_bot', np.uint16),
  ('echo_pc_ctrl', np.uint32),
  ('echo_pc_ylines', np.uint16),
  ('echo_pc_primary_yfirst', np.uint16),
  ('echo_pc_reverse_yfirst', np.uint16),
  ('echo_pc_zlines', np.uint16),
  ('echo_pc_yxfitorder', np.uint16),
  ('channel_combine_filter_type', np.int16),
  ('mavric_control', np.uint32),
  ('mavric_ImageType', np.uint32),
  ('mavric_bin_separation', np.int32),
  ('mavric_b0_offset(id)', np.float32,40),
  ('channel_combine_filter_width', np.float32),
  ('channel_combine_filter_beta', np.float32),
  ('low_pass_nex_filter_width', np.float32),
  ('aps_tg_status', np.uint32),
  ('cal_pass_set_vector', np.int32),
  ('cal_nex_vector', np.int32),
  ('cal_weight_vector', np.int32),
  ('pure_filtering_mode', np.int32),
  ('pure_lambda', np.float32),
  ('pure_tuning_factor_surface', np.float32),
  ('pure_tuning_factor_body', np.float32),
  ('pure_derived_cal_fraction', np.float32),
  ('pure_derived_cal_reapodization', np.float32),
  ('noncart_grid_factor', np.float32),
  ('noncart_dual_traj', np.int32),
  ('noncart_traj_kmax_ratio', np.int32),
  ('noncart_traj_merge_start', np.int32),
  ('noncart_traj_merge_end', np.int32),
  ('oversamplingfactor', np.float32),
  ('nspokes_highres', np.int32),
  ('nspokes_lowres', np.int32),
  ('nrefslices', np.int32),
  ('psmde_pixel_offset', np.int32),
  ('off_grad_data', np.int32),
  ('hoecc', np.int32),
  ('hoec_fit_order', np.int32),
  ('esp', np.int32),
  ('excess', np.int16,322),  
  ])

_DV24_image_header = np.dtype([
  ('autoSubParam_seriesUidToSubtract', '<i1',32),
  ('autoSubParam_imageNoToSubtract', '<i4'),
  ('autoSubParam_destSeriesNo', '<i4'),
  ('autoSubParam.destImageNo', '<i4'),
  ('autoSubParam.dummy', '<i4'),
  ('double_padding(id)', np.float64,32),
  ('dfov', np.float32),
  ('dfov_rect', np.float32),
  ('sctime', np.float32),
  ('slthick', np.float32),
  ('scanspacing', np.float32),
  ('loc', np.float32),
  ('tbldlta', np.float32),
  ('nex', np.float32),
  ('reptime', np.float32),
  ('saravg', np.float32),
  ('sarpeak', np.float32),
  ('pausetime', np.float32),
  ('vbw', np.float32),
  ('user0', np.float32),
  ('user1', np.float32),
  ('user2', np.float32),
  ('user3', np.float32),
  ('user4', np.float32),
  ('user5', np.float32),
  ('user6', np.float32),
  ('user7', np.float32),
  ('user8', np.float32),
  ('user9', np.float32),
  ('user10', np.float32),
  ('user11', np.float32),
  ('user12', np.float32),
  ('user13', np.float32),
  ('user14', np.float32),
  ('user15', np.float32),
  ('user16', np.float32),
  ('user17', np.float32),
  ('user18', np.float32),
  ('user19', np.float32),
  ('user20', np.float32),
  ('user21', np.float32),
  ('user22', np.float32),
  ('proj_ang', np.float32),
  ('concat_sat', np.float32),
  ('user23', np.float32),
  ('user24', np.float32),
  ('x_axis_rot', np.float32),
  ('y_axis_rot', np.float32),
  ('z_axis_rot', np.float32),
  ('ihtagfa', np.float32),
  ('ihtagor', np.float32),
  ('ihbspti', np.float32),
  ('rtia_timer', np.float32),
  ('fps', np.float32),
  ('vencscale', np.float32),
  ('dbdt', np.float32),
  ('dbdtper', np.float32),
  ('estdbdtper', np.float32),
  ('estdbdtts', np.float32),
  ('saravghead', np.float32),
  ('neg_scanspacing', np.float32),
  ('user25', np.float32),
  ('user26', np.float32),
  ('user27', np.float32),
  ('user28', np.float32),
  ('user29', np.float32),
  ('user30', np.float32),
  ('user31', np.float32),
  ('user32', np.float32),
  ('user33', np.float32),
  ('user34', np.float32),
  ('user35', np.float32),
  ('user36', np.float32),
  ('user37', np.float32),
  ('user38', np.float32),
  ('user39', np.float32),
  ('user40', np.float32),
  ('user41', np.float32),
  ('user42', np.float32),
  ('user43', np.float32),
  ('user44', np.float32),
  ('user45', np.float32),
  ('user46', np.float32),
  ('user47', np.float32),
  ('user48', np.float32),
  ('RegressorVal', np.float32),
  ('SliceAsset', np.float32),
  ('PhaseAsset', np.float32),
  ('sarValues(id)', np.float32, 4),
  ('shim_fov(id)', np.float32,2),
  ('shim_ctr_R(id)', np.float32,2),
  ('shim_ctr_A(id)', np.float32,2),
  ('shim_ctr_S(id)', np.float32,2),
  ('dim_X', np.float32),
  ('dim_Y', np.float32),
  ('pixsize_X', np.float32),
  ('pixsize_Y', np.float32),
  ('ctr_R', np.float32),
  ('ctr_A', np.float32),
  ('ctr_S', np.float32),
  ('norm_R', np.float32),
  ('norm_A', np.float32),
  ('norm_S', np.float32),
  ('tlhc_R', np.float32),
  ('tlhc_A', np.float32),
  ('tlhc_S', np.float32),
  ('trhc_R', np.float32),
  ('trhc_A', np.float32),
  ('trhc_S', np.float32),
  ('brhc_R', np.float32),
  ('brhc_A', np.float32),
  ('brhc_S', np.float32),
  ('menc', np.float32),
  ('normal_L', np.float32),
  ('normal_P', np.float32),
  ('normal_S', np.float32),
  ('osf', np.float32),
  ('fermi_radius', np.float32),
  ('fermi_width', np.float32),
  ('fermi_ecc', np.float32),
  ('float_padding', np.float32,25),
  ('cal_fldstr', np.uint32),
  ('user_usage_tag', np.uint32),
  ('user_fill_mapMSW', np.uint32),
  ('user_fill_mapLSW', np.uint32),
  ('im_archived', np.int32),
  ('im_complete', np.int32),
  ('int_padding1(id)', np.int32,34),
  ('im_datetime', np.int32),
  ('im_actual_dt', np.int32),
  ('tr', np.int32),
  ('ti', np.int32),
  ('te', np.int32),
  ('te2', np.int32),
  ('tdel', np.int32),
  ('mindat', np.int32),
  ('obplane', np.int32),
  ('slocfov', np.int32),
  ('obsolete1', np.int32),
  ('obsolete2', np.int32),
  ('user_bitmap', np.int32),
  ('iopt', np.int32),
  ('psd_datetime', np.int32),
  ('rawrunnum', np.int32),
  ('intr_del', np.int32),
  ('im_lastmod', np.int32),
  ('im_pds_a', np.int32),
  ('im_pds_c', np.int32),
  ('im_pds_u', np.int32),
  ('thresh_min1', np.int32),
  ('thresh_max1', np.int32),
  ('thresh_min2', np.int32),
  ('thresh_max2', np.int32),
  ('numslabs', np.int32),
  ('locsperslab', np.int32),
  ('overlaps', np.int32),
  ('slop_int_4', np.int32),
  ('dfax', np.int32),
  ('fphase', np.int32),
  ('offsetfreq', np.int32),
  ('b_value', np.int32),
  ('iopt2', np.int32),
  ('ihtagging', np.int32),
  ('ihtagspc', np.int32),
  ('ihfcineim', np.int32),
  ('ihfcinent', np.int32),
  ('num_seg', np.int32),
  ('oprtarr', np.int32),
  ('averages', np.int32),
  ('station_index', np.int32),
  ('station_total', np.int32),
  ('iopt3', np.int32),
  ('delAcq', np.int32),
  ('rxmbloblen', np.int32),
  ('rxmblob_pad', np.int32),
  ('im_no', np.int32),
  ('imgrx', np.int32),
  ('temp_phases', np.int32),
  ('driver_freq', np.int32),
  ('driver_amp', np.int32),
  ('driverCyc_Trig', np.int32),
  ('MEG_dir', np.int32),
  ('rescan_time', np.int32),
  ('spokesPerSeg', np.int32),
  ('recoveryTime', np.int32),
  ('t2PrepTE', np.int32),
  ('hoecc', np.int32),
  ('user_bitmap2', np.int32),
  ('int padding2', np.int32,20),
  ('imatrix_X', np.int16),
  ('imatrix_Y', np.int16),
  ('im_exno', np.uint16),
  ('img_window', np.uint16),
  ('img_level', np.int16),
  ('numecho', np.int16),
  ('echonum', np.int16),
  ('im_uniq', np.int16),
  ('im_seno', np.int16),
  ('contmode', np.int16),
  ('serrx', np.int16),
  ('screenformat', np.int16),
  ('plane', np.int16),
  ('im_compress', np.int16),
  ('im_scouttype', np.int16),
  ('contig', np.int16),
  ('hrtrate', np.int16),
  ('trgwindow', np.int16),
  ('imgpcyc', np.int16),
  ('obsolete3', np.int16),
  ('obsolete4', np.int16),
  ('obsolete5', np.int16),
  ('mr_flip', np.int16),
  ('cphase', np.int16),
  ('swappf', np.int16),
  ('pauseint', np.int16),
  ('obsolete6', np.int16),
  ('obsolete7', np.int16),
  ('obsolete8', np.int16),
  ('not_used_1', np.int16),
  ('imode', np.int16),
  ('pseq', np.int16),
  ('pseqmode', np.int16),
  ('ctyp', np.int16),
  ('surfctyp', np.int16),
  ('surfcext', np.int16),
  ('supp_tech', np.int16),
  ('slquant', np.int16),
  ('gpre', np.int16),
  ('satbits', np.int16),
  ('scic', np.int16),
  ('satxloc1', np.int16),
  ('satxloc2', np.int16),
  ('satyloc1', np.int16),
  ('satyloc2', np.int16),
  ('satzloc1', np.int16),
  ('satzloc2', np.int16),
  ('satxthick', np.int16),
  ('satythick', np.int16),
  ('satzthick', np.int16),
  ('flax', np.int16),
  ('venc', np.int16),
  ('thk_disclmr', np.int16),
  ('obsolete9', np.int16),
  ('obsolete10', np.int16),
  ('image_type', np.int16),
  ('vas_collapse', np.int16),
  ('proj_alg', np.int16),
  ('echo_trn_len', np.int16),
  ('frac_echo', np.int16),
  ('prep_pulse', np.int16),
  ('cphasenum', np.int16),
  ('var_echo', np.int16),
  ('scanactno', np.int16),
  ('vasflags', np.int16),
  ('integrity', np.int16),
  ('freq_dir', np.int16),
  ('vas_mode', np.int16),
  ('pscopts', np.int16),
  ('obsolete11', np.int16),
  ('obsolete12', np.int16),
  ('obsolete13', np.int16),
  ('unoriginal', np.int16),
  ('interleaves', np.int16),
  ('effechospace', np.int16),
  ('viewsperseg', np.int16),
  ('rbpm', np.int16),
  ('rtpoint', np.int16),
  ('rcvrtype', np.int16),
  ('sarMode', np.int16),
  ('dBdtMode', np.int16),
  ('govBody', np.int16),
  ('sarDefinition', np.int16),
  ('no_shimvol', np.int16),
  ('shim_vol_type', np.int16),
  ('current_phase', np.int16),
  ('art_level', np.int16),
  ('slice_group_number', np.int16),
  ('number_of_slice_groups', np.int16),
  ('show_in_autoview', np.int16),
  ('slice_number_inGroup', np.int16),
  ('specnuc', np.int16),
  ('label_duration', np.uint16),
  ('ihbsoffsetfreq', np.int16),
  ('scale_factor', np.int16),
  ('volume_prop', np.int16),
  ('excitation_mode', np.int16),
    ])


_DV24_exam_header = np.dtype([
  ('firstaxtime', np.float64),
  ('double_padding(id)', np.float64,31),
  ('zerocell', np.float32),
  ('cellspace', np.float32),
  ('srctodet', np.float32),
  ('srctoiso', np.float32),
  ('float_padding(id)', np.float32,32),
  ('ex_delta_cnt', np.int32),
  ('ex_complete', np.int32),
  ('ex_seriesct', np.int32),
  ('ex_numarch', np.int32),
  ('ex_numseries', np.int32),
  ('ex_numunser', np.int32),
  ('ex_toarchcnt', np.int32),
  ('ex_prospcnt', np.int32),
  ('ex_modelnum', np.int32),
  ('ex_modelcnt', np.int32),
  ('patCheckSum', np.int32),
  ('int_padding1(id)', np.int32,31),
  ('numcells', np.int32),
  ('magstrength', np.int32),
  ('patweight', np.int32),
  ('ex_datetime', np.int32),
  ('ex_lastmod', np.int32),
  ('patChecksumType', np.int32),
  ('int_padding2(id)', np.int32,26),
  ('ex_no', np.uint16),
  ('ex_uniq', np.int16),
  ('detect', np.int16),
  ('tubetyp', np.int16),
  ('dastyp', np.int16),
  ('num_dcnk', np.int16),
  ('dcn_len', np.int16),
  ('dcn_density', np.int16),
  ('dcn_stepsize', np.int16),
  ('dcn_shiftcnt', np.int16),
  ('patage', np.int16),
  ('patian', np.int16),
  ('patsex', np.int16),
  ('ex_format', np.int16),
  ('trauma', np.int16),
  ('protocolflag', np.int16),
  ('study_status', np.int16),
  ('short_padding(id)', np.int16,11),
  ('hist', np.int32,  61),
  ('reqnum', np.int32,  13),
  ('refphy', np.int32,33),
  ('diagrad', np.int32,  33),
  ('op', np.int32, 4),
  ('ex_desc', np.int32,  65),
  ('ex_typ', np.int32, 3),
  ('ex_sysid', np.int32, 9),
  ('ex_alloc_key', np.int32,  13),
  ('ex_diskid', np.int8),
  ('hospname', np.int32,  33),
  ('patid', np.int32, 13),
  ('patname', np.int32,10),
  
  ])       
       
       
_DV24_gradient_header = np.dtype([
  ('diffusion_grad_amp(id1,id2)', np.float32,[3,512]),
  ('hoec_bases.hoec_coef', np.float32),
  ('hoec_bases.hoec_xorder', np.int32),
  ('hoec_bases.hoec_yorder', np.int32),
  ('hoec_bases.hoec_zorder', np.int32)
])
      
_DV24_series_header = np.dtype([
  ('double_padding(id)', np.float64,32),
  ('se_pds_a', np.float32),
  ('se_pds_c', np.float32),
  ('se_pds_u', np.float32),
  ('lmhor', np.float32),
  ('start_loc', np.float32),
  ('end_loc', np.float32),
  ('echo1_alpha', np.float32),
  ('echo1_beta', np.float32),
  ('echo2_alpha', np.float32),
  ('echo2_beta', np.float32),
  ('echo3_alpha', np.float32),
  ('echo3_beta', np.float32),
  ('echo4_alpha', np.float32),
  ('echo4_beta', np.float32),
  ('echo5_alpha', np.float32),
  ('echo5_beta', np.float32),
  ('echo6_alpha', np.float32),
  ('echo6_beta', np.float32),
  ('echo7_alpha', np.float32),
  ('echo7_beta', np.float32),
  ('echo8_alpha', np.float32),
  ('echo8_beta', np.float32),
  ('landmark', np.float32),
  ('tablePosition', np.float32),
  ('pure_lambda', np.float32),
  ('pure_tuning_factor_surface', np.float32),
  ('pure_tuning_factor_body', np.float32),
  ('pure_derived_cal_fraction', np.float32),
  ('pure_derived_cal_reapodization', np.float32),
  ('float_padding(id)', np.float32,25),
  ('se_complete', np.int32),
  ('se_numarch', np.int32),
  ('se_imagect', np.int32),
  ('se_numimages', np.int32),
  ('se_delta_cnt', np.int32),
  ('se_numunimg', np.int32),
  ('se_toarchcnt', np.int32),
  ('int_padding1(id)', np.int32,33),
  ('se_datetime', np.int32),
  ('se_actual_dt', np.int32),
  ('position', np.int32),
  ('entry', np.int32),
  ('se_lndmrkcnt', np.int32),
  ('se_lastmod', np.int32),
  ('ExpType', np.int32),
  ('TrRest', np.int32),
  ('TrActive', np.int32),
  ('DumAcq', np.int32),
  ('ExptTimePts', np.int32),
  ('cal_pass_set_vector', np.int32),
  ('cal_nex_vector', np.int32),
  ('cal_weight_vector', np.int32),
  ('pure_filtering_mode', np.int32),
  ('int_padding2(id)', np.int32,29),
  ('se_exno', np.uint16),
  ('echo1_window', np.uint16),
  ('echo2_window', np.uint16),
  ('echo3_window', np.uint16),
  ('echo4_window', np.uint16),
  ('echo5_window', np.uint16),
  ('echo6_window', np.uint16),
  ('echo7_window', np.uint16),
  ('echo8_window', np.uint16),
  ('echo8_level', np.int16),
  ('echo7_level', np.int16),
  ('echo6_level', np.int16),
  ('echo5_level', np.int16),
  ('echo4_level', np.int16),
  ('echo3_level', np.int16),
  ('echo2_level', np.int16),
  ('echo1_level', np.int16),
  ('se_no', np.int16),
  ('se_typ', np.int16),
  ('se_source', np.int16),
  ('se_plane', np.int16),
  ('scan_type', np.int16),
  ('se_uniq', np.int16),
  ('se_contrast', np.int16),
  ('se_pseq', np.int16),
  ('se_sortorder', np.int16),
  ('se_nacq', np.int16),
  ('xbasest', np.int16),
  ('xbaseend', np.int16),
  ('xenhst', np.int16),
  ('xenhend', np.int16),
  ('table_entry', np.int16),
  ('SwingAngle', np.int16),
  ('LateralOffset', np.int16),
  ('GradientCoil', np.int16),
  ('se_subtype', np.int16),
  ('BWRT', np.int16),
  ('assetcal_serno', np.int16),
  ('assetcal_scnno', np.int16),
  ('content_qualifn', np.int16),
  ('purecal_serno', np.int16),
  ('purecal_scnno', np.int16),
  ('ideal', np.int16),
  ('verify_corners', np.int16),
  ('asset_cal_type', np.int16),
  ('pure_compatible', np.int16),
  ('purecal_type', np.int16),
  ('short_padding(id)', np.int16,29),
  ('se_verscre(id)', np.int16,2),
  ('se_verscur(id)', np.int16,2),
  ('se_suid(id)', np.int16,4),
  ('se_alloc_key', np.int32, 1),
  ('c_key', np.int32, 13),
  ('se_diskid', np.int8),
  ('se_desc', np.int32, 65),
  ('pr_sysid', np.int32,  9),
  ('pansysid', np.int32,  9),
  ('anref', np.int32, 3 ),
  ('prtcl', np.int32, 25 ),
  ('start_ras', '<u1'),
  ('end_ras', '<u1'),
])      
   

_DV24_data_acq_tab = np.dtype([
  ('pass_number', np.int16),
  ('slice_in_pass', np.int16),
  ('gw_point1(id)', np.float32,3),
  ('gw_point2(id)', np.float32,3),
  ('gw_point3(id)', np.float32,3),
  ('transpose', np.int16),
  ('rotate', np.int16),
  ('coilConfigUID', np.uint32),
])
   
   
   
def fix_header(header=None):
    if header == None:
      header = '''
      a.ps_scalei = fread(my_file, 1, 'float32');
      a.ps_scaleq = fread(my_file, 1, 'float32');
      a.ps_snr_warning = fread(my_file, 1, 'int32');
      a.ps_aps_or_mps = fread(my_file, 1, 'int32');
      a.ps_mps_bitmap = fread(my_file, 1, 'int32');
      a.ps_powerspec = freadc(my_file, 256);
      a.ps_filler1 = fread(my_file, 1, 'int32');
      a.ps_filler2 = fread(my_file, 1, 'int32');
      for id = 1 : 16
        a.obsolete1(id) = fread(my_file, 1, 'float32');
      end
      '''   
    out = header.replace('a.','(\'')
    out = out.replace(' = fread(','\', np.')
    out = out.replace('\');','),')
    out = out.replace(' = freadc(my_file,','\', np.int32,REPLACEEND')
    def move_n(line=''):
        loc_n = line.find('np.my_file, ')
        loc_n_end = line.find(', \'',loc_n)
        n = line[loc_n+12:loc_n_end]
        #line = line[:loc_n]+line[loc_n_end:]
        #line.replace(', , \'',', np.')
        print line, 'REPLACEEND' in line,
        if 'REPLACEEND' in line:
            line = line.replace('REPLACEEND','').replace(';',',')
        else:
            line = line.replace('my_file, '+ str(n) + ', \'','')            
        if n != '1':
            line = line.replace('),',n+'),')
        print line
        return n,line
    splits = out.split('\n')
    endstr = ''
    for line in splits:
        endstr = endstr+move_n(line)[1]+'\n'
    
    #('int_padding1(id)', np.int32,31),
    return endstr
   
   
_DV24_sizes = np.array([157276,4096,16384,16384,98304,2052,2052,2048,1500,1960,2560,2448,7488])


header_sizes = '''
#7
39984
2048
4096
4096
20480
2052
2052
2048
0
1040
1028
1044
#end

#8
60464
2048
4096
4096
40960
2052
2052
2048
0
1040
1028
1044
#end

#9
61464
2048
4096
4096
40960
2052
2052
2048
0
1040
1536
1536
#end

#10
65560
2048
4096
4096
45056
2052
2052
2048
0
1040
1536
1536
#end

#11
66072
2048
4096
4096
45056
2052
2052
2048
0
1040
2048
1536
#end

#14
135704
2048
16384
16384
90112
2052
2052
2048
0
1040
2048
1536
#end

#14.1
135704
2048
16384
16384
90112
2052
2052
2048
0
1040
2048
1536
#end

#14.2
142356
2048
16384
16384
98304
2052
2052
2048
0
1040
2048
2048
#end

#14.300000
145908
2048
16384
16384
98304
2052
2052
2048
1500
1040
2048
2048
#end

#15.000
145908
2048
16384
16384
98304
2052
2052
2048
1500
1040
2048
2048
#end

#15.001
145908
2048
16384
16384
98304
2052
2052
2048
1500
1040
2048
2048
#end

#16.000
145908
2048
16384
16384
98304
2052
2052
2048
1500
1040
2048
2048
#end

#20.001
145908
2048
16384
16384
98304
2052
2052
2048
1500
1040
2048
2048
#end

#20.002
145932
2048    
16384
16384
98304
2052
2052
2048
1500
1040
2048
2048
#end


#20.003
145932
2072
16384
16384
98304
2052
2052
2048
1500
1040
2048
2048
#end


#20.004
146564
2704
16384
16384
98304
2052
2052
2048
1500
1040
2048
2048
#end


#20.005
149788
4096
16384
16384
98304
2052
2052
2048
1500
1960
2560
2448
#end


#20.006
149788
4096
16384
16384
98304
2052
2052
2048
1500
1960
2560
2448
#end


#20.007
149788
4096
16384
16384
98304
2052
2052
2048
1500
1960
2560
2448
#end


#24.000
157276
4096
16384
16384
98304
2052
2052
2048
1500
1960
2560
2448
7488
#end
'''


            
    







if __name__ == "__main__":
    list_directory()

