import os
import time
import numpy as np
from .getSpiralParams import *
from .getRadialParams import *
from .readPhilipsExports import *

class PhilipsData():
    
    def __init__(self, fname):        
        self.fname = fname
        # check that the path actually exists
        if not os.path.exists(fname):
            return
        # get base filename and extension
        self.base = os.path.splitext(fname)[0]
        self.ext = os.path.splitext(fname)[1]
        
        # define flag for type of data to be read
        self.readType = 0  # .data/.list is default
        if str(self.ext).lower() in ['.list', '.data', '.txt']:
            self.readType = 0
        elif str(self.ext).lower() in ['.par', '.rec', '.xml']:
            self.readType = 1
        elif str(self.ext).lower() in ['.lab', '.raw', '.sin']:
            self.readType = 2
        elif str(self.ext).lower() in ['.cpx']:
            self.readType = 3
        else:
            return
        
        # set up SW and encoding detection
        self._scanner = Scanner('')
        if self.readType == 2:
            self.isMira = False
            sin = readSin(self.base+'.sin')
            self.isMira = sin['isMira']
            lab = readLab(self.base+'.lab', self.isMira)
            std_test = lab['label_type'] == b'LABEL_TYPE_STANDARD'
            coord_test = lab['control'] == b'CTRL_TRAJ_DATA'
            ind_test = (std_test & coord_test).nonzero()[0]
            if(sum(ind_test) != 0): # this will be zero for pre-R56 data
                self._scanner.ver = 'R56'
                
        
        self.chopON = False
        self.chopkzON = False
        self.selcoil = False
        self.cur_coil = -1
        self.selslice = False
        self.cur_loc = -1
        self.readParamOnly = False
        self.ApplyDownsampling = False
        self.dataDownsamplingFlag = False
        self.paddingWidth = 0
        self.averageWinWidth = 0
        self.hanningWidth = 0
        self.hanningWindow = 0
        self.retainedDatLen = 0
        self.propGRASE = 0
        self.propTSE = 0
        self.downsampleMode = 0
        self.rescale_type = 0
        self.raw_corr = 0
        self.trajtype = 0 # Use for user added orderings of trajectories
        self.delay = math.nan # Use for manual gr delays

                    

    def compute(self):

        # initialize data label
        dat_label = []

        header = dict()

        # check that the path actually exists
        if not os.path.exists(self.fname):
            return 0


        if self.readType == 0:
            if self.selcoil:
                cur_coil = -1 #self.getVal('Coil')
            else:
                cur_coil = -1  # read all the data

            if self.selslice:
                cur_loc = -1 #self.getVal('Slice')
            else:
                cur_loc = -1  # read all the data

            if os.path.exists(filename_extcase(self.base+'.data')) and os.path.exists(filename_extcase(self.base+'.list')) \
                    and self.readParamOnly is False:
                # read the data
                (dat, dat_noi, dat_phc, hdr, dat_label, noi_label, phc_label) \
                = readData(
                    filename_extcase(self.base+'.data'),
                    filename_extcase(self.base+'.list'),
                    self.chopON, self.propTSE, self.propGRASE, self.cur_coil, self.cur_loc)
                dat = dat.astype(np.float32)
                try:
                    dat_noi = dat_noi.astype(np.float32)
                except:
                    pass

                header['list'] = hdr
                header['headerType'] = 'list'
                # YCC start
                try:
                    dat_phc = dat_phc.astype(np.float32)
                except:
                    pass
                # YCC achieve
            elif self.readParamOnly is True:
                pass
            else:
                return 0

            # reset max of padding width for data downsampling
            if os.path.exists(self.base+'.data') and os.path.exists(self.base+'.list') \
                    and (self.readParamOnly is False) and (self.downsampleMode == 2):
                if self.paddingWidth > dat.shape[-2]//2:
                    self.paddingWidth = dat.shape[-2]//2
            else:
                if self.paddingWidth > 256:
                    self.paddingWidth = 64
                    
##############################################################################                    

        elif self.readType == 1:
            # read the data
            cur_loc = -1  # read all the data
            dat, hdr, dat_label = readRec(
                self.base, self.cur_loc, self.rescale_type)
            dat = dat.astype(np.float32)
            header = hdr
            
##############################################################################            

        elif self.readType == 2:
            if self.readParamOnly is False:
                hanningWidth = float(self.hanningWidth)/100.0
                hanningWindow = float(self.hanningWindow)/100.0

            if self.selcoil:
                cur_coil = -1 #self.getVal('Coil')
            else:
                cur_coil = -1  # read all the data

            if self.selslice:
                cur_loc = -1 #self.getVal('Slice')
            else:
                cur_loc = -1  # read all the data

            # read the data
            if os.path.exists(filename_extcase(self.base+'.lab')) and os.path.exists(filename_extcase(self.base+'.raw')) \
                    and os.path.exists(filename_extcase(self.base+'.sin')) \
                    and self.readParamOnly is False:
                (dat, dat_noi, dat_phc, rrs, rtops, hdr, dat_label, noi_label,
                 phc_label) = readRaw(
                    self.base, self.raw_corr, self.chopON, self.cur_coil, self.cur_loc)
                dat = dat.astype(np.complex64)
                try:
                    rrs = rrs.squeeze()
                    rtops = rtops.squeeze()
                    rrPlusRtop = np.stack((rrs, rtops))
                except:
                    pass
                try:
                    dat_phc = dat_phc.astype(np.complex64)
                except:
                    pass
                try:
                    dat_noi = dat_noi.astype(np.complex64)
                except:
                    pass
                header.update(hdr)

                # reset max of padding width for data downsampling
                if os.path.exists(self.base+'.data') and os.path.exists(self.base+'.list') \
                        and (self.readParamOnly is False) and (self.downsampleMode == 2):
                    if self.paddingWidth > dat.shape[-2]//2:
                        self.paddingWidth = dat.shape[-2]//2
                else:
                    if self.paddingWidth > 256:
                        self.paddingWidth = 64
            else:
                pass # read param only

        elif self.readType == 3:
            # read the data
            cur_coil = -1  # read all the data
            cur_loc = -1
            dat, hdr, dat_label = readCpx(
                self.base, cur_coil, cur_loc)
            dat = dat.astype(np.complex64)
            header = hdr
        # parm file

        if (self.readType in [0, 2]):
            if self._scanner.ver != 'R56':
                if os.path.exists(filename_extcase(self.base+".txt")):
                    # read the parm file
                    hdr = readParms(filename_extcase(self.base+".txt"))

                    # spOVERSAMPLING is set to either 1 or 2
                    if 'spOVERSAMPLING' in hdr:
                        if int(float(hdr['spOVERSAMPLING'][0])) == 2 and \
                                self.dataDownsamplingFlag is True:
                            doDownsampling = True

                    # This only supprt readType 0
                    if doDownsampling is True:
                        hdr['spDWELL'][0] = float(hdr['spDWELL'][0])*2
                        hdr['spREADPTS'][0] = float(hdr['spREADPTS'][0])//2
                        header['BNIspiral'] = hdr
                    else:
                        header = hdr
                elif os.path.exists(filename_extcase(self.base+".sin")) and \
                        os.path.exists(filename_extcase(self.base+".lab")):
                    if self.readParamOnly is True:
                        self.isMira = False
                        sin = readSin(filename_extcase(self.base+".sin"))
                        self.isMira = sin['isMira']
                        lab = readLab(filename_extcase(self.base+".lab"), self.isMira)
                        header = dict()
                        header['headerType'] = 'lab-sin'
                        header['sin'] = sin
                        header['lab'] = lab
                    else:
                        doDownsampling = self.dataDownsamplingFlag


            elif os.path.exists(filename_extcase(self.base+".sin")) and \
                    os.path.exists(filename_extcase(self.base+".lab")):
                if self.readParamOnly is True:
                    sin = readSin(filename_extcase(self.base+".sin"))
                    lab = readLab(filename_extcase(self.base+".lab"), self.isMira)
                    header = dict()
                    header['headerType'] = 'lab-sin'
                    header['sin'] = sin
                    header['lab'] = lab
                elif self.dataDownsamplingFlag is True:
                    doDownsampling = True
                    if 'sample_time_interval' in header:
                        header['sin']['sample_time_interval'][0][0] = float(header['sin']['sample_time_interval'][0][0])*2
                    if 'spiral_nr_grd_smpls_rmp_dn' in header:
                        header['sin']['spiral_nr_grd_smpls_rmp_dn'][0][0] = int(header['sin']['spiral_nr_grd_smpls_rmp_dn'][0][0])//2
                    if ('non_cart_min_enconding_nrs' in header) and ('non_cart_min_enconding_nrs' in header):
                        header['sin']['non_cart_min_encoding_nrs'][0][0] = int(header['sin']['non_cart_min_encoding_nrs'][0][0])//4
                        header['sin']['non_cart_max_encoding_nrs'][0][0] = \
                                int(header['sin']['non_cart_min_encoding_nrs'][0][0]) + \
                                int(header['sin']['non_cart_max_encoding_nrs'][0][0])//2

        # Convert from 2vec float to complex64
        if self.readType == 0 and self.readParamOnly is False:
            shape = list(dat.shape)
            shape.pop(-1)
            dat = np.frombuffer(dat.tobytes(), np.complex64)
            dat.shape = shape

            try:
                shape = list(dat_phc.shape)
                shape.pop(-1)
                dat_phc = np.frombuffer(dat_phc.tobytes(), np.complex64)
                dat_phc.shape = shape
            except:
                pass

            try:
                shape = list(dat_noi.shape)
                shape.pop(-1)
                dat_noi = np.frombuffer(dat_noi.tobytes(), np.complex64)
                dat_noi.shape = shape
            except:
                pass

        if self.readType in [0, 2] and self.chopkzON == 1 and self.readParamOnly is False:
            for i in range(dat.shape[-3]):
                if i % 2 == 0:
                    dat[..., i, :, :] *= -1.0

        # show some file stats
        fstats = os.stat(self.fname)
        # creation
        ctime = time.strftime('%d/%m/20%y',
                              time.localtime(fstats.st_ctime))
        # mod time
        mtime = time.strftime('%d/%m/20%y',
                              time.localtime(fstats.st_mtime))
        # access time
        atime = time.strftime('%d/%m/20%y',
                              time.localtime(fstats.st_atime))
        # filesize
        fsize = fstats.st_size
        # user id
        uid = fstats.st_uid
        # group id
        gid = fstats.st_gid

        if self.readParamOnly is False:
            d1 = list(dat.shape)
            info = "created: "+str(ctime)+"\n" \
                "accessed: "+str(atime)+"\n" \
                "modified: "+str(mtime)+"\n" \
                "UID: "+str(uid)+"\n" \
                "GID: "+str(gid)+"\n" \
                "file size (bytes): "+str(fsize)+"\n" \
                "Standard Data Info:\n" \
                "  dim labels: "+str(dat_label)+"\n" \
                "  dimensions: "+str(d1)+"\n" \
                "  type: "+str(dat.dtype)+"\n"
            try:
                d1 = list(dat_noi.shape)
                info = info+"Noise Data Info:\n" \
                    "  dim labels: "+str(noi_label)+"\n" \
                    "  dimensions: "+str(d1)+"\n" \
                    "  type: "+str(dat.dtype)+"\n"
            except:
                pass
            try:
                d1 = list(dat_phc.shape)
                info = info+"Phase Correction Data Info:\n" \
                    "  dim labels: "+str(phc_label)+"\n" \
                    "  dimensions: "+str(d1)+"\n" \
                    "  type: "+str(dat.dtype)+"\n"
            except:
                pass
            try:
                ht = hdr['headerType']
                info = info+"Header Type: '"+str(ht)+"'\n"
            except:
                pass
            self.info = info
        else:
            self.info = ''

        if self.readParamOnly is False:
            if self.dataDownsamplingFlag is True:
                dat_type = dat.dtype

                len_orig = dat.shape[-1]

                # prepare data for downsampling
                if self.downsampleMode == 0:  # FFT
                    pass
                elif self.downsampleMode == 1:  # extension + FFT
                    # generate the 1D Hanning window
                    # paddingWidth accounts for half of the additional data
                    # the full filter width is 4xpaddingWidth
                    radius = abs(2*np.linspace(-self.paddingWidth*2,
                                 (self.paddingWidth*4-1-self.paddingWidth*2),
                                 self.paddingWidth*4) /
                                 int(self.paddingWidth*4*hanningWidth))
                    hanningFilter = np.zeros(self.paddingWidth*4)
                    windIdx = radius <= 1.0
                    passIdx = radius <= (1.0 - hanningWindow)
                    func = 0.5 * (1.0 - np.cos(np.pi * (1.0 -
                                  radius[windIdx]) / hanningWindow))
                    hanningFilter[windIdx] = func
                    hanningFilter[passIdx] = 1.0
                    # take half of the filter
                    hanningFilter = hanningFilter[0:self.paddingWidth*2]

                    # save part of the raw Data (every other point at the
                    # begining and end used in the output
                    retainedDat1 = dat[..., 0:self.retainedDatLen*2:2]
                    retainedDat2 = dat[..., -self.retainedDatLen*2::2]
                    # extract data from beginning and end
                    datL = np.repeat(np.expand_dims(np.mean(
                        dat[..., -self.averageWinWidth:], axis=-1), axis=-1),
                        self.paddingWidth*2, axis=-1)
                    datLR = np.repeat(np.expand_dims(np.mean(
                        dat[..., 0:self.averageWinWidth], axis=-1), axis=-1),
                        self.paddingWidth*2, axis=-1)
                    # create data for padding
                    datLR = datLR-datL
                    datLR = datLR.reshape([-1, datL.shape[-1]])
                    # apply the hanning filter
                    for i in range(datLR.shape[0]):
                        datLR[i, :] = datLR[i, :]*hanningFilter
                    datLR = datLR.reshape(list(datL.shape)) + datL
                    # padding
                    dat = np.concatenate((dat, datLR), axis=-1)
                else:  # flipping padding
                    dat = np.concatenate((dat, dat[...,
                                         :-self.paddingWidth-1:-1]), axis=-1)
                    dat = np.concatenate((dat, dat[...,
                                         self.paddingWidth-1::-1]), axis=-1)

                # downsampling
                L = dat.shape[-1]//4
                R = L + dat.shape[-1]//2
                dat = np.fft.fftshift(np.fft.ifft(np.fft.ifftshift(dat,
                                      axes=-1), axis=-1), axes=-1). \
                    astype(dat_type)
                dat = np.fft.fftshift(np.fft.fft(np.fft.ifftshift(
                    dat[..., L:R], axes=-1), axis=-1), axes=-1). \
                    astype(dat_type)

                # extract the desired data
                if self.downsampleMode > 0:
                    dat = dat[..., 0:len_orig//2]
                if self.downsampleMode == 1:
                    dat[..., 0:self.retainedDatLen] = retainedDat1
                    dat[..., -self.retainedDatLen::] = retainedDat2

            self.data = dat
        try:
            self.noise = dat_noi
        except:
            pass
        try:
            self.phc = dat_phc
        except:
            pass
        try:
            self.rr_rtop = rrPlusRtop
        except:
            pass
        spparams = processSpiralParams(header, self.base, self._scanner)
        if len(spparams) > 1:
            self.spparams = spparams
        radparams = processRadialParams(header, self.trajtype, self.delay)
        if len(radparams) > 2:
            self.radparams = radparams
        try:
            self.header = header
        except:
            pass
        
        
