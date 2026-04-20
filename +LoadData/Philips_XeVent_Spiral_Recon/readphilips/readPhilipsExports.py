#!/usr/bin/env python
"""This module is a library of reader/parser routines for
the Philips export formats.
- raw/lab: earliest raw data type and label file
- list/data: raw data exported after several recon steps
  (eg. 1D gridding, phase correct, etc..)
- par/rec or xml/rec: reconstructed data and header info
"""

import os
import sys
import numpy as np
import re
import struct
from collections import OrderedDict


# match the correct filename case
#   * the existence check is case insensitive
#   * need the exact filename for opening a file
#   * the basename is chosen by the user so we just need to check for ext case
#   * this should be used when the file extension is guessed
def filename_extcase(fn):
    pn, ext = os.path.splitext(fn)
    bn = os.path.basename(pn)
    if bn+ext.lower() in os.listdir(os.path.dirname(fn)):
        return pn+ext.lower()
    elif bn+ext.upper() in os.listdir(os.path.dirname(fn)):
        return pn+ext.upper()
    return ''


# order preserving set generation
# pulled off the web (http://www.peterbe.com/plog/uniqifiers-benchmark)
def oset(seq, idfun=None):
    # order preserving
    if idfun is None:
        def idfun(x):
            if isinstance(x, (tuple, list)):
                return tuple(x)
            else:
                return(x)
    seen = {}
    result = []
    for item in seq:
        marker = idfun(item)
        # in old Python versions:
        # if seen.has_key(marker)
        # but in new ones:
        if marker in seen.keys():
            continue
        seen[marker] = 1
        result.append(item)
    return result


# sortd version of above function
def oset_sorted(seq, idfun=None):
    # order preserving
    if idfun is None:
        def idfun(x):
            if isinstance(x, (tuple, list)):
                return tuple(x)
            else:
                return(x)
    seen = {}
    result = []
    for item in seq:
        marker = idfun(item)
        # in old Python versions:
        # if seen.has_key(marker)
        # but in new ones:
        if marker in seen.keys():
            continue
        seen[marker] = 1
        result.append(item)
    result = sorted(result, key=lambda x: float(x))
    return result


# Support routines for .lab reading

# implemented from LabelTypeEnum.java
def LabelTypeEnum(i):

    if type(i) == int:
        hex_i = hex(int(i))

        if '0x7f00' == hex_i:
            return 'LABEL_TYPE_MIN'
        elif '0x7f01' == hex_i:
            return 'LABEL_TYPE_STANDARD'
        elif '0x7f02' == hex_i:
            return 'LABEL_TYPE_IGEO'
        elif '0x7f03' == hex_i:
            return 'LABEL_TYPE_IGEO_PMC'
        elif '0x7f04' == hex_i:
            return 'LABEL_TYPE_COIL_POS'
        elif '0x7f05' == hex_i:
            return 'LABEL_TYPE_NON_LIN'
        else:
            return 'LABEL_TYPE_MAX'

    if type(i) == np.ndarray:
        n = len(i)
        s = 'LABEL_TYPE_MAX ' * n
        label_type = np.array(s.split(), dtype='|S20')

        label_type[i == 0x7f00] = 'LABEL_TYPE_MIN'
        label_type[i == 0x7f01] = 'LABEL_TYPE_STANDARD'
        label_type[i == 0x7f02] = 'LABEL_TYPE_IGEO'
        label_type[i == 0x7f03] = 'LABEL_TYPE_IGEO_PMC'
        label_type[i == 0x7f04] = 'LABEL_TYPE_COIL_POS'
        label_type[i == 0x7f05] = 'LABEL_TYPE_NON_LIN'

        return label_type


def CtrlEnum(i):

    if type(i) == int:
        i = int(i)

        if i == -1:
            return 'CTRL_MIN'
        elif i == 0:
            return 'CTRL_NORMAL_DATA'
        elif i == 1:
            return 'CTRL_DC_OFFSET_DATA'
        elif i == 2:
            return 'CTRL_JUNK_DATA'
        elif i == 3:
            return 'CTRL_ECHO_PHASE_DATA'
        elif i == 4:
            return 'CTRL_NO_DATA'
        elif i == 5:
            return 'CTRL_NEXT_PHASE'
        elif i == 6:
            return 'CTRL_SUSPEND'
        elif i == 7:
            return 'CTRL_RESUME'
        elif i == 8:
            return 'CTRL_TOTAL_END'
        elif i == 9:
            return 'CTRL_INVALIDATION'
        elif i == 10:
            return 'CTRL_TYPE_NR_END'
        elif i == 11:
            return 'CTRL_VALIDATION'
        elif i == 12:
            return 'CTRL_NO_OPERATION'
        elif i == 13:
            return 'CTRL_DYN_SCAN_INFO'
        elif i == 14:
            return 'CTRL_SELECTIVE_END'
        elif i == 15:
            return 'CTRL_FRC_CH_DATA'
        elif i == 16:
            return 'CTRL_FRC_NOISE_DATA'
        elif i == 17:
            return 'CTRL_REFERENCE_DATA'
        elif i == 18:
            return 'CTRL_DC_FIXED_DATA'
        elif i == 19:
            return 'CTRL_NAVIGATOR_DATA'
        elif i == 20:
            return 'CTRL_FLUSH'
        elif i == 21:
            return 'CTRL_RECON_END'
        elif i == 22:
            return 'CTRL_IMAGE_STATUS'
        elif i == 23:
            return 'CTRL_TRACKING'
        elif i == 24:
            return 'CTRL_FLUOROSCOPY_TOGGLE'
        elif i == 25:
            return 'CTRL_REJECTED_DATA'
        elif i == 26:
            return 'CTRL_PROGRESS_INFO'
        elif i == 27:
            return 'CTRL_END_PREP_PHASE'
        elif i == 28:
            return 'CTRL_CHANNEL_DEFINITION'
        elif i == 29:
            return 'CTRL_START_SCAN_INFO'
        elif i == 30:
            return 'CTRL_DYN_PH_NAV_DATA'
        # for R56 DHW
        elif i == 31:
            return 'CTRL_TRAJ_DATA'
        else:
            return 'CTRL_MAX'

    if type(i) == np.ndarray:
        n = len(i)
        s = 'CTRL_MAX ' * n
        ctrl_type = np.array(s.split(), dtype='|S25')

        ctrl_type[i == -1] = 'CTRL_MIN'
        ctrl_type[i == 0] = 'CTRL_NORMAL_DATA'
        ctrl_type[i == 1] = 'CTRL_DC_OFFSET_DATA'
        ctrl_type[i == 2] = 'CTRL_JUNK_DATA'
        ctrl_type[i == 3] = 'CTRL_ECHO_PHASE_DATA'
        ctrl_type[i == 4] = 'CTRL_NO_DATA'
        ctrl_type[i == 5] = 'CTRL_NEXT_PHASE'
        ctrl_type[i == 6] = 'CTRL_SUSPEND'
        ctrl_type[i == 7] = 'CTRL_RESUME'
        ctrl_type[i == 8] = 'CTRL_TOTAL_END'
        ctrl_type[i == 9] = 'CTRL_INVALIDATION'
        ctrl_type[i == 10] = 'CTRL_TYPE_NR_END'
        ctrl_type[i == 11] = 'CTRL_VALIDATION'
        ctrl_type[i == 12] = 'CTRL_NO_OPERATION'
        ctrl_type[i == 13] = 'CTRL_DYN_SCAN_INFO'
        ctrl_type[i == 14] = 'CTRL_SELECTIVE_END'
        ctrl_type[i == 15] = 'CTRL_FRC_CH_DATA'
        ctrl_type[i == 16] = 'CTRL_FRC_NOISE_DATA'
        ctrl_type[i == 17] = 'CTRL_REFERENCE_DATA'
        ctrl_type[i == 18] = 'CTRL_DC_FIXED_DATA'
        ctrl_type[i == 19] = 'CTRL_NAVIGATOR_DATA'
        ctrl_type[i == 20] = 'CTRL_FLUSH'
        ctrl_type[i == 21] = 'CTRL_RECON_END'
        ctrl_type[i == 22] = 'CTRL_IMAGE_STATUS'
        ctrl_type[i == 23] = 'CTRL_TRACKING'
        ctrl_type[i == 24] = 'CTRL_FLUOROSCOPY_TOGGLE'
        ctrl_type[i == 25] = 'CTRL_REJECTED_DATA'
        ctrl_type[i == 26] = 'CTRL_PROGRESS_INFO'
        ctrl_type[i == 27] = 'CTRL_END_PREP_PHASE'
        ctrl_type[i == 28] = 'CTRL_CHANNEL_DEFINITION'
        ctrl_type[i == 29] = 'CTRL_START_SCAN_INFO'
        ctrl_type[i == 30] = 'CTRL_DYN_PH_NAV_DATA'
        ctrl_type[i == 31] = 'CTRL_TRAJ_DATA'

        return ctrl_type


# for parsing .list files
def readList(filename):
    '''
        Parse a *.list header accompanying the *.data data file.
    '''
    # open the list file
    if not type(filename) is str:
        print("Input filename is not a string.")
        sys.exit(1)

    if os.path.splitext(filename)[1] not in ['.list', '.LIST']:
        print("Input filename is not a .list file")
        sys.exit(1)

    # opens the file
    try:
        fil = open(filename_extcase(filename), 'r')
    except IOError:
        print('cannot open .list file ', filename)
        sys.exit(1)

    # read in all the lines at one time for processing
    lines = fil.readlines()
    fil.close()

    # data loc lines start with ' '
    loc = [line for line in lines if(line[0] == ' ')]
    loc = [line.strip() for line in loc]
    loc = [re.split(' +', line) for line in loc]

    # copied table header from a *.list file
    tab = "typ mix   dyn   card  echo  loca  chan  extr1 extr2 ky    kz" \
        + "    n.a.  aver  sign  rf    grad  enc   rtop  rr    size   offset"
    tab = re.split(' +', tab)

    info = dict()
    for label_idx in range(len(tab)):
        label = tab[label_idx]
        if label_idx < len(loc[0]):
            vals = [line[label_idx] for line in loc]
            info[label] = np.array(vals)

    # return header params and data locations
    return info


# read from the data file
def readData(filename_data, filename_list, chop_ky,
             prop_TSE, prop_GRASE, cur_coil, cur_loc):
    '''
        read from the data file using header data from the list file
    '''

    # header params from list file
    hdr = readList(filename_extcase(filename_list))

    # open the data file
    if not type(filename_data) is str:
        print("Input filename is not a string.")
        sys.exit(1)

    if os.path.splitext(filename_data)[1] not in ['.data', '.DATA']:
        print("Input filename is not a .data file")
        sys.exit(1)

    # opens the file
    try:
        fil = open(filename_extcase(filename_data), 'rb')
    except IOError:
        print('cannot open .data file ', filename_data)
        sys.exit(1)

    # generate ordered index sets
    kyvals = oset_sorted(hdr['ky'])
    kzvals = oset_sorted(hdr['kz'])
    mixes = oset_sorted(hdr['mix'])
    dynamics = oset_sorted(hdr['dyn'])
    phases = oset_sorted(hdr['card'])
    echoes = oset_sorted(hdr['echo'])
    extra1vals = oset_sorted(hdr['extr1'])
    extra2vals = oset_sorted(hdr['extr2'])
    grad_echoes = oset_sorted(hdr['grad'])
    rf_echoes = oset_sorted(hdr['rf'])
    signs = oset_sorted(hdr['sign'])

    # if a specific coil/loc is chosen to be read
    coil_case = np.ones(len(hdr['typ'])) > 0  # create array of True values
    loc_case = np.ones(len(hdr['typ'])) > 0  # create array of True values
    if cur_coil != -1:
        coil_case = hdr['chan'] == np.str(cur_coil)
    if cur_loc != -1:
        loc_case = hdr['loca'] == np.str(cur_loc)
    read_ind = (coil_case & loc_case).nonzero()[0]
    locations = oset_sorted(hdr['loca'][read_ind])

    # calculate the number of samples, channels, and averages for data output
    # standard, phc, and noise data can have different channels
    # and different numbers of samples
    std_case = hdr['typ'] == 'STD'
    std_ind = (std_case & coil_case & loc_case).nonzero()[0]
    noi_ind = (hdr['typ'] == 'NOI').nonzero()[0]
    phc_case = hdr['typ'] == 'PHC'
    phx_case = hdr['typ'] == 'PHX'
    phc_ind = (phc_case | phx_case).nonzero()[0]
    channels = oset_sorted(hdr['chan'][std_ind])
    noi_channels = oset_sorted(hdr['chan'][noi_ind])
    samples = oset_sorted(hdr['size'][std_ind])
    noi_samples = oset_sorted(hdr['size'][noi_ind])
    phc_samples = oset_sorted(hdr['size'][phc_ind])

    nchan = len(channels)
    nchan_noi = len(noi_channels)
    nsamp = max(np.array(samples, np.int64)) // 8
    if(len(noi_samples) > 0):
        nsamp_noi = max(np.array(noi_samples, np.int64)) // 8
    if len(phc_samples) > 0:
        nsamp_phc = max(np.array(phc_samples, np.int64)) // 8
    measurements = oset_sorted(hdr['aver'][std_ind])
    phc_measurements = oset_sorted(hdr['aver'][phc_ind])
    nmeas = len(measurements)
    nmeas_phc = len(phc_measurements)

    # number of each index
    nky = len(kyvals)
    nkz = len(kzvals)
    nmix = len(mixes)
    ndyn = len(dynamics)
    ncard = len(phases)
    necho = len(echoes)
    nloc = len(locations)
    nextra1 = len(extra1vals)
    nextra2 = len(extra2vals)
    ngrad = len(grad_echoes)
    nrfech = len(rf_echoes)

    # if ngrad is equal to necho, then we have multi-echo mDixon Propeller data
    if (ngrad == necho) & prop_TSE:
        ngrad = 1

    if nrfech > 1 and (prop_TSE or prop_GRASE):
        multi_shot = True
    else:
        multi_shot = False

    # allocate arrays
    # get the shape according to MPF loop order
    outshape = [nchan, nmix, ndyn, ncard, nextra1, nextra2, necho,
                nmeas, nloc, nkz, nky, nsamp, 2]
    outshape_string = np.array(['ch', 'mix', 'dyn', 'card', 'ex1', 'ex2',
                                'echo', 'meas', 'loc', 'kz', 'ky', 'samp'])

    # multi-shot TSE sequences (propeller)
    nshots = 1
    if multi_shot:

        # nblades = (total lines acquired) / (number of lines per blade)
        nshots = nky // nrfech // ngrad

        # final shape
        mshot_shape = [nchan, nmix, ndyn, ncard, nextra1, nextra2, necho,
                       nmeas, nloc, nkz, nrfech*ngrad, nshots, nsamp, 2]
        outshape_string = np.array(['ch', 'mix', 'dyn', 'card', 'ex1', 'ex2',
                                    'echo', 'meas', 'loc', 'kz', 'rf*grad',
                                    'shots', 'samp'])

    data_concat = np.zeros(outshape, dtype=np.float32)

    if len(phc_case.nonzero()[0]) > 0:
        outshape_phc = [nchan, nmix, ncard, necho, nmeas_phc, nloc, nkz,
                        nrfech, ngrad, nsamp_phc, 2]
        phc_string = np.array(['ch', 'mix', 'card', 'echo', 'meas', 'loc',
                              'kz', 'rf', 'grad', 'samp'])
        data_ph_concat = np.zeros(outshape_phc, dtype=np.float32)
    elif len(phx_case.nonzero()[0]) > 0:
        outshape_phx = [nchan, nmix, ncard, necho, nloc, nkz, nrfech, 2,
                        nsamp_phc, 2]
        phc_string = np.array(['ch', 'mix', 'card', 'echo', 'loc', 'pz', 'rf',
                               'sign', 'samp'])
        data_ph_concat = np.zeros(outshape_phx, dtype=np.float32)
    else:
        data_ph_concat = 0

    if(len(noi_samples) > 0):
        noi_shape = [nchan_noi, nsamp_noi, 2]
        noi_string = np.array(['ch', 'samp'])
        data_noi = np.zeros(noi_shape, dtype=np.float32)
    else:
        data_noi = 0

    # 1) get the data typ
    # 2) get the offset
    # 3) get the chunk size
    # read 4-byte float data
    for idx in read_ind:

        # load cur indices from sets
        ky = kyvals.index(hdr['ky'][idx])
        kz = kzvals.index(hdr['kz'][idx])
        mix = mixes.index(hdr['mix'][idx])
        dyn = dynamics.index(hdr['dyn'][idx])
        card = phases.index(hdr['card'][idx])
        echo = echoes.index(hdr['echo'][idx])
        loc = locations.index(hdr['loca'][idx])
        extr1 = extra1vals.index(hdr['extr1'][idx])
        extr2 = extra2vals.index(hdr['extr2'][idx])
        if hdr['typ'][idx] == 'PHC':
            meas = phc_measurements.index(hdr['aver'][idx])
        else:
            meas = measurements.index(hdr['aver'][idx])
        if hdr['typ'][idx] == 'NOI':
            noi_ch = noi_channels.index(hdr['chan'][idx])
        else:
            ch = channels.index(hdr['chan'][idx])
        sign = signs.index(hdr['sign'][idx])
        rf = rf_echoes.index(hdr['rf'][idx])
        grad = grad_echoes.index(hdr['grad'][idx])

        # get info of current segment
        offset = int(hdr['offset'][idx])
        size = int(hdr['size'][idx])
        samples = size // 8  # divide by complex float size

        # jump to offset
        fil.seek(offset)

        # read the data into an NPY
        data = fil.read(samples*8)
        data = np.frombuffer(data, dtype=np.float32)

        if hdr['typ'][idx] == 'STD':
            # apply chop if appropriate
            if chop_ky:
                data.shape = [samples, 2]
                if (int(hdr['ky'][idx]) % 2) == 0:
                    data *= -1.0
                data.shape = [samples*2]

        # load segment into array
        try:
            data.shape = [samples, 2]
            if hdr['typ'][idx] == 'NOI':
                data_noi[noi_ch, 0:samples, :] = data
            elif hdr['typ'][idx] == 'STD':
                data_concat[ch, mix, dyn, card, extr1, extr2, echo, meas, loc,
                            kz, ky, 0:samples, :] = data
            elif hdr['typ'][idx] == 'PHC':
                data_ph_concat[ch, mix, card, echo, meas, loc, kz, rf, grad,
                               0:samples, :] = data
            elif hdr['typ'][idx] == 'PHX':
                data_ph_concat[ch, mix, card, echo, loc, kz, rf, sign,
                               0:samples, :] = data

        except:
            # if an error occurs then print out the indexing
            # information, then raise exception
            print("\n\tERROR: index out of range:")
            if hdr['typ'][idx] == 'NOI':
                print("\tdata_noi->dimensions:"+str(data_noi.shape))
                print("\tcurrent channel  :"+str(noi_ch))
            if hdr['typ'][idx] == 'STD':
                print("\tdata->dimensions:"+str(data_concat.shape))
                print("\tcurrent index   :"+str([ch, mix, dyn, card, extr1,
                                                extr2, echo, meas, loc, kz,
                                                ky]))
            elif hdr['typ'][idx] == 'PHC':
                print("\tdata->dimensions:"+str(data_ph_concat.shape))
                print("\tcurrent index   :"+str([ch, mix, card, echo, meas,
                                                loc, kz, rf, grad]))
            elif hdr['typ'][idx] == 'PHX':
                print("\tdata->dimensions:"+str(data_ph_concat.shape))
                print("\tcurrent index   :"+str([ch, mix, card, echo, loc, kz,
                                                rf, sign]))
            print("\tnshots:"+str(nshots))
            print("\tmulti_shot:"+str(multi_shot))
            print("\tnsamp[]:"+str(nsamp))
            print(" ")
            raise

    # close the lab binary file
    fil.close()

    # reshape to include shots
    if multi_shot:
        print("\treadData(): transpose shots and rf echos (slots 10 & 11)...")
        data_concat.shape = mshot_shape
        data_concat = data_concat.transpose((
            0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 10, 12, 13))

    # setup data label strings
    data_labels = outshape_string[(np.array(
        data_concat.shape[0:len(outshape_string)]) > 1).nonzero()[0]]
    if len(phc_ind) > 0:
        phc_labels = phc_string[(np.array(
            data_ph_concat.shape[0:len(phc_string)]) > 1).nonzero()[0]]
    else:
        phc_labels = 0
    if(len(noi_samples) > 0):
        noi_labels = noi_string[(np.array(
            data_noi.shape[0:len(noi_string)]) > 1).nonzero()[0]]
    else:
        noi_labels = 0

    # return dictionary
    return(data_concat, data_noi, data_ph_concat, hdr,
           data_labels, noi_labels, phc_labels)


# for parsing .par files
def readPar(filename):
    '''
        Parse a *.PAR header accompanying the *.REC image file.
    '''

    filename = filename_extcase(filename)

    # open the par file
    if not type(filename) is str:
        print("Input filename is not a string.")
        sys.exit(1)

    if os.path.splitext(filename)[1] not in ['.PAR', '.par']:
        print("Input filename is not a .PAR file")
        sys.exit(1)

    # opens the file
    try:
        fil = open(filename, 'r')
    except IOError:
        print('cannot open .par file ', filename)
        sys.exit(1)

    # read in all the lines at one time for processing
    lines = fil.readlines()
    fil.close()

    # data loc lines do not start with '#' or '.'
    loc = [line for line in lines if(line[0] not in ['#', '.'] and
                                     len(line) > 2)]
    loc = [line.strip() for line in loc]
    loc = [re.split(' +', line) for line in loc]

    # modified table header from a V4.2 *.par file
    tab = "Slice,Echo,Dynamic,Phase,Type,Sequence,Index,Pixel Size,"\
        + "Scan Percentage,Resolution X,Resolution Y,Rescale Intercept,"\
        + "Rescale Slope,Scale Slope,Window Center,Window Width,"\
        + "Angulation AP,Angulation FH,Angulation RL,Offcenter AP,"\
        + "Offcenter FH,Offcenter RL,Slice Thickness,Slice Gap,"\
        + "Display Orientation,Slice Orientation,fMRI Status Indication,"\
        + "Image Type Ed Es,Pixel Spacing_0,Pixel Spacing_1,Echo Time,"\
        + "Dyn Scan Begin Time,Trigger Time,Diffusion B Factor,No Averages,"\
        + "Image Flip Angle,Cardiac Frequency,Min RR Interval,"\
        + "Max RR Interval,TURBO Factor,Inversion Delay,BValue,Grad Orient,"\
        + "Contrast Type,Diffusion Anisotropy Type,Diffusion AP,Diffusion FH,"\
        + "Diffusion RL,Label Type"
    tab = re.split(',', tab)

    info = dict()
    info['headerType'] = '.par'
    for label_idx in range(len(tab)):
        label = tab[label_idx]
        if label_idx < len(loc[0]):
            vals = [line[label_idx] for line in loc]
            info[label] = np.array(vals)

    # return header params and data locations
    return info


# for parsing *.xml files
def readXML(filename):
    '''
        Parse a *.xml header accompanying the *.rec image file.
    '''
    import xml.etree.ElementTree as et

    filename = filename_extcase(filename)

    # open the xml file
    if not type(filename) is str:
        print("Input filename is not a string.")
        sys.exit(1)

    if os.path.splitext(filename)[1] not in ['.XML', '.xml']:
        print("Input filename is not a .XML file")
        sys.exit(1)

    # opens the file
    try:
        tree = et.parse(filename)
        root = tree.getroot()
    except IOError:
        print('cannot open', filename)

    # get table headers directly from xml file
    tab = [line.attrib['Name'] for line in root[1][0].iter('Attribute')]

    text = str([line.text for line in root[1].iter('Attribute')]).strip('[]')
    text = re.split(', +', text)
    text = [line.strip('\'') for line in text]
    try:
        loc = list(zip(*[iter(text)]*len(tab)))
    except:
        print("coding error: defined table length not equal to number of "
              "image attributes")

    info = dict()
    info['headerType'] = '.xml'
    for label_idx in range(len(tab)):
        label = tab[label_idx]
        if label_idx < len(loc[0]):
            vals = [line[label_idx] for line in loc]
            info[label] = np.array(vals)

    # return header params and data locations
    return info


# read from the *.rec file
def readRec(filename, cur_loc, rescale_type):
    '''
        read from the REC file using header data from the PAR or XML file
    '''

    filename_par = filename_extcase(filename+'.PAR')
    filename_xml = filename_extcase(filename+'.XML')
    filename_rec = filename_extcase(filename+'.REC')
    headerType = 0

    if os.path.exists(filename_rec) and os.path.exists(filename_xml):
        # header params from xml file
        hdr = readXML(filename_xml)
        headerType = 0
    elif os.path.exists(filename_rec) and os.path.exists(filename_par):
        # header params from PAR file
        hdr = readPar(filename_par)
        headerType = 1
    else:
        print("ReadPhilips():readRec():ERROR: PAR, XML, and/or REC path"
              "doesn't exist.")
        return 1  # WARNING

    # open the REC file
    if not type(filename_rec) is str:
        print("Input filename is not a string.")
        sys.exit(1)

    if os.path.splitext(filename_rec)[1] not in ['.REC', '.rec']:
        print("Input filename is not a .REC file")
        sys.exit(1)

    # opens the file
    try:
        fil = open(filename_rec, 'rb')
    except IOError:
        print('cannot open .rec file ', filename_rec)
        sys.exit(1)

    # Find number of unique entries for each label.
    # The PAR header is not a reliable source of this information so find
    # it by calculating the number of unique values in the PAR or XML file.
    # we want most of these sorted in ascending order
    dynamics = oset_sorted(hdr['Dynamic'])
    phases = oset_sorted(hdr['Phase'])
    echoes = oset_sorted(hdr['Echo'])

    loc_case = np.ones(len(hdr['Slice'])) > 0  # create array of True values
    # if a specific loc is chosen to be read
    if cur_loc != -1:
        loc_case = hdr['Slice'] == np.str(cur_loc)
    read_ind = loc_case.nonzero()[0]
    locations = oset_sorted(hdr['Slice'][read_ind])

    imtypes = oset(hdr['Type'])
    seqtypes = oset(hdr['Sequence'])
    ndyn = len(dynamics)
    ncard = len(phases)
    necho = len(echoes)
    nlocs = len(locations)
    nimtype = len(oset(imtypes))
    nsqtype = len(oset(seqtypes))

    # calculate the number of samples, channels, and averages for data output
    # initialize PAR version to '4'
    version = '4'

    # try to read in diffusion and ASL information
    # determine the PAR version accordingly
    try:
        diffvalues = oset_sorted(hdr['BValue'])
        gradorients = oset_sorted(hdr['Grad Orient'])
        ndiff = len(diffvalues)
        ngrad = len(gradorients)
        # must be at least version '4.1'
        version = '4.1'
    except:
        pass

    try:
        labtypes = oset(hdr['Label Type'])
        nlabel = len(labtypes)
        # must be version '4.2'
        version = '4.2'
    except:
        pass

    # Sometimes images are appended to the REC file (e.g. T2 maps,
    # B0 maps, etc.). We can only be certain that the total number of images
    # is divisible by the dynamics, cardiac phases, and slice locations. This
    # variable keeps track of the number of different image contrasts.
    ncontrasts = int(np.ceil(len(hdr['Index']) / float(ndyn * ncard * nlocs)))

    # Images are sometimes appended to the set.
    # Need to flatten slower dimensions in these cases.
    expected_contrasts = necho * nimtype * nsqtype
    if (headerType == 0 or version in ['4.1', '4.2']):
        expected_contrasts = expected_contrasts * ndiff * ngrad
    if (headerType == 0 or version == '4.2'):
        expected_contrasts = expected_contrasts * nlabel
    flatten = False
    if (expected_contrasts != ncontrasts or ncontrasts % expected_contrasts
            != 0):
        flatten = True

    # set up ordered list of unique sets of contrasts, diffusion weightings....
    if flatten:
        if (headerType == 0 or version == '4.2'):
            contrasts = list(map(list, zip(*[hdr['Echo'], hdr['Type'],
                                             hdr['Sequence'], hdr['BValue'],
                                             hdr['Grad Orient'],
                                             hdr['Label Type']])))
        elif (version == '4.1'):
            contrasts = list(map(list, zip(*[hdr['Echo'], hdr['Type'],
                                             hdr['Sequence'], hdr['BValue'],
                                             hdr['Grad Orient']])))
        else:
            contrasts = list(map(list, zip(*[hdr['Echo'], hdr['Type'],
                                             hdr['Sequence']])))
        contrasts = oset(contrasts)

    # get the image dimensions
    recx = int(hdr['Resolution X'][0])
    recy = int(hdr['Resolution Y'][0])

    # set shape of output data
    if flatten:
        outshape = [ncontrasts, ndyn, ncard, nlocs, recy, recx]
        data_string = np.array(['cntrst', 'dyn', 'card', 'loc', 'x', 'y'])
    elif headerType == 0 or version == '4.2':
        outshape = [nimtype, nsqtype, nlabel, ngrad, ndiff, ndyn, ncard,
                    necho, nlocs, recy, recx]
        data_string = np.array(['imtype', 'seq', 'lab', 'grad', 'bval', 'dyn',
                                'card', 'echo', 'loc', 'x', 'y'])
    elif version == '4.1':
        outshape = [nimtype, nsqtype, ngrad, ndiff, ndyn, ncard, necho, nlocs,
                    recy, recx]
        data_string = np.array(['imtype', 'seq', 'grad', 'bval', 'dyn', 'card',
                                'echo', 'loc', 'x', 'y'])
    else:
        outshape = [nimtype, nsqtype, ndyn, ncard, necho, nlocs, recy, recx]
        data_string = np.array(['imtype', 'seq', 'dyn', 'card', 'echo', 'loc',
                                'x', 'y'])

    # allocate output array
    data_concat = np.zeros(outshape, dtype=np.float32)

    # read the REC file data for each image index
    for idx in hdr['Index'][read_ind]:

        # load values from dictionary
        idx = int(idx)
        dyn = dynamics.index(hdr['Dynamic'][idx])
        card = phases.index(hdr['Phase'][idx])
        echo = echoes.index(hdr['Echo'][idx])
        loc = locations.index(hdr['Slice'][idx])
        if headerType == 0 or version in ['4.1', '4.2']:
            bval = diffvalues.index(hdr['BValue'][idx])
            grad = gradorients.index(hdr['Grad Orient'][idx])
        if headerType == 0 or version == '4.2':
            lab = labtypes.index(hdr['Label Type'][idx])

        # get index of the current image and sequence type
        typ = imtypes.index(hdr['Type'][idx])
        seq = seqtypes.index(hdr['Sequence'][idx])

        # data type for image data
        pixel_bytes = int(hdr['Pixel Size'][idx]) // 8
        data_type = np.uint16
        if pixel_bytes == 1:
            data_type = np.uint8
        elif pixel_bytes == 2:
            data_type = np.uint16
        size_bytes = recx * recy * pixel_bytes

        # get amplitude scaling information
        ri = float(hdr['Rescale Intercept'][idx])
        rs = float(hdr['Rescale Slope'][idx])
        ss = float(hdr['Scale Slope'][idx])

        # contrast number for flattened output
        # TODO: verify the looping order
        if flatten:
            current_contrast = [hdr['Echo'][idx], hdr['Type'][idx],
                                hdr['Sequence'][idx]]
            if headerType == 0 or version in ['4.1', '4.2']:
                current_contrast.append(hdr['BValue'][idx])
                current_contrast.append(hdr['Grad Orient'][idx])
            if headerType == 0 or version == '4.2':
                current_contrast.append(hdr['Label Type'][idx])
            cont = contrasts.index(current_contrast)
            if ncontrasts != len(contrasts):
                cont = int(idx // (ndyn * ncard)) % ncontrasts

        # get offset of current image
        offset = int(idx * size_bytes)

        # jump to offset
        fil.seek(offset)

        # read the data into an NPY array
        data = fil.read(size_bytes)
        data = np.frombuffer(data, dtype=data_type)

        # Re-scale data according to selected option
        if rescale_type == 0:
            data = (data * rs + ri) / (rs * ss)  # floating point
        elif rescale_type == 1:
            data = (data * rs + ri)  # displayed value from scanner
        # else no scaling (use the stored value in the REC file)

        # load image into array
        try:
            data.shape = [recx, recy]
            if flatten:
                data_concat[cont, dyn, card, loc, 0:recx, 0:recy] = data
            elif headerType == 0 or version == '4.2':
                data_concat[typ, seq, lab, grad, bval, dyn, card, echo, loc,
                            0:recx, 0:recy] = data
            elif version == '4.1':
                data_concat[typ, seq, grad, bval, dyn, card, echo, loc,
                            0:recx, 0:recy] = data
            else:
                data_concat[typ, seq, dyn, card, echo, loc, 0:recx,
                            0:recy] = data

        except:
            # if an error occurs then print out the indexing
            # information, then raise exception
            print("\n\tERROR: index out of range:")
            print("\tdata->dimensions:"+str(data_concat.shape))
            if flatten:
                print("\tcurrent index   :"+str([cont, dyn, card, loc]))
            elif headerType == 0 or headerType == '4.2':
                print("\tcurrent index   :"+str([typ, seq, lab, grad, bval,
                                                dyn, card, echo, loc]))
            elif headerType == '4.1':
                print("\tcurrent index   :"+str([typ, seq, grad, bval, dyn,
                                                card, echo, loc]))
            else:
                print("\tcurrent index   :"+str([typ, seq, dyn, card, echo,
                                                loc]))
            print(" ")
            raise

    # close the rec binary file
    fil.close()

    # setup the data labels
    data_labels = data_string[(np.array(
        data_concat.shape[0:len(data_string)]) > 1).nonzero()[0]]

    # return dictionary
    return(data_concat, hdr, data_labels)


# label file parsing
def readLab(filename, isMira):
    '''
        Parse a *.lab header accompanying the *.raw data file.
    '''

    # open the lab file
    if not type(filename) is str:
        print("Input filename is not a string.")
        sys.exit(1)

    if os.path.splitext(filename)[1] not in ['.lab', '.LAB']:
        print("input filename is not a .lab file")
        sys.exit(1)

    # opens the file
    try:
        fil = open(filename_extcase(filename), 'rb')
    except IOError:
        print('cannot open lab file', filename)
        sys.exit(1)

    LABEL_SIZE = 64  # bytes
    lab = dict()
    lab_index = 0

    # first calculate the number of useful labels and get the non-linear labels
    fmt1 = '<'+'2If4H22h'  # non-lin label struct 
    num_labels = 0
    channel_nr = []
    nr_breakpoints = []
    nl_indices = []
    while True:
        fil.seek(lab_index)
        line = fil.read(LABEL_SIZE)

        # make sure label line exists
        if len(line) == LABEL_SIZE:
            label = struct.unpack(fmt1, line)
        else:
            break

        # only read the valid label lines
        if (LabelTypeEnum(label[4]) == 'LABEL_TYPE_MAX' or
                (LabelTypeEnum(label[4]) == 'LABEL_TYPE_STANDARD' and
                    CtrlEnum(label[5]) == 'CTRL_TOTAL_END')):
            break
        else:
            num_labels += 1
            lab_index += LABEL_SIZE
            if (LabelTypeEnum(label[4]) == 'LABEL_TYPE_NON_LIN'):
                channel_nr.append(label[5])
                nr_breakpoints.append(label[6])
                nl_indices.append(num_labels)
    lab['non_lin_info'] = dict()
    lab['non_lin_info']['channel_nr'] = np.array(channel_nr).astype(np.uint16)
    lab['non_lin_info']['nr_breakpoints'] = np.array(nr_breakpoints).astype(np.uint16)

    # read in the useful part of the file
    fil.seek(0)
    unparsed_labels = fil.read(num_labels * LABEL_SIZE)

    # this format string is used to parse the entire .lab file
    if isMira:
        fmt1 = '<'+'2If2H6B19hi' * num_labels
    else:
        fmt1 = '<'+'i6h6B19hi' * num_labels
    parsed_labels = struct.unpack(fmt1, unparsed_labels[0:
                                  LABEL_SIZE*num_labels])
    labelData = np.array(parsed_labels)
    if isMira:
        labelData.shape = (num_labels, 31)
    else:
        labelData.shape = (num_labels, 33)

    lab['data_size'] = labelData[:, 0].astype(np.uint32)
    if isMira:
        lab['coded_data_size'] = labelData[:, 1].astype(np.uint32)
        lab['normalization_factor'] = labelData[:, 2]
        baseind = 3
    else:
        lab['leading_dummies_size'] = labelData[:, 1]
        lab['trailing_dummies_size'] = labelData[:, 2]
        lab['src_code'] = labelData[:, 3]
        lab['dst_code'] = labelData[:, 4]
        baseind = 5

    lab['seq_nr'] = labelData[:, baseind].astype(np.uint16)
    lab['label_type'] = LabelTypeEnum(labelData[:, baseind+1])
    lab['control'] = CtrlEnum(labelData[:, baseind+2])
    lab['monitoring_flag'] = labelData[:, baseind+3].astype(np.uint8)
    lab['measurement_phase'] = labelData[:, baseind+4].astype(np.uint8)
    lab['measurement_sign'] = labelData[:, baseind+5].astype(np.uint8)
    lab['gain_setting_index'] = labelData[:, baseind+6].astype(np.uint8)
    lab['raw_format'] = labelData[:, baseind+7].astype(np.uint8)
    lab['spur_phase_increment'] = labelData[:, baseind+8].astype(np.uint16)
    lab['progress_cnt'] = labelData[:, baseind+9].astype(np.int16)
    lab['mix_nr'] = labelData[:, baseind+10].astype(np.uint16)
    lab['dynamic_scan_nr'] = labelData[:, baseind+11].astype(np.uint16)
    lab['cardiac_phase_nr'] = labelData[:, baseind+12].astype(np.uint16)
    lab['echo_nr'] = labelData[:, baseind+13].astype(np.uint16)
    lab['location_nr'] = labelData[:, baseind+14].astype(np.uint16)
    lab['row_nr'] = labelData[:, baseind+15].astype(np.uint16)
    lab['extra_attr_nr'] = labelData[:, baseind+16].astype(np.uint16)
    lab['measurement_nr'] = labelData[:, baseind+17].astype(np.uint16)
    lab['e1_profile_nr'] = labelData[:, baseind+18].astype(np.uint16)
    lab['e2_profile_nr'] = labelData[:, baseind+19].astype(np.uint16)
    lab['e3_profile_nr'] = labelData[:, baseind+20].astype(np.uint16)
    lab['rf_echo_nr'] = labelData[:, baseind+21].astype(np.uint16)
    lab['grad_echo_nr'] = labelData[:, baseind+22].astype(np.uint16)
    lab['enc_time'] = labelData[:, baseind+23].astype(np.int16)
    lab['random_phase'] = labelData[:, baseind+24].astype(np.int16)
    lab['rr_interval'] = labelData[:, baseind+25].astype(np.uint16)
    lab['rtop_offset'] = labelData[:, baseind+26].astype(np.uint16)
    lab['channels_active'] = labelData[:, baseind+27]

    # close the lab binary file
    fil.close()

    # return dictionary
    return(lab)


# read/parse .sin text file
def readSin(filename):
    '''
        Parse a *.sin header accompanying the *.raw data file.
    '''

    # open the sin file
    if not type(filename) is str:
        print("Input filename is not a string.")
        sys.exit(1)

    # .sin file check
    if os.path.splitext(filename)[1] not in ['.sin', '.SIN']:
        print("input filename is not a .sin file")
        sys.exit(1)

    # opens the file if it exists
    try:
        fil = open(filename, 'r')
    except IOError:
        print('cannot open', filename)
        sys.exit(1)

    # read in all the lines at one time for processing
    lines = fil.readlines()
    fil.close()

    # remove lines w/o colon delimiters
    lines = [line for line in lines if(line.find(' : ') != -1)]

    # split on colons first
    #   -strip out colons at the end of lines
    lines = [re.split(' +: +|: +| +:|:', line.strip(
        ':\r\n').strip()) for line in lines]

    # throw out single sub-array elements
    lines = [line for line in lines if(len(line) > 1)]

    # load dictionary based on second element (ie. the item name)
    sin = OrderedDict()
    for line in lines:
        # if key already exists, extend rather than overwrite
        try:
            sin[line[1]][0].extend(re.split(' +', line[2]))
            sin[line[1]][1].extend([re.split(' +', line[0])])
        except:
            sin[line[1]] = [re.split(' +', line[2]), [re.split(' +', line[0])]]

    # some logic from RawLabData.java

    # frequency encode samples:
    try:
        sin['nr_samples'] = int(sin['max_encoding_numbers'][0][0]) \
            - int(sin['min_encoding_numbers'][0][0]) + 1
    except:
        print("no max_encoding_numbers found, skipping")

    # if dc_max_encoding_number exists, then recalc
    try:
        sin['nr_samples'] = int(sin['dc_max_encoding_number'][0][0]) \
            - int(sin['dc_min_encoding_number'][0][0]) + 1
    except:
        print("no dc_max_encoding_number found, skipping")

    # if dr_max_encoding_numbers exists, then recalc
    try:
        sin['nr_samples'] = int(sin['dr_max_encoding_numbers'][0][0]) \
            - int(sin['dr_min_encoding_numbers'][0][0]) + 1
        sin['spectro'] = True
    except:
        print("no dr_max_encoding_numbers found, skipping")
        sin['spectro'] = False

    # if non_cart_max_encoding_numbers exists, then recalc
    try:
        sin['nr_samples'] = int(sin['non_cart_max_encoding_nrs'][0][0]) \
            - int(sin['non_cart_min_encoding_nrs'][0][0]) + 1
    except:
        print("no non_cart_max_encoding_nrs found, skipping")

    # use channel_names to determine if data is Mira
    if any("DCC" or "MN" in name for name in sin['channel_names'][0]):
        sin['isMira'] = True
    else:
        #R56 FLORET - start - use previous Mira logic for now
        try:
            if int(sin['spiral_trajectory_shape'][0][0]) == 3: #FLORET
                try:
                    sin['relative_fear_bandwidth']
                    sin['isMira'] = True
                except:
                # if enable_pda exists then this is not a Mira file
                    try:
                        sin['enable_pda']
                        sin['isMira'] = False
                    except:
                        sin['isMira'] = True
            else:
                sin['isMira'] = False
        except:
            #R56 FLORET -end
            sin['isMira'] = False


    # try to determine if data is Mira (logic not correct for pre R5 non-Cart)
    # try:
    #     sin['relative_fear_bandwidth']
    #     sin['isMira'] = True
    # except:
    #     # if enable_pda exists then this is not a Mira file
    #     try:
    #         sin['enable_pda']
    #         sin['isMira'] = False
    #     except:
    #         sin['isMira'] = True

    # phase encode lines
    sin['nr_e1_profiles'] = int(sin['max_encoding_numbers'][0][1]) \
        - int(sin['min_encoding_numbers'][0][1]) + 1
    sin['nr_e2_profiles'] = int(sin['max_encoding_numbers'][0][2]) \
        - int(sin['min_encoding_numbers'][0][2]) + 1
    sin['nr_e3_profiles'] = int(sin['max_encoding_numbers'][0][3]) \
        - int(sin['min_encoding_numbers'][0][3]) + 1

    return(sin)


# set up the general Philips raw -> list data corrections
def setupCorrections(fil_raw, lab, sin, nchan, isMira):
    if 'enable_dc_offset_corr' in sin.keys():
        enable_dc_corr = bool(sin['enable_dc_offset_corr'][0])
    else:
        enable_dc_corr = False
    basic_corr = dict()

    try:
        coil_factor = float(sin['scale_per_channel_arr'][0][1])
        coil_factor *= np.exp(-2j * float(sin['phase_per_channel_arr'][0][1]))
    except:
        coil_factor = 1.0
    basic_corr['coil_factor'] = coil_factor

    try:
        pda_factors = np.array(sin['pda_ampl_factors'])
        pda_factors.shape = [nchan, 24]
    except:
        pda_factors = np.ones([nchan, 24])
    basic_corr['pda_factors'] = pda_factors

    # if R5 or later
    try:
        relative_fear_bandwidth = float(sin['relative_fear_bandwidth'][0][0])
    except:
        relative_fear_bandwidth = 0.0
    basic_corr['relative_fear_bandwidth'] = relative_fear_bandwidth

    basic_corr['nchan'] = nchan

    # get non-linearity correction values
    channels = lab['non_lin_info']['channel_nr']
    nr_breakpoints = lab['non_lin_info']['nr_breakpoints']
    try:
        # find indices and seek offsets for non-linearity data
        seek_offset = 512
        nl_case = lab['label_type'] == b'LABEL_TYPE_NON_LIN'
        nl_ind = (nl_case).nonzero()[0][0]

        format_vals = lab['raw_format'][0:nl_ind] == 4
        format_vals |= lab['raw_format'][0:nl_ind] == 6
        if isMira:
            seek_val = np.sum(np.where(format_vals, 
                                   lab['coded_data_size'][0:nl_ind],
                                   lab['data_size'][0:nl_ind]))
        else:
            seek_val = np.sum(lab['data_size'][0:nl_ind])
        seek_offset = int(seek_val + seek_offset)
        for channel in channels:
            fil_raw.seek(seek_offset)
            nbp_chan = nr_breakpoints[channel]
            data_size = int(12 * nbp_chan)
            bytebuff = fil_raw.read(data_size)
            bp_raw = np.frombuffer(bytebuff, dtype=np.float32)
            bp_raw.shape = [nbp_chan, 3]
            # pad breakpoints by 2 to account for levels above and below range
            breakpoints = np.zeros([nbp_chan + 2, 3])
            breakpoints[1:nbp_chan+1, :] = bp_raw
            breakpoints[0, 0] = 0.
            breakpoints[nbp_chan+1, 0] = 10**10
            breakpoints[0, 1::] = breakpoints[1, 1::]
            breakpoints[nbp_chan+1, 1::] = breakpoints[nbp_chan, 1::]
            basic_corr['bp_chan'+str(channel)] = breakpoints
            seek_offset += data_size
    except:
        print('No RXE Non-linear correction data present.')

    if enable_dc_corr:
        # find indices and seek offsets for dc fixed data
        std_case = lab['label_type'] == b'LABEL_TYPE_STANDARD'
        dc_fixed_case = lab['control'] == b'CTRL_DC_FIXED_DATA'
        dc_ind = (std_case & dc_fixed_case).nonzero()[0][0]

        # find seek position in raw file
        seek_offset = 512
        format_vals = lab['raw_format'][0:dc_ind] == 4
        format_vals |= lab['raw_format'][0:dc_ind] == 6
        if isMira:
            seek_val = np.sum(np.where(format_vals, 
                                   lab['coded_data_size'][0:dc_ind],
                                   lab['data_size'][0:dc_ind]))
        else:
            seek_val = np.sum(lab['data_size'][0:dc_ind])
        seek_offset = int(seek_val + seek_offset)
        
        # read dc fixed data
        data_size = int(lab['data_size'][dc_ind])
        fil_raw.seek(seek_offset)
        bytebuff = fil_raw.read(data_size)
        temp_data = np.frombuffer(bytebuff, dtype=np.int32)
        dc_fixed_vals = temp_data.shape[0] // 2
        temp_data.shape = [dc_fixed_vals, 2]
        dc_fixed_arr = temp_data[0:nchan, 0] + 1j*temp_data[0:nchan, 1]

        # normalize dc fixed array
        try:
            norm_factor = lab['normalization_factor'][dc_ind]
            norm_factor = norm_factor / ((2**15 - 1) * 10000.)
        except:
            norm_factor = 1. / ((2**15 - 1) * 10000.)
        dc_fixed_arr *= norm_factor
        basic_corr['dc_fixed_arr'] = dc_fixed_arr

    return basic_corr


# apply the basic Philips raw -> list data corrections
def basicCorrections(data, lab, sin, lab_pos, basic_corr, isMira):
    if 'enable_dc_offset_corr' in sin.keys():
        enable_dc_corr = bool(sin['enable_dc_offset_corr'][0])
    else:
        enable_dc_corr = False

    # general correction terms
    coil_factor = basic_corr['coil_factor']
    pda_factors = basic_corr['pda_factors']
    relative_fear_bandwidth = basic_corr['relative_fear_bandwidth']
    nchan = basic_corr['nchan']

    # apply RXE non-linearity correction
    try:
        norm_factor = int(lab['normalization_factor'][lab_pos])
    except:
        norm_factor = 1.
    for chan in range(data.shape[0]):
        data[chan, ...] *= norm_factor
        # get breakpoints for current channel
        try:
            breakpoints = basic_corr['bp_chan'+str(chan)]
            if breakpoints.shape[0] > 1:
                # bin the data based on the breakpoint levels
                bins = breakpoints[:, 0]
                x = np.abs(data[chan, ...])
                binned = np.digitize(x, bins)

                # now we need to linearly interpolate between breakpoints
                corr = breakpoints[:, 1] + 1j*breakpoints[:, 2]
                x0 = np.where(binned > 0, bins[binned - 1], x)
                x1 = bins[binned]
                y0 = np.where(binned > 0, corr[binned - 1], 1.)
                y1 = corr[binned]
                y = y0 + (x - x0) * ((y1 - y0) / (x1 - x0))
                data[chan, ...] *= y / 10000.
            else:
                data[chan, ...] /= 10000.
        except:
            data[chan, ...] /= 10000.

    # DC correction
    if enable_dc_corr:
        dc_fixed_arr = basic_corr['dc_fixed_arr']
        dc_fixed_arr.shape = [dc_fixed_arr.shape[0], 1]
        data += dc_fixed_arr


    # measurement phase, PDA, REAR, and FEAR corrections
    random_phase = 2.*np.pi*float(lab['random_phase'][lab_pos]) / (2**16 - 1)
    measurement_phase = 0.5*np.pi*lab['measurement_phase'][lab_pos]
    gain_setting_index = lab['gain_setting_index'][lab_pos]
    if isMira:
        pda_corr = 1
    else:
        pda_corr = pda_factors[:, gain_setting_index]
        pda_corr.shape = [nchan, 1]
    echo_nr = int(lab['echo_nr'][lab_pos])
    measurement_sign = bool(lab['measurement_sign'][lab_pos])
    sign = bool(sin['spectrum_signs'][0][echo_nr]) != measurement_sign

    if not isMira:
        random_phase *= -1.0

    phase_factor = np.exp(1j * (random_phase - measurement_phase))

    corr_factor = pda_corr * phase_factor * coil_factor

    if relative_fear_bandwidth > 0:
        fear_corr = np.exp(1j * np.arange(data.shape[-1]) *
                           random_phase * relative_fear_bandwidth)
    else:
        fear_corr = 1.0

    if sign:
        data *= fear_corr * corr_factor
    else:
        data[:, :] = data[:, ::-1] * fear_corr * corr_factor

    return data


# read from the raw file
def readRaw(filename, raw_corr, chop_ky, cur_coil, cur_loc):
    maxRetroCardPhases = 100

    try:
        from .readMira import decodeMira
    except:
        print("Missing library for reading in R4 and R5 Philips "
              "lab/raw/sin data. Please contact your local Philips "
              "representative.")
        raise

    filename_lab = filename+'.lab'
    filename_raw = filename+'.raw'
    filename_sin = filename+'.sin'

    # read sin file
    if os.path.exists(filename_extcase(filename_sin)):
        sin = readSin(filename_extcase(filename_sin))
    else:
        print("Sin file not found!")

    # open the raw file
    if not type(filename_raw) is str:
        print("Input filename is not a string.")
        sys.exit(1)

    if os.path.splitext(filename_raw)[1] not in ['.raw', '.RAW']:
        print("input filename is not a .raw file")
        sys.exit(1)

    # opens the file
    try:
        fil_raw = open(filename_extcase(filename_raw), 'rb')
    except IOError:
        print('cannot open', filename_raw)
        sys.exit(1)

    # label file block size
    lab_pos = 0  # seek position
    isMira = sin['isMira']
    print("Mira Data:",isMira)

    # get number of channels
    try:
        nchan = int(sin['max_measured_channels'][0][0])
    except:
        try:
            nchan = int(sin['nr_measured_channels'][0][0])
        except:
            nchan = 1
    try:
        nchan_phc = int(sin['ph_nr_measured_channels'][0][0])
    except:
        nchan_phc = nchan

    # read the lab file
    lab = readLab(filename_lab, isMira)
    lab_length = len(lab['data_size'])

    # get the general correction terms from the sin file
    if raw_corr:
       basic_corr = setupCorrections(fil_raw, lab, sin, nchan, isMira)

    # find indices and seek offsets for data to be read
    std_case = lab['label_type'] == b'LABEL_TYPE_STANDARD'
    norm_case = lab['control'] == b'CTRL_NORMAL_DATA'
    echo_phase_case = lab['control'] == b'CTRL_ECHO_PHASE_DATA'
    noise_case = lab['control'] == b'CTRL_FRC_NOISE_DATA'
    norm_ind = (std_case & norm_case).nonzero()[0]
    phc_ind = (std_case & echo_phase_case).nonzero()[0]
    noi_ind = (std_case & noise_case).nonzero()[0]

    # generate ordered index sets
    mixes = oset_sorted(lab['mix_nr'][norm_ind])
    dynamics = oset_sorted(lab['dynamic_scan_nr'][norm_ind])
    phases = oset_sorted(lab['cardiac_phase_nr'][norm_ind])
    rtops = oset_sorted(lab['rtop_offset'][norm_ind])
    echoes = oset_sorted(lab['echo_nr'][norm_ind])
    e2vals = oset_sorted(lab['e2_profile_nr'][norm_ind])
    e1vals = oset_sorted(lab['e1_profile_nr'][norm_ind])
    e3vals = oset_sorted(lab['e3_profile_nr'][norm_ind])
    extras = oset_sorted(lab['extra_attr_nr'][norm_ind])
    rows = oset_sorted(lab['row_nr'][norm_ind])
    grad_echoes = oset_sorted(lab['grad_echo_nr'][norm_ind])
    rf_echoes = oset_sorted(lab['rf_echo_nr'][norm_ind])
    measurements = oset_sorted(lab['measurement_nr'][norm_ind])
    phc_measurements = oset_sorted(lab['measurement_nr'][phc_ind])
    # if a specific loc is chosen to be read
    loc_case = np.ones(lab_length) > 0  # create array of True values
    if cur_loc != -1:
        loc_case = lab['location_nr'] == cur_loc
    loc_ind = (loc_case & std_case & norm_case).nonzero()[0]
    locations = oset_sorted(lab['location_nr'][loc_ind])

    # determine at which indices data will be read
    if len(phc_ind) > 1:
        read_ind = (std_case & ((norm_case & loc_case) | noise_case |
                                echo_phase_case)).nonzero()[0]
    else:
        read_ind = (std_case & ((norm_case & loc_case) |
                    noise_case)).nonzero()[0]
    format_vals = lab['raw_format'] == 4
    format_vals |= lab['raw_format'] == 6
    if isMira:
        seek_array = np.cumsum(np.where(format_vals, lab['coded_data_size'],
                               lab['data_size']))
    else:
        seek_array = np.cumsum(lab['data_size'])
    seek_array = np.insert(seek_array, 0, 0) + 512
    seek_vals = seek_array[read_ind]

    # calculate number of samples to read
    if isMira:
        nsamp = int(max(lab['data_size'][norm_ind])) // (8 * nchan)
        if len(phc_ind) > 1:
            nsamp_phc = int(max(lab['data_size'][phc_ind])) // (8 * nchan_phc)
        nsamp_noi = int(max(lab['data_size'][noi_ind])) // (8 * nchan)
    else:
        nsamp = int(max(lab['data_size'][norm_ind])) // (4 * nchan)
        if len(phc_ind) > 1:
            nsamp_phc = int(max(lab['data_size'][phc_ind])) // (4 * nchan)
        nsamp_noi = int(max(lab['data_size'][noi_ind])) // (4 * nchan)

    # number of each index
    necho = len(echoes)
    nmeas = len(measurements)
    ndyn = len(dynamics)
    nloc = len(locations)
    nmix = len(mixes)
    ne2 = len(e2vals)
    ne1 = len(e1vals)
    ne3 = len(e3vals)
    nextra = len(extras)
    nrow = len(rows)
    ngrad = len(grad_echoes)
    nrfech = len(rf_echoes)
    nmeas_phc = len(phc_measurements)
    
    retroCardiac = False
    ncard = len(phases)
    if (ncard == 1) and (len(rtops) > 1):
        print("Retro cardiac mode active (beta)")
        retroCardiac = True
        ncard = maxRetroCardPhases # arbitrary, cut down later

    # allocate data arrays
    data_string = np.array(['chan', 'mix', 'dyn', 'card', 'echo', 'row',
                            'extra', 'loc', 'e3', 'meas', 'e2', 'e1', 'samp'])
    if cur_coil != -1:
        arr = np.zeros([1, nmix, ndyn, ncard, necho, nrow, nextra, nloc, ne3, 
                        nmeas, ne2, ne1, nsamp], dtype=np.complex64)
    else:
        arr = np.zeros([nchan, nmix, ndyn, ncard, necho, nrow, nextra, nloc, 
                        ne3, nmeas, ne2, ne1, nsamp], dtype=np.complex64)
    if len(phc_ind) > 1:
        phc_string = np.array(['chan', 'mix', 'card', 'echo', 'loc', 'e3',
                              'meas', 'e2', 'rf', 'grad', 'samp'])
        echo_arr = np.zeros([nchan_phc, nmix, ncard, necho, nloc, ne3,
                            nmeas_phc, ne2, nrfech, ngrad, nsamp_phc],
                            dtype=np.complex64)
    else:
        echo_arr = 0

    if retroCardiac:
        # 1D in the channel, ncard, samples dimensions
        cardTrackArr = np.zeros_like(arr[0,:,:,0,...,0], dtype = np.int16)
        rrArr = np.zeros_like(arr[0,...,0], dtype = np.int16)
        rtopArr = np.zeros_like(arr[0,...,0], dtype = np.int16)
        #Initialize rr/rtop to -1 to show where data were not loaded
        rrArr[:] = -1
        rtopArr[:] = -1
    else:
        rrArr = np.zeros([0], dtype = np.int16)
        rtopArr = np.zeros([0], dtype = np.int16)

    noi_string = np.array(['chan', 'samp'])
    noi_arr = np.zeros([nchan, nsamp_noi], dtype=np.complex64)

    # main data reading loop
    for seek_ind in range(len(seek_vals)):
        seek_offset = int(seek_vals[seek_ind])
        lab_pos = read_ind[seek_ind]

        # number of samples and channels depends on data type
        if lab['control'][lab_pos] == b'CTRL_FRC_NOISE_DATA':
            samples = nsamp_noi
            channels = nchan
            control_type = 'noise'
        elif lab['control'][lab_pos] == b'CTRL_ECHO_PHASE_DATA':
            samples = nsamp_phc
            channels = nchan_phc
            control_type = 'phc'
        elif lab['control'][lab_pos] == b'CTRL_NORMAL_DATA':
            samples = nsamp
            channels = nchan
            control_type = 'data'

        # if the data is from the digitized coil receivers then decode it
        act_data_size = int(lab['data_size'][lab_pos])
        if isMira:
            temp_data = np.zeros(act_data_size//4, dtype=np.int32)

            fil_raw.seek(seek_offset)
            bytebuff = fil_raw.read(int(
                lab['coded_data_size'][lab_pos]))

            try:  # c++ PyFI implementation of Mira decoding (fast)
                array_dims = np.array(temp_data.shape, np.int64)
                bb = np.frombuffer(bytebuff, np.uint8)
                temp_data = decodeMira(bb, array_dims, act_data_size)
            except:  # pure python implementation of Mira decoding (very slow)
                temp_data = decodeMira(bytebuff, temp_data, act_data_size)

        else:
            fil_raw.seek(seek_offset)
            bytebuff = fil_raw.read(act_data_size)
            temp_data = np.frombuffer(bytebuff, dtype=np.int16)

        try:
            if isMira:
                samp = int(lab['data_size'][lab_pos]) // (8 * nchan)
            else:
                samp = int(lab['data_size'][lab_pos]) // (4 * nchan)
            temp_data.shape = [channels, samp, 2]
            if samp < samples:
                pad = np.zeros([channels, samples-samp, 2])
                padded_data = np.append(temp_data, pad, 1)
            else:
                padded_data = temp_data
            cpxData = padded_data[:, :, 0] + 1j*padded_data[:, :, 1]
        except:
            print("Wrong data size. Skipping line...")
            continue

        if raw_corr:
            cxpData = basicCorrections(cpxData, lab, sin, lab_pos, basic_corr, isMira)

        # load cur indices from sets
        if control_type == 'phc':
            meas = phc_measurements.index(lab['measurement_nr'][lab_pos])
        else:
            meas = measurements.index(lab['measurement_nr'][lab_pos])
        if control_type != 'noise':
            echo = echoes.index(lab['echo_nr'][lab_pos])
            dyn = dynamics.index(lab['dynamic_scan_nr'][lab_pos])
            loc = locations.index(lab['location_nr'][lab_pos])
            mix = mixes.index(lab['mix_nr'][lab_pos])
            e2 = e2vals.index(lab['e2_profile_nr'][lab_pos])
            e1 = e1vals.index(lab['e1_profile_nr'][lab_pos])
            e3 = e3vals.index(lab['e3_profile_nr'][lab_pos])
            grad = grad_echoes.index(lab['grad_echo_nr'][lab_pos])
            rf = rf_echoes.index(lab['rf_echo_nr'][lab_pos])
            extra = extras.index(lab['extra_attr_nr'][lab_pos])
            row = rows.index(lab['row_nr'][lab_pos])
            if retroCardiac:
                # Check the label in the next position to determine if this profile should be rejected
                card_lab_pos = lab_pos + 1
                if lab['control'][card_lab_pos] == b'CTRL_INVALIDATION':
                    continue
                # Get current repetition number from the tracking array
                card = cardTrackArr[mix, dyn, echo, row, extra, loc,
                    e3, meas, e2, e1]
                if card >= maxRetroCardPhases:
                    print("Too many retrospective cardiac phases (max 250), aborting")
                    raise
                cardTrackArr[mix, dyn, echo, row, extra, loc,
                    e3, meas, e2, e1] += 1
                # Fill the rr and rtop arrays directly from labels
                rrArr[mix, dyn, card, echo, row, extra, loc,
                    e3, meas, e2, e1] = lab['rr_interval'][card_lab_pos]
                rtopArr[mix, dyn, card, echo, row, extra, loc,
                    e3, meas, e2, e1] = lab['rtop_offset'][card_lab_pos]
            else:
                card = phases.index(lab['cardiac_phase_nr'][lab_pos])

            # apply chop if appropriate
            if chop_ky and (e1 % 2) == 0:
                cpxData *= -1.0

        # accumulate data
        if control_type == 'noise':
            noi_arr[:, :] = cpxData
        elif control_type == 'phc':
            echo_arr[:, mix, card, echo, loc, e3, meas,
                     e2, rf, grad, :] = cpxData
        elif control_type == 'data':
            if cur_coil != -1:
                arr[0, mix, dyn, card, echo, row, extra, loc,
                    e3, meas, e2, e1, :] = cpxData[cur_coil, :]
            else:
                arr[:, mix, dyn, card, echo, row, extra, loc,
                    e3, meas, e2, e1, :] = cpxData

    print("finished reading data")
    data_labels = data_string[(np.array(
        arr.shape[0:len(data_string)]) > 1).nonzero()[0]]
    noi_labels = noi_string[(np.array(
        noi_arr.shape[0:len(data_string)]) > 1).nonzero()[0]]
    if len(phc_ind) > 1:
        phc_labels = phc_string[(np.array(
            echo_arr.shape[0:len(data_string)]) > 1).nonzero()[0]]
    else:
        phc_labels = 0

    # close the raw binary file
    fil_raw.close()

    hdr = dict()
    hdr['headerType'] = 'lab-sin'
    hdr['lab'] = lab
    hdr['sin'] = sin

    # Fix the phase that is caused by the decoding logic
    # For some reason the real and imaginary channels are swapped
    arr = np.conj(arr)
    arr *= np.exp(1j * np.pi/2.0)
    
    if retroCardiac:
        # Discard any unused phases
        maxCardNum = np.amax(cardTrackArr) - 1
        print("Maximum number of cardiac phases detected: ", maxCardNum)
        arr = arr[:,:,:,0:maxCardNum,...]
        rrArr = rrArr[:,:,0:maxCardNum,...]
        rtopArr = rtopArr[:,:,0:maxCardNum,...]

    # convert to numpy array
    return (arr, noi_arr, echo_arr, rrArr, rtopArr, hdr, data_labels, noi_labels, phc_labels)


# read a cpx file
def readCpx(filename, cur_coil, cur_loc):
    '''
        Parse and read a .cpx file.
    '''

    filename = filename_extcase(filename+'.cpx')

    # open the cpx file
    if not type(filename) is str:
        print("Input filename is not a string.")
        sys.exit(1)

    if os.path.splitext(filename)[1] not in ['.cpx', '.CPX']:
        print("input filename is not a .cpx file")
        sys.exit(1)

    # opens the file
    try:
        fil = open(filename, 'rb')
    except IOError:
        print('cannot open .cpx file ', filename)
        sys.exit(1)

    hdr = dict()

    index = 0
    offset = 0
    # pre-allocate info array
    info = np.zeros([10000, 127])
    fmt = '<'+'15i2f15i10f1q84i'
    # generate hdr table
    while True:
        fil.seek(offset)
        line = fil.read(512)

        # make sure line exists
        if len(line) == 512:
            info[index, :] = struct.unpack(fmt, line)
        else:
            break

        # calculate position of next header
        nx = info[index, 10]
        ny = info[index, 11]
        compression_factor = info[index, 13]

        if compression_factor == 0:
            break
        else:
            index += 1
            offset += 512 + int(nx*ny*8 // compression_factor)

        # put info into dictionary
        key = 'hdr_'+str(index)
        hdr[key] = info[index, :]
    hdr['headerType'] = 'cpx'

    # truncate info array
    info = info[:index, :]

    num_images = index

    # pre-allocate data array
    mixes = oset(info[:, 0])
    locs = oset(info[:, 1])
    echoes = oset(info[:, 4])
    phases = oset(info[:, 5])
    dynamics = oset(info[:, 6])
    rows = oset(info[:, 7])
    x_size = oset(info[:, 10])
    y_size = oset(info[:, 11])
    coils = oset(info[:, 21])

    nmix = len(mixes)
    nloc = len(locs)
    necho = len(echoes)
    ncard = len(phases)
    ndyn = len(dynamics)
    nrow = len(rows)
    nx = int(np.max(x_size))
    ny = int(np.max(y_size))
    nchan = len(coils)

    data_string = np.array(['chan', 'mix', 'dyn', 'card', 'echo', 'row',
                            'loc', 'y', 'x'])
    data = np.zeros([nchan, nmix, ndyn, ncard, necho, nrow, nloc, ny, nx],
                    dtype=np.complex64)

    # read in the cpx file
    offset = 512
    for index in range(num_images):

        fil.seek(offset)

        mix = mixes.index(info[index, 0])
        loc = locs.index(info[index, 1])
        echo = echoes.index(info[index, 4])
        card = phases.index(info[index, 5])
        dyn = dynamics.index(info[index, 6])
        row = rows.index(info[index, 7])
        coil = coils.index(info[index, 21])
        nx = int(info[index, 10])
        ny = int(info[index, 11])
        compression_factor = info[index, 13]
        size_bytes = int(8 * nx * ny // compression_factor)
        offset += 512 + size_bytes

        unparsed_data = fil.read(size_bytes)
        if compression_factor == 1:
            temp_data = np.frombuffer(unparsed_data, dtype=np.float32)
        elif compression_factor == 2:
            temp_data = np.frombuffer(unparsed_data, dtype=np.int16)
        elif compression_factor == 4:
            temp_data = np.frombuffer(unparsed_data, dtype=np.int8)
        temp_data.shape = [ny, nx, 2]
        complex_data = temp_data[:, :, 0] + 1j*temp_data[:, :, 1]
        data[coil, mix, dyn, card, echo, row, loc, :, :] = complex_data

    # setup the data labels
    data_labels = data_string[(np.array(
        data.shape[0:len(data_string)]) > 1).nonzero()[0]]

    return (data, hdr, data_labels)
