#! /usr/bin/env python
"""
Created on Mon Feb 12 16:45:57 2018

GE Headers
Find information about a GE PFile

@author: Josh Kaggie
"""
import numpy as np


DV24_sizes = np.array([157276,4096,16384,16384,98304,2052,2052,2048,1500,1960,2560,2448,7488])
DV26_sizes = np.array([163868,4096,16384,16384,98304,2052,2052,2048,2524,1960,2560,2448,7488,5568])
DV26_1_sizes = np.array([164532,4096,16384,16384,98304,2052,2052,2048,3188,1960,2560,2448,7488,5568])
DV26_2_sizes = np.array([213684,4096,16384,16384,147456,2052,2052,2048,3188,1960,2560,2448,7488,5568])
_header_sizes = {'26.000': ['26.000', 163868, 4096, 16384, 16384, 98304, 2052, 2052, 2048, 2524, 1960, 2560, 2448, 7488, 5568], '20.002': ['20.002', 145932, 2048, 16384, 16384, 98304, 2052, 2052, 2048, 1500, 1040, 2048, 2048], '20.003': ['20.003', 145932, 2072, 16384, 16384, 98304, 2052, 2052, 2048, 1500, 1040, 2048, 2048], '20.001': ['20.001', 145908, 2048, 16384, 16384, 98304, 2052, 2052, 2048, 1500, 1040, 2048, 2048], '20.006': ['20.006', 149788, 4096, 16384, 16384, 98304, 2052, 2052, 2048, 1500, 1960, 2560, 2448], '20.007': ['20.007', 149788, 4096, 16384, 16384, 98304, 2052, 2052, 2048, 1500, 1960, 2560, 2448], '20.004': ['20.004', 146564, 2704, 16384, 16384, 98304, 2052, 2052, 2048, 1500, 1040, 2048, 2048], '20.005': ['20.005', 149788, 4096, 16384, 16384, 98304, 2052, 2052, 2048, 1500, 1960, 2560, 2448], '24.000': ['24.000', 157276, 4096, 16384, 16384, 98304, 2052, 2052, 2048, 1500, 1960, 2560, 2448, 7488], '7': ['7', 39984, 2048, 4096, 4096, 20480, 2052, 2052, 2048, 0, 1040, 1028, 1044], '9': ['9', 61464, 2048, 4096, 4096, 40960, 2052, 2052, 2048, 0, 1040, 1536, 1536], '8': ['8', 60464, 2048, 4096, 4096, 40960, 2052, 2052, 2048, 0, 1040, 1028, 1044], '14.300000': ['14.300000', 145908, 2048, 16384, 16384, 98304, 2052, 2052, 2048, 1500, 1040, 2048, 2048], '16.000': ['16.000', 145908, 2048, 16384, 16384, 98304, 2052, 2052, 2048, 1500, 1040, 2048, 2048], '21.001': ['21.001', 150336, 4096, 16384, 16384, 98304, 2052, 2052, 2048, 2048, 1960, 2560, 2448], '11': ['11', 66072, 2048, 4096, 4096, 45056, 2052, 2052, 2048, 0, 1040, 2048, 1536], '10': ['10', 65560, 2048, 4096, 4096, 45056, 2052, 2052, 2048, 0, 1040, 1536, 1536], '15.000': ['15.000', 145908, 2048, 16384, 16384, 98304, 2052, 2052, 2048, 1500, 1040, 2048, 2048], '15.001': ['15.001', 145908, 2048, 16384, 16384, 98304, 2052, 2052, 2048, 1500, 1040, 2048, 2048], '14.1': ['14.1', 135704, 2048, 16384, 16384, 90112, 2052, 2052, 2048, 0, 1040, 2048, 1536], '14': ['14', 135704, 2048, 16384, 16384, 90112, 2052, 2052, 2048, 0, 1040, 2048, 1536], '14.2': ['14.2', 142356, 2048, 16384, 16384, 98304, 2052, 2052, 2048, 0, 1040, 2048, 2048], '25.003': ['25.003', 167196, 4096, 16384, 16384, 98304, 2052, 2052, 2048, 2524, 1960, 2560, 2448, 7488, 8896], '25.002': ['25.002', 166172, 4096, 16384, 16384, 98304, 2052, 2052, 2048, 1500, 1960, 2560, 2448, 7488, 8896], '25.001': ['25.001', 166172, 4096, 16384, 16384, 98304, 2052, 2052, 2048, 1500, 1960, 2560, 2448, 7488, 8896],
                 '26.001': DV26_1_sizes, '26.002': DV26_2_sizes}



header_dict = {
    24.000: DV24_sizes,
    26.000: DV26_sizes,
    26.001: DV26_1_sizes,
    26.002: DV26_2_sizes,
    }

#f['image']['mr_flip'], f['image']['tr'], f['image']['te']
#f['rdb']['te'], f['rdb']['fov']
#f['rdb']['ps_mps_tg']

## 
def get_multiple_file_info(files = [], infotype='rdb',keyin = 'ps_mps_tg',path='/media/jk636/d/Projects/MRF Repeat Brain/data/*/*/'):
    if files == []:
        import glob
        files = glob.glob(path+'P*7')
    for filename in files:
        print(filename + str(get_rdb_info(filename=filename)[keyin]))



def get_user_cvs_path(path='/media/jk636/d/Projects/MRF Renal/new_fidall/',pattern='P*7'):
    import glob
    path = fix_path(path)    
    files = glob.glob(path+pattern)    
    for filename in files:
        print(filename[-8:],get_user_cvs(filename))


#Read Database Header 
def get_rdb_info(filename = '/media/jk636/d/Projects/MRF Repeat Brain/data/3TCam/martin/P66560.7'):  
    version = get_rdb_version(filename)        
    fp = open(filename,'rb')        
    fp.seek(0)
    if version >26.0:
        rdbhdr = np.fromfile(fp, DV26_2_rdbhdr,1)
    elif version == 24.0:
        rdbhdr = np.fromfile(fp, DV24_rdbhdr,1)        
    fp.close()
    return rdbhdr

##Image Header
def get_image_info(filename = '/media/jk636/d/Projects/MRF Repeat Brain/data/3TCam/martin/P66560.7'):  
    version = get_rdb_version(filename)        
    fp = open(filename,'rb')        
    fp.seek(0)
    if version >26.0:
        rdbhdr = np.fromfile(fp, DV26_2_rdbhdr,1)
        off_image = rdbhdr['off_image']        
        fp.seek(int(off_image))
        image_info = np.fromfile(fp, DV26_2_image_header,1)
    fp.close()
    return image_info


def get_MRF_string(filename = '/media/jk636/d/Projects/MRF Repeat Brain/data/3TCam/martin/P66560.7'):  
    cvs = get_user_cvs(filename)
    infos=get_rdb_exam_image(filename)
    if cvs[20] == 0:
        tr = infos['image']['tr']
    else:
        tr = cvs[20]    
    tot_str = 'TE'+str(cvs[16])+'TR'+str(tr)+'V21'+str(cvs[20])+'V24'+str(cvs[24])
    tot_str = tot_str.replace('.','_').replace('[','').replace(']','')
    return tot_str

##Image Header
def get_exam_info(filename = '/media/jk636/d/Projects/MRF Repeat Brain/data/3TCam/martin/P66560.7'):  
    version = get_rdb_version(filename)        
    fp = open(filename,'rb')        
    fp.seek(0)
    if version >26.0:
        rdbhdr = np.fromfile(fp, DV26_2_rdbhdr,1)
        offset = rdbhdr['off_exam']        
        fp.seek(int(offset))
        image_info = np.fromfile(fp, DV26_2_exam_header,1)
    fp.close()
    return image_info


###Raw data
def get_data(filename = '/media/jk636/d/Projects/MRF Repeat Brain/data/3TCam/martin/P66560.7',asraw=False):  
    rdb_info = get_rdb_info(filename)    
    off_data = rdb_info['off_data']
    version = get_rdb_version(filename)        
    fp = open(filename,'rb')        
    fp.seek(off_data)    
    data = np.fromfile(fp,np.int32)    
    coils = rdb_info['dab(id)'][0][1]-rdb_info['dab(id)'][0][0]+1
    xres = rdb_info['da_xres'][0]
    yres = rdb_info['da_yres'][0]
    nslices = rdb_info['nslices'][0]    
    if asraw:
        return data
    #data = data.reshape([2,xres,coils,yres,nslices])
    data = data[::2] + 1j* data[1::2]
    data = data.reshape([nslices, yres, coils, xres])   
    #header.get_exam_info()['magstrength']
    return data



def get_akwav_number(filename = '/media/jk636/d/Projects/MRF Repeat Brain/data/3TCam/martin/P66560.7'):    
    version = get_rdb_version(filename)        
    fp = open(filename,'rb')        
    fp.seek(0)
    if version >26.0:
        rdbhdr = np.fromfile(fp, DV26_2_rdbhdr,1)
        off_image = rdbhdr['off_image']        
        fp.seek(int(off_image))
        image_info = np.fromfile(fp, DV26_2_image_header,1)
    fp.close()
    return image_info['user3'][0]


def get_user_cvs(filename = '/media/jk636/d/Projects/MRF Repeat Brain/data/3TCam/martin/P66560.7'):    
    version = get_rdb_version(filename)        
    fp = open(filename,'rb')        
    fp.seek(0)
    if version >26.0:
        rdbhdr = np.fromfile(fp, DV26_2_rdbhdr,1)
        off_image = rdbhdr['off_image']        
        fp.seek(int(off_image))
        image_info = np.fromfile(fp, DV26_2_image_header,1)
    fp.close()
    out_array = []
    for xi in xrange(25):
        out_array.append(image_info['user'+str(xi)][0])
    return out_array


def get_rdb_exam_image(filename = '/media/jk636/d/Projects/MRF Repeat Brain/data/3TCam/martin/P66560.7'):  
    rdb_info = get_rdb_info(filename)    
    exam_info = get_exam_info(filename)
    image_info = get_image_info(filename)
    return {'rdb':rdb_info, 'exam':exam_info, 'image':image_info}


def get_basics(filename = '/media/jk636/d/Projects/MRF Repeat Brain/data/3TCam/martin/P66560.7'):  
    rdb_info = get_rdb_info(filename)    
    exam_info = get_exam_info(filename)
    image_info = get_image_info(filename)
    
    #coils = rdb_info['dab(1)']-rdb_info['dab(0)']+1
    coils = rdb_info['dab(id)'][0][1]-rdb_info['dab(id)'][0][0]+1
    xres = rdb_info['da_xres'][0]
    yres = rdb_info['da_yres'][0]
    nslices = rdb_info['nslices'][0]
    tg = rdb_info['ps_aps_tg'][0], rdb_info['ps_mps_tg'][0]
    freq = rdb_info['ps_aps_freq'][0], rdb_info['ps_mps_freq'][0]
    broadband_select = rdb_info['broad_band_select'][0]                
    scan_time = rdb_info['scan_time'].tostring().replace('\x00','')
    scan_date = rdb_info['scan_date'].tostring().replace('\x00','')
    magstrength = exam_info['magstrength'][0]
    extras = tostring(exam_info['hospname']), tostring(exam_info['ex_desc'])  ##needs work
    psdname = tostring(image_info['psdname'])
    return psdname, scan_time, scan_date, magstrength, coils, xres, yres, nslices, tg, freq, broadband_select, extras

    #####scan_time = rdb_info[.exam['ex_desc'].tostring().replace('\x00',''),  WE NEED THIS ONE!
    #if 'dab(1)' in rdb_info.keys():
    #    print('')    
    #(1, array([1344], dtype=uint16), array([979], dtype=int16), array([22], dtype=uint16))
    #as['basics']['n_coils']*datas['basics']['n_coils']*datas['basics']['n_coils']*2
    


def get_psdname(filename = '/media/jk636/d/Projects/MRF Repeat Brain/data/3TCam/martin/P66560.7'):    
    version = get_rdb_version(filename)    
    fp = open(filename,'rb')        
    fp.seek(0)
    if version >26.0:
        rdbhdr = np.fromfile(fp, DV26_2_rdbhdr,1)
        off_image = rdbhdr['off_image']        
        fp.seek(int(off_image))
        image_info = np.fromfile(fp, DV26_2_image_header,1)
    fp.seek(int(off_image)  )
    return tostring(image_info['psdname'])


    
def tostring(number_array):
    return number_array.tostring().replace('\x00','')
    

def get_rdb_version(filename = '/media/jk636/d/Projects/MRF Repeat Brain/data/3TCam/martin/P66560.7'):    
    fp = open(filename,'rb')        
    fp.seek(0)
    version = np.fromfile(fp,np.float32,1)
    fp.close()
    return version[0]
    
def get_header_types(rdb_version=26.002):
    if rdb_version == 26.002:
        rdb_hdr = DV24_rdbhdr
        DV24_image_header 
        DV24_exam_header 
        DV24_series_header


def list_directory(path='/usr/g/mrraw',pattern='P*7',printlist=True,returnlist=False, new_method = True,argv = None):
    if len(argv)==2:
        print(argv)
        path = argv[1]    
    import glob
    path = fix_path(path)    
    files = glob.glob(path+pattern)
    
    for filename in files:
        version = get_rdb_version(filename)    
        print(filename, version)
        if version > 26.0:
            print(get_basics(filename))
            print
        else:
            print('Not >DV26.0')            
        
    #files.sort(key=os.path.getmtime)
    returnvalues = {}
    print files
    
def fix_path(path=''):
    if path[-1] != '/':
        path = path+'/'
    return path




def convert_hdr_sizes(filename='/home/jk636/fidalltest/hdr_sizes.txt'):
    f = open(filename,'r')
    collection = {}
    reading = True
    collect_version = []
    while reading:
        reading = False
        for line in f.readlines():
            #print line, reading
            line = line.strip()
            if '#end'  in line:
                print collect_version[0]
                collection[collect_version[0]] = collect_version
            elif '#' in line:
                collect_version = [line.strip().replace('#','')]
                print line,
            elif line == '#24':
                reading = False
            else:
                if len(line)>0 and line[0] != '#':
                    collect_version.append(int(line.strip()))
    f.close()
    return collection
        

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
    new_endstr = ''
    in_for_loop = False
    try:
        for line in endstr.splits('\n'):
            if 'for id' in line:
                in_for_loop = True
            if 'end' in line:
                in_for_loop = False
            if in_for_loop:
                new_endstr = ''
            else:
                new_endstr = new_endstr+line+'\n'
        #('int_padding1(id)', np.int32,31),
        return new_endstr
    except:
        return endstr









DV24_rdbhdr = np.dtype( [
  ('rdbm_rev',   np.float32),
  ('run_int',   '<i4'),
  ('scan_seq',   '<i2'),
  ('run_char', '<i1'  , 6), 
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
  ])

DV24_image_header = np.dtype([
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
    ])


DV24_exam_header = np.dtype([
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
  ('short_padding(id)', np.int16,35)
  ])       
       
       
       
      
DV24_series_header = np.dtype([
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
])      
   

DV24_data_acq_tab = np.dtype([
  ('pass_number', np.int16),
  ('slice_in_pass', np.int16),
  ('gw_point1(id)', np.float32,3),
  ('gw_point2(id)', np.float32,3),
  ('gw_point3(id)', np.float32,3),
  ('transpose', np.int16),
  ('rotate', np.int16),
  ('coilConfigUID', np.uint32),
])













DV26_2_exam_header= np.dtype( [
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
  ('short_padding(id)', np.int16,35),
  ('hist', np.int32, 257), 
  ('refphy', np.int32, 65), 
  ('diagrad', np.int32, 65),
  ('operator_new', np.int32, 65), 
  ('ex_desc', np.int32, 65),
  ('ex_typ', np.int32, 3), 
  ('ex_sysid', np.int32, 17), 
  ('ex_alloc_key', np.int32, 13), 
  ('ex_diskid', '<i1'),
  ('hospname', np.int32, 33)
])








DV26_2_rdbhdr = np.dtype( [
  ('rdbm_rev',   np.float32),
  ('off_data',   np.int32),
  ('off_per_pass',   np.int32),
  ('off_unlock_raw',   np.int32),
  ('off_data_acq_tab',   np.int32),
  ('off_nex_tab',   np.int32),
  ('off_nex_abort_tab',   np.int32),
  ('off_tool',   np.int32),
  ('off_exam',   np.int32),
  ('off_series',   np.int32),
  ('off_image',   np.int32),
  ('off_ps',   np.int32),
  ('off_grad_data',   np.int32),
  ('off_CTT_data',   np.int32),
  ('off_spare1',   np.int32),
  ('off_spare2',   np.int32),
  ('off_spare3',   np.int32),
  ('off_spare4',   np.int32),
  ('off_spare5',   np.int32),
  ('off_spare6',   np.int32),
  ('run_int',   np.int32),
  ('scan_seq',   np.int16),
  ('run_char', '<i1'  , 6), 
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
  ])



DV26_2_rdbhdr = [\
  ('rdbm_rev', np.float32),
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
  ('off_grad_data', np.int32),
  ('off_CTT_data', np.int32),
  ('off_spare1', np.int32),
  ('off_spare2', np.int32),
  ('off_spare3', np.int32),
  ('off_spare4', np.int32),
  ('off_spare5', np.int32),
  ('off_spare6', np.int32),
  ('run_int', np.int32),
  ('scan_seq', np.int16),
  ('run_char', '<i1'  , 6), 
  ('scan_date', '<i1'  , 10), 
  ('scan_time', '<i1'  , 8), 
  ('logo', '<i1'  , 10), 
  ('file_contents', np.int16),
  ('lock_mode', np.int16),
  ('dacq_ctrl', np.int16),
  ('recon_ctrl', np.uint16),
  ('exec_ctrl', np.uint16),
  ('scan_type', np.int16),
  ('data_collect_type', np.int16),
  ('data_format', np.int16),
  ('recon', np.int16),
  ('datacq', np.int16),
  ('npasses', np.int16),
  ('npomp', np.int16),
  ('nslices', np.uint16),
  ('nechoes', np.int16),
  ('navs', np.int16),
  ('nframes', np.int16),
  ('baseline_views', np.int16),
  ('hnover', np.int16),
  ('frame_size', np.uint16),
  ('point_size', np.int16),
  ('vquant', np.int16),
  ('cheart', np.int16),
  ('ctr', np.float32),
  ('ctrr', np.float32),
  ('initpass', np.int16),
  ('incrpass', np.int16),
  ('method_ctrl', np.int16),
  ('da_xres', np.uint16),
  ('da_yres', np.int16),
  ('rc_xres', np.int16),
  ('rc_yres', np.int16),
  ('im_size', np.int16),
  ('rc_zres', np.int32),
  ('fermi_radius', np.float32),
  ('fermi_width', np.float32),
  ('fermi_ecc', np.float32),
  ('clip_min', np.float32),
  ('clip_max', np.float32),
  ('default_offset', np.float32),
  ('xoff', np.float32),
  ('yoff', np.float32),
  ('nwin', np.float32),
  ('ntran', np.float32),
  ('scalei', np.float32),
  ('scaleq', np.float32),
  ('rotation', np.int16),
  ('transpose', np.int16),
  ('kissoff_views', np.int16),
  ('slblank', np.int16),
  ('gradcoil', np.int16),
  ('ddaover', np.int16),
  ('sarr', np.int16),
  ('fd_tr', np.int16),
  ('fd_te', np.int16),
  ('fd_ctrl', np.int16),
  ('algor_num', np.int16),
  ('fd_df_dec', np.int16),
  ('dab(id)', np.int16,8),
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
  ('v_type', np.int32),
  ('v_coefxa', np.float32),
  ('v_coefxb', np.float32),
  ('v_coefxc', np.float32),
  ('v_coefxd', np.float32),
  ('v_coefya', np.float32),
  ('v_coefyb', np.float32),
  ('v_coefyc', np.float32),
  ('v_coefyd', np.float32),
  ('v_coefza', np.float32),
  ('v_coefzb', np.float32),
  ('v_coefzc', np.float32),
  ('v_coefzd', np.float32),
  ('vm_coef1', np.float32),
  ('vm_coef2', np.float32),
  ('vm_coef3', np.float32),
  ('vm_coef4', np.float32),
  ('v_venc', np.float32),
  ('spectral_width', np.float32),
  ('csi_dims', np.int16),
  ('xcsi', np.int16),
  ('ycsi', np.int16),
  ('zcsi', np.int16),
  ('roilenx', np.float32),
  ('roileny', np.float32),
  ('roilenz', np.float32),
  ('roilocx', np.float32),
  ('roilocy', np.float32),
  ('roilocz', np.float32),
  ('numdwell', np.float32),
  ('ps_command', np.int32),
  ('ps_mps_r1', np.int32),
  ('ps_mps_r2', np.int32),
  ('ps_mps_tg', np.int32),
  ('ps_mps_freq', np.uint32),
  ('ps_aps_r1', np.int32),
  ('ps_aps_r2', np.int32),
  ('ps_aps_tg', np.int32),
  ('ps_aps_freq', np.uint32),
  ('ps_scalei', np.float32),
  ('ps_scaleq', np.float32),
  ('ps_snr_warning', np.int32),
  ('ps_aps_or_mps', np.int32),
  ('ps_mps_bitmap', np.int32),
  ('ps_powerspec', '<i1', 256), #// readdc
  ('ps_filler1', np.int32),
  ('ps_filler2', np.int32),
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
  ('feextra', np.int32),
  ('raw_pass_size', np.uint64),
  ('sspsave', np.uint64),
  ('udasave', np.uint64),
  ('vibrant', np.int16),
  ('asset_torso', np.int16),
  ('asset_alt_cal', np.int32),
  ('kacq_uid', np.int32),
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
  ('hoecc', np.int32),
  ('hoec_fit_order', np.int32),
  ('esp', np.int32),
  ('viewSharing3D', np.uint32),
  ('gradwarp_interp_type', np.int32),
  ('pure_blur_enable', np.int32),
  ('pure_blur', np.float32),
  ('slice_info_tab_checksum', np.uint32),
  ('mb_factor', np.uint32),
  ('mb_slice_fov_shift_factor', np.uint32),
  ('mb_slice_order_sign', np.int32),
  ('cal_options', np.uint32),
  ('cs_factor', np.float32),
  ('cs_flag', np.int32),
  ('cs_maxiter', np.int32),
  ('cs_consistency', np.float32),
  ('cs_ph_stride', np.int32),
  ('cs_sl_stride', np.int32),
  ('noncart_traj_mode', np.int32),
  ('pure_alpha', np.float32),
  ('pure_otsu_class', np.int32),
  ('pure_exp_wt', np.float32),
  ('pure_erode_dist', np.int32),
  ('pure_dilate_dist', np.int32),
  ('pure_aniso_blur', np.int32),
  ('pure_aniso_erode_dist', np.int32),
  ('pure_aniso_dilate_dist', np.int32),
  ('calmode', np.int32),
  ('cine_ctrl', np.uint16),
  ('cine_ncycles', np.uint16),
  ('cine_recon_phases', np.uint16),
  ('cine_acq_phases', np.uint16),
  ('cine_virtual_coils', np.uint16),
  ('cine_nrejects', np.uint16),
  ('cine_low_rr', np.uint16),
  ('cine_high_rr', np.uint16),
  ('cine_avg_rr', np.uint16),
  ('cine_bpm', np.uint16),
  ('kat_arc_cal_size_ky', np.uint16),
  ('kat_arc_cal_size_kz', np.uint16),
  ('kat_arc_accel_outer', np.uint16),
  ('kat_arc_accel_center', np.uint16),
  ('cardt1map_ctrl', np.uint16),
  ('cardt1map_ti', np.uint16),
  ('propellerCtrl', np.int32),
  ('daviewsPerBlade', np.int32),
  ('rotationThreshold', np.float32),
  ('shiftThreshold', np.float32),
  ('correlationThreshold', np.float32),
  ('npwfactor', np.float32),
  ('viewshd_num', np.int32),
  ('bvalstab(id)', np.float32,100),
  ('synbvalstab(id)', np.float32,10),
  ('numsynbvals', np.uint16),
  ('position_detection', np.uint16),
  ('medal_param', np.float32),
  ('clariview_filter_x_nr_val', np.uint16),
  ('clariview_filter_x_sh_val', np.uint16),
  ('clariview_filter_x_is_factory', np.uint16),
  ('clariview_filter_x_option', np.uint16),
  ('grad_intensity_thres', np.float32),
  ('anne_mask_dilation_length', np.int32),
  ('anne_intensity_thres', np.float32),
  ('anne_channel_percentage', np.float32),
  ('ec3_iterations', np.int32),
  ('combined_coil_map_thres', np.float32),
  ('coil_map_smooth_factor', np.float32),
  ('coil_map_2_filter_width', np.float32),
  ('ec3_iteration_method', np.int32),
  ('img_ctrl', np.int32),
  ('moco_ctrl', np.int32),
  ('excess', np.int16,526),
  ]
  


DV26_2_image_header   = [\
  ('autoSubParam.seriesUidToSubtract', '<i1',32),
  ('autoSubParam.imageNoToSubtract', np.int32),
  ('autoSubParam.destSeriesNo', np.int32),
  ('autoSubParam.destImageNo', np.int32),
  ('autoSubParam.dummy', np.int32),
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
  ('sarValues(id)', np.float32,4),
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
  ('CSAccel', np.float32),
  ('esp', np.float32),
  ('bwPerPix', np.float32),
  ('pixelSizeX', np.float32),
  ('pixelSizeY', np.float32),
  ('fov_freq_scale', np.float32),
  ('fov_phase_scale', np.float32),
  ('slthick_scale', np.float32),
  ('freq_loc_shift', np.float32),
  ('phase_loc_shift', np.float32),
  ('slice_loc_shift', np.float32),
  ('b1PeakLowSAR', np.float32),
  ('b1rmsLowSAR', np.float32),
  ('wbLowSAR', np.float32),
  ('headLowSAR', np.float32),
  ('seriesMaxTimeLowSAR', np.float32),
  ('float_padding(id)', np.float32,9),
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
  ('rtb0cc', np.int32),
  ('user_bitmap2', np.int32),
  ('MultiBandAccel', np.int32),
  ('tissueT1', np.int32),
  ('syn_acq_b_info', np.int32),
  ('Acquirephase', np.int32),
  ('HKTAccel', np.int32),
  ('saveOriginalImage', np.int32),
  ('cardt1map_hb_pattern', np.int32),
  ('int_padding2(id)', np.int32,12),
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
  ('fat_shift_dir', np.int16),
  ('ab1_status', np.int16),
  ('ab1_tg', np.int16),
  ('vest', np.int16),
  ('msde_dir', np.int16),
  ('put_in_database', np.int16),
  ('scout_scan', np.int16),
  ('short_padding(id)', np.int16,28),
  ('psdname', '<i1',33),
  ('proj_name', '<i1',13),
  ('psd_iname', '<i1',13),
  ('im_diskid', '<i1'), 
  ('pdid(id)', '<i1',14),
  ('im_suid(id)', '<i1',4),
  ('contrastIV', '<i1',17),
  ('contrastOral', '<i1',17),
  ('loc_ras', '<i1'),
  ('forimgrev', '<i1',4),
  ('cname', '<i1',17),]


  
'''  ('im_verscre(id)', np.char,2),]'''
'''('im_verscur(id)', np.char,2)]'''
  
'''('im_alloc_key', np.int32, 13c_key', np.int32,REPLACEEND 13)),
  ('ref_img', np.char),
  ('sum_img', np.char),
  ('filter_mode', np.int32, 16mode', np.int32,REPLACEEND 16)),
  ('slop_str_2', np.int32, 16r_2', np.int32,REPLACEEND 16)), 
  ('image_uid(id)', np.char,32),
  ('sop_uid(id)', np.char,32), 
  ('GEcname(id)', np.char,24),
  ('usedCoilData(id)', np.char,100),
  ('astcalseriesuid(id)', np.char,32),
  end
  for id = 1 : 32
('purecalseriesuid(id)', np.char),
  end
  for id = 1 : 32
('xml_psc_shm_vol(id)', np.char),
  end
  for id = 1 : 64
('rxmpath(id)', np.char),
  end
  ('psdnameannot', np.int32, 33annot', np.int32,REPLACEEND 33)),
  for id = 1 : 55
('anatomy(id)', np.char),
  end
  for id = 1 : 195
('img_hdr_padding(id)', np.char),
  end'''


if __name__ == "__main__":
    import sys
    list_directory(argv = sys.argv)