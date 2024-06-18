function hdr=write_dicom(fname,b,h,inp,dcmdict,move2dir,verb)
%WRITE_DICOM  Write GE MR compatible dicom's for importing into database
% hdr = write_dicom(fname,b,h,inp,dcmdict,move2dir,verb)
%    fname  Filename of dicom image to create
%        b  2D Images to store (max 5D)
%           multislice: slices along 5th dimension
%           3D: slices along 3rd dimension
%        h  P-file header structure
%      inp  Input structure to copy to dicom header
%  dcmdict  GE MR Dicom dictionary (optional)
%           defaults=getenv('DCMDICT') and 'gems-dicom-dict.txt'
% move2dir  Move saved files to import directory (0=off)             (1)
%           if getenv('IMPORT_IMAGE_DIR') exists: 1=move; 2=copy
%     verb  Verbose mode: 0=off, 1=normal, 2=extensive               (1)
%      hdr  return header info from template + modification
%
% 12/2020 Rolf Schulte
if (nargin<1), help(mfilename); return; end
if (nargin<3), warning('At least three input arguments required'); end

max_scale = 2^15-2;


%% check input input
if ~exist('b','var'),        b = []; end
if isempty(b), warning('b=[] -> exiting'); return; end
if ~exist('inp','var'),      inp = []; end
if ~exist('dcmdict','var'),  dcmdict = []; end
if ~exist('move2dir','var'), move2dir = []; end
if isempty(move2dir),        move2dir = 1; end
if ~exist('verb','var'),     verb = []; end
if isempty(verb),            verb = 1; end

if length(size(b))>5
    warning('Dimensions(=%g)>5D',length(size(b)));
    fprintf('Exiting write_dicom without storing\n');
    return;
end
if ~isreal(b), warning('b complex: taking magnitude'); b = abs(b); end


%% look for dicom dictionary
if isempty(dcmdict)
    dcmdict = getenv('DCMDICT');
    if ~exist(dcmdict,'file')
        % dcmdict = 'gems-dicom-dict.txt';
        % avoiding bug in matlab not finding dictionary in search path
        dcmdict = which('gems-dicom-dict.txt');
    end
end
if ~exist(dcmdict,'file')
    warning('Dicom dictionary %s not found',dcmdict);
    fprintf('Exiting write_dicom without storing\n');
    return;
end


%% load dicom header from template
if verb>0, fprintf('write_dicom ****** \n'); end
nn = size(b);
nn = [nn,ones(1,5-length(nn))];
hdr = sub_convert_header(inp,h,nn,verb);
if ~isempty(regexpi(hdr.InternalPulseSeqName,'fidall'))
    is_fidall = true;
else
    is_fidall = false;
end
if ~isempty(regexpi(hdr.InternalPulseSeqName,'3drad')) || ...
        ~isempty(regexpi(hdr.InternalPulseSeqName,'burzte'))
    is_zte = true;
else
    is_zte = false;
end


%% creating directory
if ~isempty(regexpi(fname,'\.dcm')), fname = fname(1:end-4); end
fname = [fname '_dcm'];
if verb>0, fprintf('Storing dicom data in\n\t%s\n',fname); end
if exist(fname,'dir')
    warning('directory ''%s'' exists; removing',fname);
    rmdir(fname,'s');
end
mkdir(fname);


%% adjust scaling
maxval = max(abs(b(:)));
b = int16(b/maxval*max_scale);     % actual rescaling to +- 2^15-2
hdr.UserData16 = maxval;
hdr.WindowCenter = floor(max_scale/2);
hdr.WindowWidth = floor(max_scale);
hdr.BitDepth = 16;


%% determine acquisition type
SpacingBetweenSlices = h.image.slthick + h.image.scanspacing;
threeD = false;
if ~isempty(strfind(hdr.MRAcquisitionType,'3'))
    SpacingBetweenSlices = h.rdb_hdr.fov/nn(3);
    threeD = true;
end
if ~threeD && (nn(5)>1)
    multiSlice = true;
else
    multiSlice = false;
end


%% adapting slice location
start_loc = h.series.start_loc;   % slice location first slice
end_loc = h.series.end_loc;       % slice location last slice
if abs(h.image.loc-start_loc)>1d-5
    warning('h.image.loc(=%g)~=h.series.start_loc(=%g)',h.image.loc,start_loc); 
end
ipp = sub_get_ipp(h,h.rdb_hdr.fov./nn(1:2)); % ImagePositionPatient first slice
% if  threeD || multiSlice
hdr.SpacingBetweenSlices = SpacingBetweenSlices;
vec = [-h.image.norm_R ; -h.image.norm_A ; h.image.norm_S]; % dir of next slices
if abs(sqrt(sum(vec.^2))-1)>1d-3
    warning('length of vec(=%g)~=1',sqrt(sum(vec.^2)));
end
vec = vec*SpacingBetweenSlices;    % scale to [mm]
ori = 'o';          % oblique is default orientation
av = abs(vec);
dloc = max(av);
if (av(1)>1000*av(2))&&(av(1)>1000*av(3)), ori = 's'; end
if (av(2)>1000*av(1))&&(av(2)>1000*av(3)), ori = 'c'; end
if (av(3)>1000*av(1))&&(av(3)>1000*av(2)), ori = 'a'; end
if threeD && is_fidall
    start_loc = start_loc - dloc*nn(3)/2;
    end_loc   = end_loc   + dloc*nn(3)/2;
    ipp = ipp - vec*nn(3)/2;
end
if is_zte && h.image.plane==4 % sagital
    % empirical correction
    % NOTE: ZTE product recon reverses sagital+coronal, hence z
    % inverted in recon_zte
    vec = -vec;
end

if start_loc > end_loc-1d-10
    if verb>1, fprintf('Inverting dloc and vec\n'); end
    dloc = -dloc;
    vec = -vec;
end
if verb>0
    fprintf('ImagePositionPatient = [%g,%g,%g]\n',ipp(1),ipp(2),ipp(3));
    fprintf('ImageOrientationPatient = [%g %g %g %g %g %g]\n',...
        hdr.ImageOrientationPatient);
    fprintf('start_loc = %g; end_loc = %g\n',start_loc,end_loc);
    fprintf('slice vector = [%g,%g,%g]\n',vec(1),vec(2),vec(3));
    fprintf('orientation = %s; dloc = %g\n',ori,dloc);
end

% else
if ~(threeD || multiSlice)
    hdr.ImagePositionPatient = ipp;
    hdr.SliceLocation = start_loc;
    hdr.ImageLocation = start_loc;
    hdr.RASLetterOfImageLocation = sub_get_ras(ori,start_loc);
end

% locations of first and last slices
hdr.FirstScanLocation = abs(round(start_loc*1d4)*1d-4);
hdr.FirstScanRAS = sub_get_ras(ori,start_loc);
hdr.LastScanLocation = abs(round(end_loc*1d4)*1d-4);
hdr.LastScanRAS = sub_get_ras(ori,end_loc);



%% saving
if verb>0, fprintf('Starting loop for storing\n'); end
ll = 0;
for l5=1:nn(5)
    for l4=1:nn(4)
        for l3=1:nn(3)
            if verb>1
                fprintf('dicomwrite %g/%g (l3=%g,l4=%g,l5=%g)\t',...
                    ll+1,nn(3)*nn(4)*nn(5),l3,l4,l5);
            else
                if (verb==1) && (nn(3)*nn(4)*nn(5)>20) && ...
                        (mod(ll,nn(3)*nn(4)*nn(5)/20)<1)
                    fprintf('.');
                end
            end
            if threeD, li = l3; end
            if multiSlice, li = l5; end
            if threeD || multiSlice
                hdr.ImagePositionPatient = ipp+vec*(li-1);
                loc = round((start_loc+dloc*(li-1))*1d4)*1d-4;

                hdr.SliceLocation = loc;
                hdr.ImageLocation = loc;
                hdr.RASLetterOfImageLocation = sub_get_ras(ori,loc);
                hdr.InStackPositionNumber = li;
                if verb>1
                    fprintf('\tipp=[%g,%g,%g]\tloc=%g RAS=%s\n',...
                        hdr.ImagePositionPatient,loc,...
                        hdr.RASLetterOfImageLocation);
                end
                if (li==1) && abs(loc-start_loc)>1d-3
                    warning('loc(=%g)~=start_loc(=%g)',loc,start_loc);
                end
                if ((multiSlice && (li==nn(5))) || (threeD && (li==nn(3)))) ...
                        && (abs(loc-end_loc)>1d-3)
                    warning('loc(=%g)~=end_loc(=%g)',loc,end_loc);
                end
            else    
                if verb>1, fprintf('\n'); end
            end
            ll = ll+1;
            fn = sprintf('%s/image%.4d.dcm',fname,ll);
            hdr.InstanceNumber = ll;

            dicomwrite(b(:,:,l3,l4,l5),fn,hdr,...
                'dictionary',dcmdict,'Endian','ieee-le',...
                'VR','explicit','ObjectType','MR Image Storage',...
                'WritePrivate',true);
            
            if (ll==1), warning('off'); end
        end
    end
end
if verb==1, fprintf('\n'); end
warning('on');


%% move to import directory
IMPORT_IMAGE_DIR = getenv('IMPORT_IMAGE_DIR');
if (move2dir>0) && ~isempty(IMPORT_IMAGE_DIR)
    if verb>0
        if move2dir==1
            fprintf('Moving files for import into GEMR database to\n');
        else
            fprintf('Copying files for import into GEMR database to\n');
        end
    end
    fprintf('\t%s\n',IMPORT_IMAGE_DIR);
    fn_dcmexp = [IMPORT_IMAGE_DIR '/tmp' num2str(randi(99999)) '.sdcopen'];
    if move2dir==1
        movefile(fname,fn_dcmexp);
    else
        copyfile(fname,fn_dcmexp);
    end
end

end      % main function write_dicom.m


%% sub-functions
% converting p-file to dicom header
function hdr = sub_convert_header(hdr,h,nn,verb)

cnv = sub_get_hdr_conversion;     % load conversion cell array

for l=1:length(cnv)
    if ~isfield(hdr,cnv{l}{1})    % is dicom header field already assigned?
        if isfield(h,cnv{l}{2})   % is p-file header field exsting?
            if isfield(h.(cnv{l}{2}),cnv{l}{3})
                val = h.(cnv{l}{2}).(cnv{l}{3});
                if length(cnv{l})>3
                    if strcmp(cnv{l}{4},'char')
                        if ~ischar(val), val = deblank(char(val)); end
                        cntrl = isstrprop(val,'cntrl');
                        if any(cntrl)
                            fprintf('Warning: field(=%s): val(=%s) contains cntrl char\n',...
                                cnv{l}{1},val);
                            if cntrl(1)==1
                                val(cntrl) = ' ';
                            else
                                tmp = find(cntrl);
                                val = val(1:(tmp(1)-1));
                            end
                            fprintf('\tnew val = %s\n',val);
                        end
                    else
                        warning('~strcmp(cnv{l}{4},''char'')')
                    end
                end
                hdr.(cnv{l}{1}) = val;
            else
                warning('header field ''%s.%s'' not found',...
                    cnv{l}{2},cnv{l}{3});
            end
        else
            warning('header field ''%s'' not found',cnv{l}{2});
        end
    else
        if verb>1, fprintf('hdr.%s existing\n',cnv{l}{1}); end
    end       % if ~isfield(hdr,cnv{l}{1})
end           % for l=1:length(cnv)

%% date and time
hdr = sub_add_datetime(hdr,'StudyDate','StudyTime',h.exam.ex_datetime,verb);
hdr = sub_add_datetime(hdr,'SeriesDate','SeriesTime',h.series.se_datetime,verb);
hdr = sub_add_datetime(hdr,'ImageDate','ImageTime',h.image.im_datetime,verb);
hdr = sub_add_datetime(hdr,'AcquisitionDate','AcquisitionTime',...
    h.image.im_datetime,verb);

%% add text fields
hdr = sub_add2hdr(hdr,'ImageType','ORIGINAL\PRIMARY\OTHER',verb);
hdr = sub_add2hdr(hdr,'Manufacturer','GE MEDICAL SYSTEMS',verb);
hdr = sub_add2hdr(hdr,'ProductID','SIGNA',verb);
hdr = sub_add2hdr(hdr,'PatientsWeight',h.exam.patweight/1000,verb);
hdr = sub_add2hdr_dim(hdr,'MRAcquisitionType',h.image.imode,verb);
hdr = sub_add2hdr(hdr,'AngioFlag','N',verb);
hdr = sub_add2hdr(hdr,'RepetitionTime',h.image.tr/1000,verb);
hdr = sub_add2hdr(hdr,'EchoTime',h.rdb_hdr.te/1000,verb);
hdr = sub_add2hdr(hdr,'InversionTime',h.image.ti/1000,verb);
hdr = sub_add2hdr(hdr,'ImagingFrequency',h.rdb_hdr.ps_mps_freq*1d-7,verb);
hdr = sub_add2hdr_nuc(hdr,'ImagedNucleus',h.image.specnuc,verb);
b0 = round(h.rdb_hdr.ps_mps_freq/gyrogamma(h.image.specnuc)*2*pi)/10;
hdr = sub_add2hdr(hdr,'MagneticFieldStrength',b0,verb);
if isfield(h.rdb_hdr,'oversamplingfactor')
    hdr = sub_add2hdr(hdr,'PercentSampling',h.rdb_hdr.oversamplingfactor*100,verb);
end
hdr = sub_add2hdr(hdr,'StudyID',num2str(h.exam.ex_no),verb);
% hdr = sub_add2hdr(hdr,'ScannerStudyID',num2str(h.exam.ex_no),verb);
hdr = sub_add2hdr(hdr,'ManufacturersModelName','MNS Research Pack',verb);
hdr = sub_add2hdr(hdr,'SpecificCharacterSet','ISO_IR 100',verb);
hdr = sub_add2hdr(hdr,'PatientsAge',sprintf('%.3dY',h.exam.patage),verb);


%% adjust geometry + #images
hdr = sub_add2hdr_orientation(hdr,h,verb);
hdr = sub_add2hdr(hdr,'Rows',nn(1),verb);
hdr = sub_add2hdr(hdr,'Columns',nn(2),verb);
hdr = sub_add2hdr(hdr,'ImagesInAcquisition',nn(3)*nn(4)*nn(5),verb);
hdr = sub_add2hdr(hdr,'LocationsInAcquisition',1,verb);
hdr = sub_add2hdr(hdr,'NumberOfAcquisitions',1,verb);
hdr = sub_add2hdr(hdr,'LastImageNumberUsed',nn(3)*nn(4)*nn(5),verb);
hdr = sub_add2hdr(hdr,'ImageDimensionX',nn(1),verb);
hdr = sub_add2hdr(hdr,'ImageDimensionY',nn(2),verb);
hdr = sub_add2hdr(hdr,'PixelSpacing',h.rdb_hdr.fov./nn(1:2),verb);
if strcmpi(hdr.MRAcquisitionType,'3')
    slthick = h.rdb_hdr.fov/nn(3);
    spacing = slthick;
else
    slthick = h.image.slthick;
    spacing = h.image.slthick+h.image.scanspacing;
end
hdr = sub_add2hdr(hdr,'SliceThickness',slthick,verb);
hdr = sub_add2hdr(hdr,'SpacingBetweenSlices',spacing,verb);
hdr = sub_add2hdr_sex(hdr,h);


%% add/convert dicom unique identifiers
if false
    rootuid = uncompress_uid(h.series.equipmnt_uid,21,verb);
    rootuid = rootuid(1:21);
    if ~strcmp(rootuid(end),'.')
        warning('rootuid (=''%s'') not ending with .; adding',rootuid);
        rootuid = [rootuid '.'];
    end
    if ~strcmp(rootuid(1:15),'1.2.840.113619.')
        warning('rootuid (=''%s'') not a GE root ID (=''1.2.840.113619.'')',...
            rootuid)
    end
end

hdr = sub_add2hdr_uid(hdr,'StudyInstanceUID',h.exam.study_uid,verb);
hdr = sub_add2hdr_uid(hdr,'FrameOfReferenceUID',h.series.landmark_uid,verb);
hdr = sub_add2hdr_uid(hdr,'EquipmentUID',h.series.equipmnt_uid,verb);
% hdr = sub_add2hdr_uid(hdr,'SeriesInstanceUID',h.series.series_uid,verb);
if ~isfield(hdr,'SeriesInstanceUID')
    tmp = uncompress_uid(h.series.series_uid,verb);
    hdr.SeriesInstanceUID = [tmp '.' num2str(randi(999))];
end

end           % sub_convert_header


%% cell array with info converting p-file to dicom header
function hdr_conversion = sub_get_hdr_conversion()

% fieldnames: dicom header - p-file header
hdr_conversion = {...
    {'StudyDescription','exam','ex_desc','char'},...
    {'SeriesDescription','series','se_desc','char'},...
    {'Modality','exam','ex_typ','char'},...
    {'InstitutionName','exam','hospname','char'},...
    {'StationName','exam','ex_sysid','char'},...
    {'SuiteID','exam','ex_suid','char'},...
    {'ServiceID','exam','service_id','char'},...
    {'ActualSeriesDataTimeStamp','series','se_actual_dt'},...
    {'PatientsName','exam','patnameff','char'},...
    {'PatientID','exam','patidff','char'},...
    {'PatientsBirthDate','exam','dateofbirth','char'},...
    ...%{'PatientsSex','exam','patsex','char',},...
    {'ImageActualDate','image','im_datetime'},...
    {'NumberOfAverages','image','nex'},...
    {'EchoNumbers','rdb_hdr','nechoes'},...
    {'EchoTrainLength','rdb_hdr','etl'},...
    {'DeviceSerialNumber','exam','uniq_sys_id','char'},...
    {'ProtocolName','series','prtcl','char'},...
    {'HeartRate','image','hrtrate'},...
    {'TriggerWindow','image','trgwindow'},...
    {'ReconstructionDiameter','image','dfov'},...
    {'ReceivingCoil','image','cname','char'},...
    {'FlipAngle','image','mr_flip'},...
    {'SAR','image','saravg'},...
    {'SeriesPlane','series','se_plane'},...
    {'DisplayFieldOfView','image','dfov'},...
    {'AcquisitionDuration','image','sctime'},...
    {'NumberOfEchoes','rdb_hdr','nechoes'},...
    {'APSCenterFrequency','rdb_hdr','ps_aps_freq'},...
    {'InternalPulseSeqName','image','psd_iname','char'},...
    {'RawDataRunNumber','image','rawrunnum'},...
    {'SeriesNumber','series','se_no'},...
    {'PlaneType','image','plane'},...
    };
end


%% add to header if field not yet existing
function hdr = sub_add2hdr(hdr,fname,value,verb)
if ~isfield(hdr,fname)
    hdr.(fname) = value; 
else
    if verb>1, fprintf('hdr.%s existing\n',fname); end
end

end           % sub_add2header


%% add/convert time
function [hdr] = sub_add_datetime(hdr,dfname,tfname,dt,verb)
% function datetime available >=R2014b
% d_str = char(datetime(dt,'ConvertFrom','posixtime','Format','yyyyMMdd'));
% t_str = char(datetime(dt,'ConvertFrom','posixtime','Format','HHmmss'));
d_str = datestr(dt/86400+datenum('01-Jan-1970'),'yyyymmdd');
t_str = datestr(dt/86400+datenum('01-Jan-1970'),'HHMMss');
if ~isfield(hdr,dfname)
    hdr.(dfname) = d_str; 
else
    if verb>1, fprintf('hdr.%s existing\n',dfname); end
end
if ~isfield(hdr,tfname)
    hdr.(tfname) = t_str; 
else
    if verb>1, fprintf('hdr.%s existing\n',tfname); end
end

end      % sub_convert_datetime


%% add/translate nucleus
function hdr = sub_add2hdr_nuc(hdr,fname,nucleus,verb)
switch nucleus
    case 1,  nuc = '1H';
    case 2,  nuc = '2D';
    case 3,  nuc = '3He';
    case 7,  nuc = '7Li';
    case 11, nuc = '11B';
    case 13, nuc = '13C';
    case 17, nuc = '17O';
    case 19, nuc = '19F';        
    case 23, nuc = '23Na';
    case 29, nuc = '29Si';
    case 31, nuc = '31P';
    case 129, nuc = '129Xe';
    otherwise
        warning('nucleus (=%g) not in list; setting to ''1H''',nucleus);
end
if ~isfield(hdr,fname)
    hdr.(fname) = nuc; 
else
    if verb>1, fprintf('hdr.%s existing\n',fname); end
end

end           % sub_add2header


%% add patient position/orientation to header
function hdr = sub_add2hdr_orientation(hdr,h,verb)

if ~isfield(hdr,'PatientPosition')
    pos_str = '';
    switch h.series.entry
        case 1, pos_str = [pos_str 'HF'];
        case 2, pos_str = [pos_str 'FF'];
        otherwise
            warning('unknown h.series.entry (=%g)',h.series.entry);
    end
     switch h.series.position
         case 1, pos_str = [pos_str 'S'];
         case 2, pos_str = [pos_str 'P'];
         case 4, pos_str = [pos_str 'DL'];
         case 8, pos_str = [pos_str 'DR'];
         otherwise
             warning('unknown h.series.position position (=%g)',...
                 h.series.position);
    end
    hdr.PatientPosition = pos_str; 
else
    if verb>1, fprintf('hdr.%s existing\n',fname); end
end

if ~isfield(hdr,'ImageOrientationPatient')
    hdr.ImageOrientationPatient = sub_get_iop(h);
else
    if verb>1, fprintf('hdr.%s existing\n',fname); end
end

end           % sub_add2header_orientation


%% add UID to header
function hdr = sub_add2hdr_uid(hdr,fname,value,verb)
if ~isfield(hdr,fname)
    hdr.(fname) = uncompress_uid(value,[],verb); 
else
    if verb>1, fprintf('hdr.%s existing\n',fname); end
end

end           % sub_add2header_uid


%% add dimension to header
function hdr = sub_add2hdr_dim(hdr,fname,imode,verb)
dim_str = '2D';
switch imode
    case {1,4,6} 
    case {2,9}, dim_str = '3D';
    otherwise, warning('unknown imode(=%g): using ''2D''',imode);
end
if ~isfield(hdr,fname)
    hdr.(fname) = dim_str;
else
    if verb>1, fprintf('hdr.%s existing\n',fname); end
end

end           % sub_add2header_dim


%% create dicom UID
function uid = sub_create_uid(rootuid)
tmp = dicomuid;
uid = [rootuid tmp(24:length(tmp))];

end           % sub_get_uid


%% get RAS letter for orientation
function RAS = sub_get_ras(ori,loc)
switch ori
    case 'a', if loc>=0, RAS = 'S'; else, RAS = 'I'; end
    case 'c', if loc>=0, RAS = 'A'; else, RAS = 'P'; end
    case 's', if loc>=0, RAS = 'R'; else, RAS = 'L'; end
    case 'o', if loc>=0, RAS = '+'; else, RAS = '-'; end
end
end           % sub_get_ras


%% ImageOrientationPatient (from gemrraw2dcm.c)
function iop = sub_get_iop(h)
fov = h.rdb_hdr.fov;
iop = [...
    -(h.image.trhc_R-h.image.tlhc_R) ...
    -(h.image.trhc_A-h.image.tlhc_A) ...
    +(h.image.trhc_S-h.image.tlhc_S) ...
    -(h.image.brhc_R-h.image.trhc_R) ...
    -(h.image.brhc_A-h.image.trhc_A) ...
    +(h.image.brhc_S-h.image.trhc_S)];

% iop = +(iop>1d-10) -(iop<-1d-10);       % normalise
% iop = iop.';
iop = iop.'/fov;
iop = round(iop*1d5)/1d5;
iop(abs(iop)<1d-10) = 0;

end           % sub_get_iop


%% ImagePositionPatient (from gemrraw2dcm.c)
function ipp = sub_get_ipp(h,pixsize)
% DICOM standard is LPS not RAS -> converting
% same as ipp = [h.image.normal_L ; h.image.normal_P ; h.image.normal_S]
% but normal_X only in rdb header>=20.007

iop = sub_get_iop(h);
tlhc_C_l = (pixsize(1)*iop(1) + pixsize(2)*iop(4))/2.0 - h.image.tlhc_R;
tlhc_C_p = (pixsize(1)*iop(2) + pixsize(2)*iop(5))/2.0 - h.image.tlhc_A;
tlhc_C_s = (pixsize(1)*iop(3) + pixsize(2)*iop(6))/2.0 + h.image.tlhc_S;
ipp = [tlhc_C_l ; tlhc_C_p ; tlhc_C_s];

end     % sub_get_ipp

%% patient sex
function hdr = sub_add2hdr_sex(hdr,h)
fld = 'PatientsSex';
if ~isfield(hdr,fld)
    if isfield(h.exam,'patsex')
        patsex = h.exam.patsex;
        if isnumeric(patsex)
            switch patsex
                case 0, str = '';
                case 1, str = 'M';
                case 2, str = 'F';
                case 3, str = 'O';
                otherwise
                    str = '';
                    warning('patsex(=%g) unknown',patsex);
            end
            hdr.(fld) = str;
        else
            warning('patsex not a number');
        end
    else
        warning('field h.exam.patsex not found');
    end
end

end      % sub_add2hdr_sex
