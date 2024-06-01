function ute_to_ismrmrd(ute_file,Subj_ID,mrdfile)
%%
%   mrdfile             -> Name of MRD file to write (will default to
%                          'Subj_ID_dixon.h5 if not provided)
%
%   ute_file             -> Path to 1H UTE (anatomic) twix file (will prompt to select a file in a UI if not provided) 
%
%   Subj_ID             -> Subject ID (string) ((will prompt user for
%                          subject ID if not provided)
%

%% Make this so it can be done with just the name of the function if desired:
if nargin == 0
    [ute_file,mypath] = uigetfile('*.dat','Select Anatomic UTE Raw Data file');
    Subj_ID = inputdlg('Subject ID','Input Subject ID',[1 50]); % Or maybe better to prompt user to type in Subject ID?
    Subj_ID = Subj_ID{1};
    mrdfile = fullfile(mypath,[Subj_ID '_proton.h5']);
elseif nargin == 1
    [mypath,~,~] = fileparts(ute_file); 
    Subj_ID = inputdlg('Subject ID','Input Subject ID',[1 50]); % Or maybe better to prompt user to type in Subject ID?
    Subj_ID = Subj_ID{1};
    mrdfile = fullfile(mypath,[Subj_ID '_proton.h5']);
elseif nargin == 2
    [mypath,~,~] = fileparts(ute_file); 
    mrdfile = fullfile(mypath,[Subj_ID '_proton.h5']);
end

%Make sure file has .h5 extension
if ~contains(mrdfile,'.')
    mrdfile = [mrdfile '.h5'];
end
%% Get UTE Data (trajectories and raw)
[data,traj] = LoadData.ismrmrd.DataImport.get_ute_data(ute_file);


%% Create an empty ismrmrd dataset
if exist(mrdfile,'file')
    delete(mrdfile)
end

% tmp = h5read('OUTPUT/MRI_Raw.h5','/Kdata/KData_E0_C0');
% K = complex(tmp.real,tmp.imag); clear tmp;
% kx = h5read('OUTPUT/MRI_Raw.h5','/Kdata/KX_E0');
% ky = h5read('OUTPUT/MRI_Raw.h5','/Kdata/KY_E0');
% kz = h5read('OUTPUT/MRI_Raw.h5','/Kdata/KZ_E0');


dset = LoadData.ismrmrd.Dataset(mrdfile);


nX = size(data,1);
nY = size(data,2);
acqblock = LoadData.ismrmrd.Acquisition(nY);

acqblock.head.version(:) = 1;
acqblock.head.number_of_samples(:) = nX;
acqblock.head.active_channels(:) = size(data,3);
acqblock.head.trajectory_dimensions = 3*ones(1,nY);

% T = zeros(3,nX,nY);
% 
% T(1,:,:) = kx;
% T(2,:,:) = ky;
% T(3,:,:) = kz;

for acqno = 1:nY
    acqblock.head.scan_counter(acqno) = acqno-1;
    acqblock.head.idx.kspace_encode_step_1(acqno) = acqno-1;
    acqblock.head.idx.repetition(acqno) = 0;
    acqblock.head.idx.contrast(acqno) = 0;
    % Set the flags
    acqblock.head.flagClearAll(acqno);
    if acqno == 1
        acqblock.head.flagSet('ACQ_FIRST_IN_ENCODE_STEP1', acqno);
        acqblock.head.flagSet('ACQ_FIRST_IN_SLICE', acqno);
        acqblock.head.flagSet('ACQ_FIRST_IN_REPETITION', acqno);
    elseif acqno==size(data,2)
        acqblock.head.flagSet('ACQ_LAST_IN_ENCODE_STEP1', acqno);
        acqblock.head.flagSet('ACQ_LAST_IN_SLICE', acqno);
        acqblock.head.flagSet('ACQ_LAST_IN_REPETITION', acqno);
    end

    % fill the data
    acqblock.data{acqno} = squeeze(data(:,acqno,:));
    acqblock.traj{acqno} = squeeze(traj(:,:,acqno));
end

% Append the acquisition block
dset.appendAcquisition(acqblock);


header = [];

% Experimental Conditions (Required)
header.experimentalConditions.H1resonanceFrequency_Hz = 128000000; % 3T

twix = LoadData.ismrmrd.DataImport.mapVBVD(ute_file);
if length(twix) > 1
    twix = twix{end};
end

% Acquisition System Information (Optional)
header.acquisitionSystemInformation.receiverChannels = size(data,2);
header.acquisitionSystemInformation.systemFieldStrength_T = twix.hdr.Dicom.flMagneticFieldStrength;
header.acquisitionSystemInformation.systemVendor = twix.hdr.Dicom.Manufacturer;
header.acquisitionSystemInformation.systemModel = twix.hdr.Dicom.ManufacturersModelName;
header.acquisitionSystemInformation.institutionName = twix.hdr.Dicom.InstitutionName;


scanDate = twix.hdr.Phoenix.tReferenceImage0; 
scanDate = strsplit(scanDate,'.');
scanDate = scanDate{end};
scanDateStr = [scanDate(1:4),'-',scanDate(5:6),'-',scanDate(7:8)];
%scanDateStr = string(datetime(str2double(scanDate(1:4)),str2double(scanDate(5:6)),scanDate(7:8)));

header.studyInformation.studyDate = scanDateStr;
header.subjectInformation.patientID = Subj_ID;

imsize = twix.hdr.MeasYaps.sKSpace.lBaseResolution;
FOV = twix.hdr.Config.ReadFoV;
% The Encoding (Required)
header.encoding.trajectory = 'radial';
header.encoding.encodedSpace.fieldOfView_mm.x = FOV;
header.encoding.encodedSpace.fieldOfView_mm.y = FOV;
header.encoding.encodedSpace.fieldOfView_mm.z = FOV;
header.encoding.encodedSpace.matrixSize.x = imsize;
header.encoding.encodedSpace.matrixSize.y = imsize;
header.encoding.encodedSpace.matrixSize.z = imsize;
% Recon Space
% (in this case same as encoding space)
header.encoding.reconSpace = header.encoding.encodedSpace;
% Encoding Limits
header.encoding.encodingLimits.kspace_encoding_step_0.minimum = 0;
header.encoding.encodingLimits.kspace_encoding_step_0.maximum = size(data,2)-1;
header.encoding.encodingLimits.kspace_encoding_step_0.center = floor(size(data,2)/2);
header.encoding.encodingLimits.kspace_encoding_step_1.minimum = 0;
header.encoding.encodingLimits.kspace_encoding_step_1.maximum = size(data,1)-1;
header.encoding.encodingLimits.kspace_encoding_step_1.center = floor(size(data,1)/2);
header.encoding.encodingLimits.repetition.minimum = 0;
header.encoding.encodingLimits.repetition.maximum = 0;
header.encoding.encodingLimits.repetition.center = 0;

%% Serialize and write to the data set
xmlstring = LoadData.ismrmrd.xml.serialize(header);
dset.writexml(xmlstring);


dset.close();
