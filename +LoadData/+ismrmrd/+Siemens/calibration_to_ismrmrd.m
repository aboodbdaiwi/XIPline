function calibration_to_ismrmrd(Cal_file,Subj_ID,mrdfile)
%%
%   mrdfile             -> Name of MRD file to write (will default to
%                          'Subj_ID_cali.h5 if not provided)
%
%   Cal_file             -> Path to Calibration twix file (will prompt to select a file in a UI if not provided) 
%
%   Subj_ID             -> Subject ID (string) (will prompt user for
%                          subject ID if not provided)
%Make this so it can be done with just the name of the function if desired:
if nargin == 0
    [Cal_file,mypath] = uigetfile('*.dat','Select Calibration Raw Data file');
    Cal_file = fullfile(mypath,Cal_file);
    Subj_ID = inputdlg('Subject ID','Input Subject ID',[1 50]); % Or maybe better to prompt user to type in Subject ID?
    Subj_ID = Subj_ID{1};
    mrdfile = fullfile(mypath,[Subj_ID '_calibration.h5']);
elseif nargin == 1
    [mypath,~,~] = fileparts(Cal_file); 
    Subj_ID = inputdlg('Subject ID','Input Subject ID',[1 50]); % Or maybe better to prompt user to type in Subject ID?
    Subj_ID = Subj_ID{1};
    mrdfile = fullfile(mypath,[Subj_ID '_calibration.h5']);
elseif nargin == 2
    [mypath,~,~] = fileparts(Cal_file); 
    mrdfile = fullfile(mypath,[Subj_ID '_calibration.h5']);
end

%Make sure file has .h5 extension
if ~contains(mrdfile,'.')
    mrdfile = [mrdfile '.h5'];
end

%% Create an empty ismrmrd dataset
if exist(mrdfile,'file')
    delete(mrdfile);
end

dset = LoadData.ismrmrd.Dataset(mrdfile);

%Read in calibration data
Cal_Dat_twix = LoadData.ismrmrd.DataImport.mapVBVD(Cal_file,'ignoreSeg');

% Get various scan parameters
dwell_time = Cal_Dat_twix.hdr.MeasYaps.sRXSPEC.alDwellTime{1,1};
dwell_time=dwell_time*1E-3; %Convert from ns to us. 
% nFids = Cal_Dat_twix.hdr.Config.NRepMeas;
% nPts = Cal_Dat_twix.hdr.Config.VectorSize;
freq = Cal_Dat_twix.hdr.Dicom.lFrequency;
%imsize = Cal_Dat_twix.hdr.MeasYaps.sKSpace.lBaseResolution;
TR = (Cal_Dat_twix.hdr.MeasYaps.alTR{1}/1000);
TE = (Cal_Dat_twix.hdr.MeasYaps.alTE{1}/1000);
GasFA = Cal_Dat_twix.hdr.MeasYaps.adFlipAngleDegree{2};
DisFA = Cal_Dat_twix.hdr.MeasYaps.adFlipAngleDegree{1};
FOV = Cal_Dat_twix.hdr.Config.ReadFoV;
%Does the location of frequency offset change with different versions?
freq_offset = Cal_Dat_twix.hdr.MeasYaps.sWipMemBlock.alFree{5};


% I think in an older version of John's code, this is in the 2 spot 
nDis = Cal_Dat_twix.hdr.Phoenix.sWipMemBlock.alFree{3};
if nDis == 1
    nDis = Cal_Dat_twix.hdr.Phoenix.sWipMemBlock.alFree{2}; %It seems like we should preserve this information.
end
try
    scanDate = Cal_Dat_twix.hdr.Phoenix.tReferenceImage0; 
    scanDate = strsplit(scanDate,'.');
    scanDate = scanDate{end};
    scanDateStr = [scanDate(1:4),'-',scanDate(5:6),'-',scanDate(7:8)];
   % scanDateStr = string(datetime(str2double(scanDate(1:4)),str2double(scanDate(5:6)),scanDate(7:8)));
catch
    scanDate = Cal_Dat_twix.hdr.Meas.SeriesLOID;
    scanDate = strsplit(scanDate,'.');
    scanDate = char(scanDate(10));
    scanDateStr = [scanDate(1:4),'-',scanDate(5:6),'-',scanDate(7:8)];
   % scanDateStr = string(datetime(str2double(scanDate(1:4)),str2double(scanDate(5:6)),scanDate(7:8)));
end


%Read data
data = squeeze(double(Cal_Dat_twix.image()));

%% write data in acq blocks - I'm inclined to separate gas and dissolved into two blocks and give a different measurement ID.
nX = size(data,1); % Number of points
nY = size(data,2); % Number of FIDs
%tsp = 1e6/2/h.rdb_hdr.bw/1000;
acqblock = LoadData.ismrmrd.Acquisition(nY);

acqblock.head.version(:) = 1;
acqblock.head.number_of_samples(:) = nX;
acqblock.head.sample_time_us(:) = dwell_time; %dwell time currently in s
acqblock.head.active_channels(:) = 1;

for acqno = 1:nY
    acqblock.head.scan_counter(acqno) = acqno-1;
    acqblock.head.idx.kspace_encode_step_1(acqno) = acqno-1;
    acqblock.head.idx.repetition(acqno) = 0;
    acqblock.head.idx.measurement_uid(acqno) = 0;
    
    %Set contrasts - Gas = 0, Dissolved = 1
    if acqno <= nDis+1 % first scan is noise, then dissolved, so need to add 1 to the number to get this right
        acqblock.head.idx.contrast(acqno) = 2;
    else
        acqblock.head.idx.contrast(acqno) = 1;
    end
    
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
    acqblock.data{acqno} = squeeze(data(:,acqno));
end

% Append the acquisition block
dset.appendAcquisition(acqblock);

header = [];

% Experimental Conditions (Required)
header.experimentalConditions.H1resonanceFrequency_Hz = 128000000; % 3T

% Acquisition System Information (Optional)
header.acquisitionSystemInformation.receiverChannels = 1;
header.acquisitionSystemInformation.systemFieldStrength_T = Cal_Dat_twix.hdr.Dicom.flMagneticFieldStrength;
header.acquisitionSystemInformation.systemVendor = Cal_Dat_twix.hdr.Dicom.Manufacturer;
header.acquisitionSystemInformation.systemModel = Cal_Dat_twix.hdr.Dicom.ManufacturersModelName;
header.acquisitionSystemInformation.institutionName = Cal_Dat_twix.hdr.Dicom.InstitutionName;

% scanDate = twix_obj.hdr.Phoenix.tReferenceImage0; 
% scanDate = strsplit(scanDate,'.');
% scanDate = scanDate{end};
% scanDateStr = [scanDate(1:4),'-',scanDate(5:6),'-',scanDate(7:8)];
header.studyInformation.studyDate = scanDateStr;
header.measurementInformation.patientPosition = '';
%header.studyInformation.studyDate = ['20' h.rdb_hdr.scan_date(end-1:end) '-' h.rdb_hdr.scan_date(1:2) '-' h.rdb_hdr.scan_date(4:5)];
header.subjectInformation.patientID = Subj_ID;

header.sequenceParameters.TR(1) = TR;
header.sequenceParameters.TR(2) = TR;
header.sequenceParameters.flipAngle_deg(1) = GasFA;
header.sequenceParameters.flipAngle_deg(2) = DisFA;
%header.sequenceParameters.flipAngle_deg(2) = h.image.mr_flip;
header.sequenceParameters.TE = TE;

% The Encoding (Required)
header.encoding.trajectory = 'cartesian';
header.encoding.encodedSpace.fieldOfView_mm.x = FOV;
header.encoding.encodedSpace.fieldOfView_mm.y = FOV;
header.encoding.encodedSpace.fieldOfView_mm.z = FOV;
header.encoding.encodedSpace.matrixSize.x = size(data,1);
header.encoding.encodedSpace.matrixSize.y = size(data,2);
header.encoding.encodedSpace.matrixSize.z = 1;
% Recon Space
% (in this case same as encoding space)
header.encoding.reconSpace = header.encoding.encodedSpace;
% Encoding Limits
header.encoding.encodingLimits.kspace_encoding_step_0.minimum = 0;
header.encoding.encodingLimits.kspace_encoding_step_0.maximum = size(data,1)-1;
header.encoding.encodingLimits.kspace_encoding_step_0.center = floor(size(data,1)/2);
header.encoding.encodingLimits.kspace_encoding_step_1.minimum = 0;
header.encoding.encodingLimits.kspace_encoding_step_1.maximum = size(data,2)-1;
header.encoding.encodingLimits.kspace_encoding_step_1.center = floor(size(data,2)/2);
header.encoding.encodingLimits.repetition.minimum = 0;
header.encoding.encodingLimits.repetition.maximum = 0;
header.encoding.encodingLimits.repetition.center = 0;

% Custom trajectory parameters

up(1).name = "xe_center_frequency";
up(1).value = freq;
up(2).name = "xe_dissolved_offset_frequency";
up(2).value = (freq_offset);

% Is there a place to write how many dissolved/gas fids there are?
header.userParameters.userParameterLong = up;


%% Serialize and write to the data set
xmlstring = LoadData.ismrmrd.xml.serialize(header);
dset.writexml(xmlstring);

dset.close();
