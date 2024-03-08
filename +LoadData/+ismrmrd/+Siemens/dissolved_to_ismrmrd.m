function dissolved_to_ismrmrd(Xe_file,Subj_ID,mrdfile)
%%
%   mrdfile             -> Name of MRD file to write (will default to
%                          'Subj_ID_dixon.h5 if not provided)
%
%   Xe_file             -> Path to Dissolved Phase twix file (will prompt to select a file in a UI if not provided) 
%
%   Subj_ID             -> Subject ID (string) ((will prompt user for
%                          subject ID if not provided)
%

%Make this so it can be done with just the name of the function if desired:
if nargin == 0
    [Xe_file,mypath] = uigetfile('*.dat','Select Gas Exchange Raw Data file');
    Subj_ID = inputdlg('Subject ID','Input Subject ID',[1 50]); % Or maybe better to prompt user to type in Subject ID?
    Subj_ID = Subj_ID{1};
    mrdfile = fullfile(mypath,[Subj_ID '_dixon.h5']);
elseif nargin == 1
    [mypath,~,~] = fileparts(Xe_file); 
    Subj_ID = inputdlg('Subject ID','Input Subject ID',[1 50]); % Or maybe better to prompt user to type in Subject ID?
    Subj_ID = Subj_ID{1};
    mrdfile = fullfile(mypath,[Subj_ID '_dixon.h5']);
elseif nargin == 2
    [mypath,~,~] = fileparts(Xe_file); 
    mrdfile = fullfile(mypath,[Subj_ID '_dixon.h5']);
end

%Make sure file has .h5 extension
if ~contains(mrdfile,'.')
    mrdfile = [mrdfile '.h5'];
end

%% Read in twix file:
Xe_Dat_twix = LoadData.ismrmrd.DataImport.mapVBVD(Xe_file,'ignoreSeg');
Xe_Dat_twix.flagIgnoreSeg = 1;
Xe_Dat_twix.image.flagIgnoreSeg = 1;
Xe_Dat_twix.image.flagAverageReps = 1;
Xe_Dat_twix.flagAverageReps = 1;
%% Function to pull data from twix files of standard methods (Mugler and Niedbalski-built) at KUMC 
[Dis_Fid,Gas_Fid,Dis_Traj,Gas_Traj,Params,Post_Cal] = LoadData.ismrmrd.DataImport.pull_dis_data(Xe_file);

%% Should have all data, can move forward with writing: 
%% Create an empty ismrmrd dataset
if exist(mrdfile,'file')
    delete(mrdfile);
end

% [~,h] = read_p(dissolvedarchive);
% load(dissolvedmat);
% load(calmat,'te90');
% load(calmat,'targetAX');

dset = LoadData.ismrmrd.Dataset(mrdfile);

%% Need to handle the fact that for Niedbalski Single-Breath Method, Dis and Gas data aren't same size;
if size(Gas_Fid,1) ~= size(Dis_Fid,1)
    tmp = zeros(size(Gas_Fid));
    tmp(tmp==0) = nan;
    tmp(1:size(Dis_Fid,1),:) = Dis_Fid;
    Dis_Fid = tmp;
    clear tmp;
    
    tmp = zeros(size(Gas_Traj));
    tmp(tmp==0) = nan;
    tmp(:,1:size(Dis_Traj,2),:) = Dis_Traj;
    Dis_Traj = tmp;
    clear tmp;
end


nX = size(Gas_Fid,1);
nY = size(Gas_Fid,2);
tsp = Params.Dwell * 1e6; %need this in us %1e6/bw;
acqblock = LoadData.ismrmrd.Acquisition(2*nY);

acqblock.head.version(:) = 1;
acqblock.head.number_of_samples(:) = nX;
acqblock.head.sample_time_us(:) = tsp;
acqblock.head.active_channels(:) = 1;
acqblock.head.trajectory_dimensions = 3*ones(1,2*nY);

K = zeros(nX,2*nY);
T = zeros(3,nX,2*nY);
K(:,1:2:end) = Gas_Fid;
K(:,2:2:end) = Dis_Fid;

T(1,:,1:2:end) = Gas_Traj(1,:,:);
T(2,:,1:2:end) = Gas_Traj(2,:,:);
T(3,:,1:2:end) = Gas_Traj(3,:,:);
T(1,:,2:2:end) = Dis_Traj(1,:,:);
T(2,:,2:2:end) = Dis_Traj(2,:,:);
T(3,:,2:2:end) = Dis_Traj(3,:,:);

for acqno = 1:2*nY
    acqblock.head.scan_counter(acqno) = acqno-1;
    acqblock.head.idx.kspace_encode_step_1(acqno) = acqno-1;
    acqblock.head.idx.repetition(acqno) = 0;
    %set contrast to differentiate dissolved/gas. 0 = gas. 1 = dissolved
    if mod(acqno,2) == 0
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
    elseif acqno==size(K,2)
        acqblock.head.flagSet('ACQ_LAST_IN_ENCODE_STEP1', acqno);
        acqblock.head.flagSet('ACQ_LAST_IN_SLICE', acqno);
        acqblock.head.flagSet('ACQ_LAST_IN_REPETITION', acqno);
    end

    % fill the data
    acqblock.data{acqno} = squeeze(K(:,acqno));
    acqblock.traj{acqno} = squeeze(T(:,:,acqno));
end

% Append the acquisition block
dset.appendAcquisition(acqblock);

%% If available, upload the post-imaging spectra for calibration - This needs to be fixed. Scrub for now

% if ~isnan(sum(Post_Cal))
%     acqblock = ismrmrd.Acquisition(1);
% 
%     acqblock.head.version(:) = 1;
%     acqblock.head.number_of_samples(:) = length(Post_Cal);
%     acqblock.head.sample_time_us(:) = tsp;
%     acqblock.head.active_channels(:) = 1;
%     %acqblock.head.trajectory_dimensions = 3*ones(1,2*nY);
%     acqblock.head.measurement_uid = 1;
% 
%     acqblock.head.scan_counter(1) = 0;
%     acqblock.head.idx.kspace_encode_step_1(1) = 0;
%     acqblock.head.idx.repetition(1) = 0;
%     %appended spectrum acquired on dissolved frequency, so set contrast to
%     %1
%     acqblock.head.idx.contrast(1) = 1;
% 
%     acqblock.head.flagSet('ACQ_FIRST_IN_ENCODE_STEP1', 1);
%     acqblock.head.flagSet('ACQ_FIRST_IN_SLICE', 1);
%     acqblock.head.flagSet('ACQ_FIRST_IN_REPETITION', 1);
%     acqblock.head.flagSet('ACQ_LAST_IN_ENCODE_STEP1', 1);
%     acqblock.head.flagSet('ACQ_LAST_IN_SLICE', 1);
%     acqblock.head.flagSet('ACQ_LAST_IN_REPETITION', 1);
% 
%     acqblock.data{1} = squeeze(Post_Cal);
%     acqblock.traj{1} = [];
% end
% dset.appendAcquisition(acqblock);

%%
header = [];

% Experimental Conditions (Required)
header.experimentalConditions.H1resonanceFrequency_Hz = 128000000; % 3T

% Acquisition System Information (Optional)
header.acquisitionSystemInformation.receiverChannels = 1;
header.acquisitionSystemInformation.systemFieldStrength_T = Xe_Dat_twix.hdr.Dicom.flMagneticFieldStrength;
header.acquisitionSystemInformation.systemVendor = Xe_Dat_twix.hdr.Dicom.Manufacturer;
header.acquisitionSystemInformation.systemModel = Xe_Dat_twix.hdr.Dicom.ManufacturersModelName;
header.acquisitionSystemInformation.institutionName = Xe_Dat_twix.hdr.Dicom.InstitutionName;

Y = str2double(Params.scandatestr(1:4));
M = str2double(Params.scandatestr(6:7));
D = str2double(Params.scandatestr(8:9));

header.measurementInformation.scandate = string(datetime(Y,M,D));%datetimeParams.scandatestr;% ['20' h.rdb_hdr.scan_date(end-1:end) '-' h.rdb_hdr.scan_date(1:2) '-' h.rdb_hdr.scan_date(4:5)];
header.measurementInformation.patientPosition = '';

%header.studyInformation.studyDate = ['20' h.rdb_hdr.scan_date(end-1:end) '-' h.rdb_hdr.scan_date(1:2) '-' h.rdb_hdr.scan_date(4:5)];
header.subjectInformation.patientID = Subj_ID;


header.sequenceParameters.TR = Params.TR;
header.sequenceParameters.flipAngle_deg(1) = Params.GasFA;
header.sequenceParameters.flipAngle_deg(2) = Params.DisFA;
header.sequenceParameters.TE = Params.TE;%te90 + 40e-6;

% The Encoding (Required) - Hardcode to standards
header.encoding.trajectory = 'radial';
header.encoding.encodedSpace.fieldOfView_mm.x = Params.GE_FOV;
header.encoding.encodedSpace.fieldOfView_mm.y = Params.GE_FOV;
header.encoding.encodedSpace.fieldOfView_mm.z = Params.GE_FOV;
header.encoding.encodedSpace.matrixSize.x = Params.imsize;
header.encoding.encodedSpace.matrixSize.y = Params.imsize;
header.encoding.encodedSpace.matrixSize.z = Params.imsize;
% Recon Space
% (in this case same as encoding space)
header.encoding.reconSpace = header.encoding.encodedSpace;
% Encoding Limits
header.encoding.encodingLimits.kspace_encoding_step_0.minimum = 0;
header.encoding.encodingLimits.kspace_encoding_step_0.maximum = size(K,2)-1;
header.encoding.encodingLimits.kspace_encoding_step_0.center = floor(size(K,2)/2);
header.encoding.encodingLimits.kspace_encoding_step_1.minimum = 0;
header.encoding.encodingLimits.kspace_encoding_step_1.maximum = size(K,1)-1;
header.encoding.encodingLimits.kspace_encoding_step_1.center = floor(size(K,1)/2);
header.encoding.encodingLimits.repetition.minimum = 0;
header.encoding.encodingLimits.repetition.maximum = 0;
header.encoding.encodingLimits.repetition.center = 0;

% Custom trajectory parameters

a(1).name = "rampTime";
a(1).value = Xe_Dat_twix.hdr.MeasYaps.sWipMemBlock.alFree{8};
up(1).name = "xe_center_frequency";
up(1).value = Xe_Dat_twix.hdr.Dicom.lFrequency;
up(2).name = "xe_dissolved_offset_frequency";
up(2).value = Xe_Dat_twix.hdr.MeasYaps.sWipMemBlock.alFree{5};

header.encoding.trajectoryDescription.identifier = "Custom Trajectory Info";
header.encoding.trajectoryDescription.userParameterDouble = a;
header.encoding.trajectoryDescription.userParameterLong = [];
header.userParameters.userParameterLong = up;

%% Serialize and write to the data set
xmlstring = LoadData.ismrmrd.xml.serialize(header);
dset.writexml(xmlstring);

dset.close();
