function calibration_to_ismrmrd(mrdfile,MainInput)
%%
%   mrdfile -> Name of MRD file to write
%
%   calfile -> Scan Archive file generated from calibration scan
%               (e.g. ScanArchive*.h5)

%% Create an empty ismrmrd dataset
% if exist(mrdfile,'file')
%     error(['File ' mrdfile ' already exists.  Please remove first'])
% end
cd(MainInput.CalDataLocation)
dset = LoadData.ismrmrd.Dataset(mrdfile);

[data,h] = LoadData.ismrmrd.GE.Functions.read_p(MainInput.CalFileName);

cf = h.psc.mps_freq/10;

% The path here should be the path to the local calibration waveform file
directory = 'C:\XIPline\GE\waveforms\xe_calibration';
% Define file patterns to search for
filePatterns = {'*freq.fdl', '*TR.fdl'};
filePaths = {};
for i = 1:length(filePatterns)
    fileList = dir(fullfile(directory, filePatterns{i}));
    if isempty(fileList)
        fprintf('No file matching %s found in the directory.\n', filePatterns{i});
    else
        % Append each matching file's full path to the list
        for j = 1:length(fileList)
            filePaths{end+1} = fullfile(directory, fileList(j).name);
        end
    end
end

% Display the full paths
if isempty(filePaths)
    disp('No matching files found.');
else
    fprintf('Matching files found:\n');
    disp(filePaths');
end

% fdl_freq_file = 'D:\Github\XIPline\+LoadData\+ismrmrd\+GE\calibration3D_129xe_fov400_nsp256_intlv600_kdt40_dur10p2cal_freq.fdl';
fdl_freq_file = filePaths{1};
tmp = LoadData.ismrmrd.GE.Functions.read_fdl(fdl_freq_file);
freq_off = tmp(1); 
nDis = sum(tmp > 0); clear tmp;

% fdl_tr_file = 'D:\Github\XIPline\+LoadData\+ismrmrd\+GE\calibration3D_129xe_fov400_nsp256_intlv600_kdt40_dur10p2cal_TR.fdl';
fdl_tr_file = filePaths{2};
tmp = LoadData.ismrmrd.GE.Functions.read_fdl(fdl_tr_file);
tr_gas = tmp(end); clear tmp;

nX = size(data,2);
nY = size(data,1);
tsp = 1e6/2/h.rdb_hdr.bw/1000;
acqblock = LoadData.ismrmrd.Acquisition(nY);

acqblock.head.version(:) = 1;
acqblock.head.number_of_samples(:) = nX;
acqblock.head.sample_time_us(:) = tsp;
acqblock.head.active_channels(:) = 1;

for acqno = 1:nY
    acqblock.head.scan_counter(acqno) = acqno-1;
    acqblock.head.idx.kspace_encode_step_1(acqno) = acqno-1;
    acqblock.head.idx.repetition(acqno) = 0;
    
    acqblock.head.measurement_uid(acqno) = 0;  % no bonus spectra, set to 1 if bonus spectra
    
    if acqno <= nDis 
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
    elseif acqno==size(data,1)
        acqblock.head.flagSet('ACQ_LAST_IN_ENCODE_STEP1', acqno);
        acqblock.head.flagSet('ACQ_LAST_IN_SLICE', acqno);
        acqblock.head.flagSet('ACQ_LAST_IN_REPETITION', acqno);
    end

    % fill the data
    acqblock.data{acqno} = squeeze(data(acqno,:));
end

% Append the acquisition block
dset.appendAcquisition(acqblock);

header = [];

% Experimental Conditions (Required)
header.experimentalConditions.H1resonanceFrequency_Hz = 128000000; % 3T

% Acquisition System Information (Optional)
header.acquisitionSystemInformation.receiverChannels = 1;
header.acquisitionSystemInformation.systemFieldStrength_T = 3.0;
header.acquisitionSystemInformation.systemVendor = 'GE';
header.acquisitionSystemInformation.systemModel = 'Premiere';
header.acquisitionSystemInformation.institutionName = 'Robarts Research Institute';

header.measurementInformation.patientPosition = '';

header.studyInformation.studyDate = ['20' h.rdb_hdr.scan_date(end-1:end) '-' h.rdb_hdr.scan_date(1:2) '-' h.rdb_hdr.scan_date(4:5)];

header.subjectInformation.patientID = h.exam.patnameff;
header.subjectInformation.patientWeight_kg = h.exam.patweight/1000;

if ~isempty(h.exam.dateofbirth)
    header.subjectInformation.patientBirthdate = [h.exam.dateofbirth(1:4) '-' h.exam.dateofbirth(5:6) '-' h.exam.dateofbirth(7:8)];
end

header.sequenceParameters.TR(1) = h.image.tr*1e-3 + tr_gas*1e-3;
header.sequenceParameters.TR(2) = h.image.tr*1e-3;
header.sequenceParameters.flipAngle_deg(1) = 10;
header.sequenceParameters.flipAngle_deg(2) = h.image.mr_flip;
header.sequenceParameters.TE = h.rdb_hdr.te*1e-3;

% The Encoding (Required)
header.encoding.trajectory = 'cartesian';
header.encoding.encodedSpace.fieldOfView_mm.x = 400;
header.encoding.encodedSpace.fieldOfView_mm.y = 400;
header.encoding.encodedSpace.fieldOfView_mm.z = 30;
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
header.encoding.trajectoryDescription.identifier = "Custom Trajectory Info";
header.encoding.trajectoryDescription.userParameterDouble = [];
header.encoding.trajectoryDescription.userParameterLong(1).name = "ramp_time";
header.encoding.trajectoryDescription.userParameterLong(1).value = 0;

header.userParameters.userParameterLong(1).name = "xe_center_frequency";
header.userParameters.userParameterLong(1).value = cf;
header.userParameters.userParameterLong(2).name = "xe_dissolved_offset_frequency";
header.userParameters.userParameterLong(2).value = freq_off;


%% Serialize and write to the data set
xmlstring = LoadData.ismrmrd.xml.serialize(header);
dset.writexml(xmlstring);


dset.close();

