function dissolved_to_ismrmrd(mrdfile,dissolvedarchive,dissolvedmat,calmat)
%%
%   mrdfile             -> Name of MRD file to write
%
%   dissolvedarchive    -> Scan Archive file generated from Dixon scan
%                         (e.g. ScanArchive*.h5)
%
%   dissolvedmat        -> Mat file genereated from prep_dixon_gx
%                         (e.g. Dixon_ScanArchive*.mat)
%
%   calmat              -> Mat file generated from recon_calibration
%                         (e.g ScanArchive*.mat, NOT Spect_ScanArchive.mat)

%% Create an empty ismrmrd dataset
if exist(mrdfile,'file')
    error(['File ' filename ' already exists.  Please remove first'])
end

[~,h] = LoadData.ismrmrd.GE.Functions.read_p(dissolvedarchive);
load(dissolvedmat);
load(calmat,'te90');
load(calmat,'targetAX');

% freqfdlpath = MainInput.freqfdlFullPath;
% The path here should be the path to the local calibration waveform file
directory = 'C:\XIPline\GE\waveforms\xe_dissolved';
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

fdl_freq_file = filePaths{1};
tmp = LoadData.ismrmrd.GE.Functions.read_fdl(fdl_freq_file);
% tmp = read_fdl('C:/Users/stcherne/Documents/Code_Repository/GE_2_ISMRMRD/DP/radial3D_129Xe_fov400_mtx64_intlv2000_kdt20_gmax33_smax147_dur4p6_coca_G_freq.fdl');
freq_off = tmp(2); clear tmp;

dset = LoadData.ismrmrd.Dataset(mrdfile);

nX = size(k_dp,1);
nY = size(k_dp,2);
tsp = 1e6/bw;
acqblock = LoadData.ismrmrd.Acquisition(2*nY);

acqblock.head.version(:) = 1;
acqblock.head.number_of_samples(:) = nX;
acqblock.head.sample_time_us(:) = tsp;
acqblock.head.active_channels(:) = 1;
acqblock.head.trajectory_dimensions = 3*ones(1,2*nY);

K = zeros(nX,2*nY);
T = zeros(3,nX,2*nY);
K(:,1:2:end) = k_gp;
K(:,2:2:end) = k_dp;

T(1,:,1:2:end) = x_gp;
T(2,:,1:2:end) = y_gp;
T(3,:,1:2:end) = z_gp;
T(1,:,2:2:end) = x_dp;
T(2,:,2:2:end) = y_dp;
T(3,:,2:2:end) = z_dp;

for acqno = 1:2*nY
    acqblock.head.scan_counter(acqno) = acqno-1;
    acqblock.head.idx.kspace_encode_step_1(acqno) = acqno-1;
    acqblock.head.idx.repetition(acqno) = 0;    
    
    if mod(acqno,2)
        acqblock.head.idx.contrast(acqno) = 2;  % even numbered views dissolved
    else
        acqblock.head.idx.contrast(acqno) = 1;  % odd numbered views gas
    end
    
    acqblock.head.measurement_uid(acqno) = 0;  % no bonus spectra, set to 1 if bonus spectra
      
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

header = [];

% Experimental Conditions (Required)
header.experimentalConditions.H1resonanceFrequency_Hz = 128000000; % 3T

% Acquisition System Information (Optional)
header.acquisitionSystemInformation.receiverChannels = 1;
header.acquisitionSystemInformation.systemFieldStrength_T = 3.0;
header.acquisitionSystemInformation.systemVendor = 'GE';
header.acquisitionSystemInformation.systemModel = 'Premiere';
header.acquisitionSystemInformation.institutionName = 'University of Iowa';

header.measurementInformation.patientPosition = '';

header.studyInformation.studyDate = ['20' h.rdb_hdr.scan_date(end-1:end) '-' h.rdb_hdr.scan_date(1:2) '-' h.rdb_hdr.scan_date(4:5)];

header.subjectInformation.patientID = h.exam.patnameff;
header.subjectInformation.patientWeight_kg = h.exam.patweight/1000;

if ~isempty(h.exam.dateofbirth)
    header.subjectInformation.patientBirthdate = [h.exam.dateofbirth(1:4) '-' h.exam.dateofbirth(5:6) '-' h.exam.dateofbirth(7:8)];
end

header.sequenceParameters.TR(1) = h.image.tr*1e-3;
header.sequenceParameters.TR(2) = h.image.tr*1e-3;
header.sequenceParameters.flipAngle_deg(1) = 0.5;
header.sequenceParameters.flipAngle_deg(2) = h.image.mr_flip;
header.sequenceParameters.TE = te90*1e3;

% The Encoding (Required)
header.encoding.trajectory = 'radial';
header.encoding.encodedSpace.fieldOfView_mm.x = 400;
header.encoding.encodedSpace.fieldOfView_mm.y = 400;
header.encoding.encodedSpace.fieldOfView_mm.z = 400;
header.encoding.encodedSpace.matrixSize.x = 64;
header.encoding.encodedSpace.matrixSize.y = 64;
header.encoding.encodedSpace.matrixSize.z = 64;
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
header.encoding.trajectoryDescription.identifier = "Custom Trajectory Info";
header.encoding.trajectoryDescription.userParameterDouble = [];
header.encoding.trajectoryDescription.userParameterLong(1).name = "ramp_time";
header.encoding.trajectoryDescription.userParameterLong(1).value = 72;

header.userParameters.userParameterLong(1).name = "xe_center_frequency";
header.userParameters.userParameterLong(1).value = h.ps.mps_freq/10;
header.userParameters.userParameterLong(2).name = "xe_dissolved_offset_frequency";
header.userParameters.userParameterLong(2).value = freq_off;
header.userParameters.userParameterString(1).name = "orientation";
header.userParameters.userParameterString(1).value = "coronal";


%% Serialize and write to the data set
xmlstring = LoadData.ismrmrd.xml.serialize(header);
dset.writexml(xmlstring);

dset.close();
