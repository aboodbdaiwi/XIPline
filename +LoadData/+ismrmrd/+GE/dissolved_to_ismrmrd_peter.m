function dissolved_to_ismrmrd(Xe_file,Subj_ID,mrdfile)
%%
%   mrdfile             -> Name of MRD file to write
%   Xe_file             -> cell containing three entries
%   Xe_file{1}          -> Scan Archive file generated from Dixon scan
%                         (e.g. ScanArchive*.h5)
%
%   Xe_file{2}          -> Mat file genereated from prep_dixon_gx
%                         (e.g. Dixon_ScanArchive*.mat)
%
%   Xe_file{3}          -> Mat file generated from recon_calibration
%                         (e.g ScanArchive*.mat, NOT Spect_ScanArchive.mat)

%% Create an empty ismrmrd dataset
if exist(mrdfile,'file')
    error(['File ' filename ' already exists.  Please remove first'])
end

dissolvedarchive = Xe_file{1};
dissolvedmat = Xe_file{2};
calmat = Xe_file{3};

[~,h] = read_p(dissolvedarchive);
load(dissolvedmat);
load(calmat,'te90');
load(calmat,'targetAX');

%tmp = read_fdl('/Users/ahahn/Work/Iowa/mns-xe/xe_dissolved/Genentech/radial3D_129Xe_fov400_mtx64_intlv2000_kdt20_gmax33_smax147_dur1p3_coca_rew1_G_freq.fdl');

% PJN - is there something we can do here to make this more generic and
% flexible?
tmp = read_fdl('/Users/ahahn/Work/Iowa/mns-xe/xe_dissolved/Genentech/ME_Sheffield/radial3D_129Xe_fov400_mtx64_intlv2000_kdt20_gmax30_smax118_dur4_coca_rew1_G_freq.fdl');
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
    %set contrast to differentiate dissolved/gas. 0 = gas. 1 = dissolved -
    %Does this point correctly to GE data?
    if mod(acqno,2) == 0
        acqblock.head.idx.contrast(acqno) = 1;
    else
        acqblock.head.idx.contrast(acqno) = 0;
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


header = [];

% Experimental Conditions (Required)
header.experimentalConditions.H1resonanceFrequency_Hz = 128000000; % 3T

% Acquisition System Information (Optional)
header.acquisitionSystemInformation.receiverChannels = 1;
header.acquisitionSystemInformation.systemFieldStrength_T = 3.0;
header.acquisitionSystemInformation.systemVendor = 'GE';
header.acquisitionSystemInformation.systemModel = 'Premiere';
header.acquisitionSystemInformation.institutionName = 'University of Iowa';

%header.measurementInformation.scandate = ['20' h.rdb_hdr.scan_date(end-1:end) '-' h.rdb_hdr.scan_date(1:2) '-' h.rdb_hdr.scan_date(4:5)];
header.measurementInformation.scandate = ['20' h.rdb_hdr.scan_date(end-1:end) '-' h.rdb_hdr.scan_date(1:2) '-1'];
header.measurementInformation.patientPosition = '';

%header.studyInformation.studyDate = ['20' h.rdb_hdr.scan_date(end-1:end) '-' h.rdb_hdr.scan_date(1:2) '-' h.rdb_hdr.scan_date(4:5)];
%header.subjectInformation.patientID = h.exam.patnameff(1:13);
header.subjectInformation.patientID = Subj_ID;

header.sequenceParameters.TR = 0.0075;
% For my purposes, default Dissolved to flip angle 1, Gas to Flip Angle 2
header.sequenceParameters.flipAngle_deg(1) = 20;
header.sequenceParameters.flipAngle_deg(2) = 0.5;
header.sequenceParameters.TE = te90 + 40e-6;

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

up(1).name = "dwellTime";
up(1).value = 20;
up(2).name = "rampTime";
up(2).value = 72;
up(3).name = "gasExciteFrequency";
up(3).value = targetAX;
up(4).name = "dissolvedExciteFrequency";
up(4).value = (targetAX + freq_off);

header.encoding.trajectoryDescription.identifier = "Custom Trajectory Info";
header.encoding.trajectoryDescription.userParameterDouble = up;
header.encoding.trajectoryDescription.userParameterLong = [];

%% Serialize and write to the data set
xmlstring = LoadData.ismrmrd.xml.serialize(header);
dset.writexml(xmlstring);


dset.close();
