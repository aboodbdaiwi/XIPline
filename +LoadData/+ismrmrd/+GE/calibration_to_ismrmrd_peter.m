function calibration_to_ismrmrd(calfile,Subj_ID,mrdfile)
%%
%   mrdfile -> Name of MRD file to write
%
%   calfile -> Scan Archive file generated from calibration scan
%               (e.g. ScanArchive*.h5)

%% Create an empty ismrmrd dataset
if exist(mrdfile,'file')
    error(['File ' filename ' already exists.  Please remove first'])
end

dset = LoadData.ismrmrd.Dataset(mrdfile);

[data,h] = read_p(calfile);

cf = h.ps.mps_freq/10;

% The path here should be the path to the local calibration waveform file
tmp = read_fdl('/Users/ahahn/Work/Iowa/mns-xe/xe_calibration/calibration3D_129xe_fov400_nsp256_intlv600_kdt40_dur10p2cal_freq.fdl');
freq_off = tmp(1); clear tmp;

nX = size(data,2);
nY = size(data,1);
tsp = 1e6/2/h.rdb_hdr.bw/1000;
acqblock = ismrmrd.Acquisition(nY);

acqblock.head.version(:) = 1;
acqblock.head.number_of_samples(:) = nX;
acqblock.head.sample_time_us(:) = tsp;
acqblock.head.active_channels(:) = 1;

% It would be really great to mark which spectra are dissolved and which
% are gas.
nDis = 500; %PJN - Don't like the hardcoding... Any way around this?

for acqno = 1:nY
    acqblock.head.scan_counter(acqno) = acqno-1;
    acqblock.head.idx.kspace_encode_step_1(acqno) = acqno-1;
    acqblock.head.idx.repetition(acqno) = 0;

    if acqno <= nDis % first scan is noise, then dissolved, so need to add 1 to the number to get this right
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
header.acquisitionSystemInformation.institutionName = 'University of Iowa';

%header.measurementInformation.scandate = ['20' h.rdb_hdr.scan_date(end-1:end) '-' h.rdb_hdr.scan_date(1:2) '-' h.rdb_hdr.scan_date(4:5)];
header.measurementInformation.scandate = ['20' h.rdb_hdr.scan_date(end-1:end) '-' h.rdb_hdr.scan_date(1:2) '-1'];

header.measurementInformation.patientPosition = '';
%header.studyInformation.studyDate = ['20' h.rdb_hdr.scan_date(end-1:end) '-' h.rdb_hdr.scan_date(1:2) '-' h.rdb_hdr.scan_date(4:5)];
%header.subjectInformation.patientID = h.exam.patnameff(1:13);
header.subjectInformation.patientID = Subj_ID;

header.sequenceParameters.TR = h.image.tr*1e-6;
header.sequenceParameters.flipAngle_deg(1) = 10;
header.sequenceParameters.flipAngle_deg(2) = 20;
%header.sequenceParameters.flipAngle_deg(2) = h.image.mr_flip;
header.sequenceParameters.TE = h.rdb_hdr.te*1e-6;

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

up(1).name = "dwellTime";
up(1).value = 40;
up(2).name = "gasExciteFrequency";
up(2).value = cf;
up(3).name = "dissolvedExciteFrequency";
up(3).value = (cf + freq_off);

header.encoding.trajectoryDescription.identifier = "Custom Trajectory Info";
header.encoding.trajectoryDescription.userParameterDouble = up;
header.encoding.trajectoryDescription.userParameterLong = [];


%% Serialize and write to the data set
xmlstring = ismrmrd.xml.serialize(header);
dset.writexml(xmlstring);


dset.close();
