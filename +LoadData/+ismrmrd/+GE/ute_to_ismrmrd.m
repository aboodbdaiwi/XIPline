function ute_to_ismrmrd(ute_file,Subj_ID,mrdfile)

% Create an empty ismrmrd dataset
if exist(mrdfile,'file')
    error(['File ' filename ' already exists.  Please remove first'])
end

% tmp = h5read('OUTPUT/MRI_Raw.h5','/Kdata/KData_E0_C0');
% K = complex(tmp.real,tmp.imag); clear tmp;
% kx = h5read('OUTPUT/MRI_Raw.h5','/Kdata/KX_E0');
% ky = h5read('OUTPUT/MRI_Raw.h5','/Kdata/KY_E0');
% kz = h5read('OUTPUT/MRI_Raw.h5','/Kdata/KZ_E0');

tmp = h5read(ute_file,'/Kdata/KData_E0_C0');
K = complex(tmp.real,tmp.imag); clear tmp;
kx = h5read(ute_file,'/Kdata/KX_E0');
ky = h5read(ute_file,'/Kdata/KY_E0');
kz = h5read(ute_file,'/Kdata/KZ_E0');

dset = LoadData.ismrmrd.Dataset(mrdfile);

nX = size(K,1);
nY = size(K,2);
acqblock = LoadData.ismrmrd.Acquisition(nY);

acqblock.head.version(:) = 1;
acqblock.head.number_of_samples(:) = nX;
acqblock.head.active_channels(:) = 1;
acqblock.head.trajectory_dimensions = 3*ones(1,nY);

T = zeros(3,nX,nY);

T(1,:,:) = kx;
T(2,:,:) = ky;
T(3,:,:) = kz;

for acqno = 1:nY
    acqblock.head.scan_counter(acqno) = acqno-1;
    acqblock.head.idx.kspace_encode_step_1(acqno) = acqno-1;
    acqblock.head.idx.repetition(acqno) = 0;
    
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
date_hold = datestr(datetime("today"),'yyyy-mm-dd');
date_hold((end-1):end) = '01';
header.studyInformation.studyDate = date_hold; %Not ideal - would rather get the date of scan from the raw data file if possible
header.subjectInformation.patientID = Subj_ID;


% The Encoding (Required)
header.encoding.trajectory = 'radial';
header.encoding.encodedSpace.fieldOfView_mm.x = 400;
header.encoding.encodedSpace.fieldOfView_mm.y = 400;
header.encoding.encodedSpace.fieldOfView_mm.z = 400;
header.encoding.encodedSpace.matrixSize.x = 128;
header.encoding.encodedSpace.matrixSize.y = 128;
header.encoding.encodedSpace.matrixSize.z = 128;
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

%% Serialize and write to the data set
xmlstring = LoadData.ismrmrd.xml.serialize(header);
dset.writexml(xmlstring);


dset.close();
