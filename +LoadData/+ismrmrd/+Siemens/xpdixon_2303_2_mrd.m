function xpdixon_2303_2_mrd(Xe_file,Subj_ID,mrdfile)
% Function to write default xpdixon sequence (Mugler) to MRD format for
% de-identified sharing of data

%% Read data
twix_obj = DataImport.mapVBVD(Xe_file);

%I really prefer using the twix object, but there are memory issues when
%doing it this way.
%xe_data = double(squeeze(twix_obj.image()));

%test = twix_obj.image(:,:,1,:,:,:,:,:,:,:,:);
xe_data = DataImport.ReadSiemensMeasVD13_idea(Xe_file);

raw = xe_data.rawdata;
loops = xe_data.loopcounters;

%% loop through radial acquisitions and store
%Preallocate memory - probably possible to find the number of gas echo
%times and number of dissolved echo times in the twix header, but hardcode
%for now
Gas = zeros(loops(1,16),max(loops(:,1)),2); %Hardcode
Dis = zeros(loops(1,16),max(loops(:,1)),3); %Hardcode
PostDis = zeros(loops(end,16),twix_obj.hdr.MeasYaps.sWipMemBlock.adFree{6});
PostGas = zeros(loops(end,16),twix_obj.hdr.MeasYaps.sWipMemBlock.adFree{11});

% appears to store both gas echoes, then both dissolved echoes
contrast = 1;
get_post = 1; %Need to keep track of what we've stored so we put data in the right place
for i = 1:size(loops,1)
    if loops(i,16) > 64
        if get_post == 1
            PostDis(:,loops(i,7)) = transpose(raw(i,:));
            if loops(i,7) == twix_obj.hdr.MeasYaps.sWipMemBlock.adFree{6};
                get_post = 0;
            end
        else
            PostGas(:,loops(i,7)) = transpose(raw(i,:));
        end
    else
        if contrast < 3 %Hardcode
            Gas(:,loops(i,1),loops(i,5)) = transpose(raw(i,1:loops(1,16)));
        elseif contrast > 2 %Hardcode
            Dis(:,loops(i,1),loops(i,5)) = transpose(raw(i,1:loops(1,16)));
        end
        contrast = contrast + 1;
        if contrast > 5 % Hardcode
            contrast = 1;
        end
    end
end
    
% Get Trajectories - Same trajectories for all contrasts
Traj = DataImport.gx_traj(size(Dis,1),size(Dis,2),twix_obj);

%% Create an empty ismrmrd dataset
if exist(mrdfile,'file')
    delete(mrdfile);
end

dset = ismrmrd.Dataset(mrdfile);

%% Prepare Imaging Acquisition Block and write to dataset
nX = size(Gas,1);
nY = size(Gas,2);
Dwell = twix_obj.hdr.MeasYaps.sRXSPEC.alDwellTime{1,1}*1e-9;
tsp = Dwell * 1e6; %need this in us %1e6/bw;
acqblock = ismrmrd.Acquisition((size(Gas,3) + size(Dis,3))*nY); 
acqblock.head.version(:) = 1;
acqblock.head.number_of_samples(:) = nX;
acqblock.head.sample_time_us(:) = tsp;
acqblock.head.active_channels(:) = 1;

Tot_Echoes = size(Gas,3) + size(Dis,3);

acqblock.head.trajectory_dimensions = 3*ones(1,Tot_Echoes*nY);
K = zeros(nX,Tot_Echoes*nY);
T = zeros(3,nX,Tot_Echoes*nY);

for i = 1:Tot_Echoes
    if i < size(Gas,3)+1
        K(:,i:Tot_Echoes:end) = Gas(:,:,i);
    else
        K(:,i:Tot_Echoes:end) = Dis(:,:,i-size(Gas,3));
    end
    T(:,:,i:Tot_Echoes:end) = Traj;
end

echo = 1; %Lots of good ways to keep track of what echo we're on... This isn't a particularly good one, but it'll work
for acqno = 1:Tot_Echoes*nY
    acqblock.head.scan_counter(acqno) = acqno-1;
    line = floor((acqno-1)/Tot_Echoes);
    acqblock.head.idx.kspace_encode_step_1(acqno) = line;
    acqblock.head.idx.repetition(acqno) = 0;
    %set contrast to differentiate dissolved/gas. 1 = gas. 2 = dissolved
    if echo < size(Gas,3) + 1 %mrd doesn't have a place for echo, so call it set instead
        acqblock.head.idx.set(acqno) = echo;
        acqblock.head.idx.contrast(acqno) = 1;
    else
        acqblock.head.idx.set(acqno) = echo - size(Gas,3);
        acqblock.head.idx.contrast(acqno) = 2;
    end
    echo = echo + 1;
    if echo > Tot_Echoes
        echo = 1;
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

%% Prepare Post Dissolved spectroscopy and add to acquisition block
acqblock = ismrmrd.Acquisition(size(PostDis,2));

acqblock.head.version(:) = 1;
acqblock.head.number_of_samples(:) = length(PostDis);
acqblock.head.sample_time_us(:) = twix_obj.hdr.MeasYaps.sWipMemBlock.adFree{10}*1e-6/2; % I think this should actually be twix_obj.hdr.MeasYaps.sWipMemBlock.adFree{10}, but not sure
acqblock.head.active_channels(:) = 1;
%acqblock.head.trajectory_dimensions = 3*ones(1,2*nY);


for i = 1:size(PostDis,2)
    acqblock.head.scan_counter(i) = i-1;
    acqblock.head.idx.kspace_encode_step_1(i) = i-1;
    acqblock.head.idx.repetition(i) = 0;
    acqblock.head.measurement_uid(i) = 1; %Identify the spectra easily
    %appended spectrum acquired on dissolved frequency, so set contrast to
    %2
    acqblock.head.idx.contrast(i) = 2;

    if i == 1
        acqblock.head.flagSet('ACQ_FIRST_IN_ENCODE_STEP1', 1);
        acqblock.head.flagSet('ACQ_FIRST_IN_SLICE', 1);
        acqblock.head.flagSet('ACQ_FIRST_IN_REPETITION', 1);
    end
    if i == size(PostDis,2)
        acqblock.head.flagSet('ACQ_LAST_IN_ENCODE_STEP1', 1);
        acqblock.head.flagSet('ACQ_LAST_IN_SLICE', 1);
        acqblock.head.flagSet('ACQ_LAST_IN_REPETITION', 1);
    end
    acqblock.data{i} = squeeze(PostDis(:,i));
    acqblock.traj{i} = [];
end
dset.appendAcquisition(acqblock);

%% Prepare Post Gas spectroscopy and add to acquisition block
acqblock = ismrmrd.Acquisition(size(PostGas,2));

acqblock.head.version(:) = 1;
acqblock.head.number_of_samples(:) = length(PostGas);
acqblock.head.sample_time_us(:) = twix_obj.hdr.MeasYaps.sWipMemBlock.adFree{15}*1e-6/2;%tsp; % I think this should actually be twix_obj.hdr.MeasYaps.sWipMemBlock.adFree{15}, but not sure
acqblock.head.active_channels(:) = 1;
%acqblock.head.trajectory_dimensions = 3*ones(1,2*nY);

for i = 1:size(PostGas,2)
    acqblock.head.measurement_uid(i) = 1; %Identify the spectra easily
    acqblock.head.scan_counter(i) = i-1;
    acqblock.head.idx.kspace_encode_step_1(i) = i-1;
    acqblock.head.idx.repetition(i) = 0;
    %appended spectrum acquired on gas frequency, so set contrast to
    %1
    acqblock.head.idx.contrast(i) = 1;

    if i == 1
        acqblock.head.flagSet('ACQ_FIRST_IN_ENCODE_STEP1', 1);
        acqblock.head.flagSet('ACQ_FIRST_IN_SLICE', 1);
        acqblock.head.flagSet('ACQ_FIRST_IN_REPETITION', 1);
    end
    if i == size(PostDis,2)
        acqblock.head.flagSet('ACQ_LAST_IN_ENCODE_STEP1', 1);
        acqblock.head.flagSet('ACQ_LAST_IN_SLICE', 1);
        acqblock.head.flagSet('ACQ_LAST_IN_REPETITION', 1);
    end
    acqblock.data{i} = squeeze(PostGas(:,i));
    acqblock.traj{i} = [];
end
dset.appendAcquisition(acqblock);

%% Create Header
header = [];

% Experimental Conditions (Required)
header.experimentalConditions.H1resonanceFrequency_Hz = 128000000; % 3T

% Acquisition System Information (Optional)
header.acquisitionSystemInformation.receiverChannels = 1;
header.acquisitionSystemInformation.systemFieldStrength_T = twix_obj.hdr.Dicom.flMagneticFieldStrength;
header.acquisitionSystemInformation.systemVendor = twix_obj.hdr.Dicom.Manufacturer;
header.acquisitionSystemInformation.systemModel = twix_obj.hdr.Dicom.ManufacturersModelName;
header.acquisitionSystemInformation.institutionName = twix_obj.hdr.Dicom.InstitutionName;

header.measurementInformation.patientPosition = '';
scanDate = twix_obj.hdr.Phoenix.tReferenceImage0; 
scanDate = strsplit(scanDate,'.');
scanDate = scanDate{end};
scanDateStr = [scanDate(1:4),'-',scanDate(5:6),'-',scanDate(7:8)];
header.studyInformation.studyDate = scanDateStr;
%header.studyInformation.studyDate = ['20' h.rdb_hdr.scan_date(end-1:end) '-' h.rdb_hdr.scan_date(1:2) '-' h.rdb_hdr.scan_date(4:5)];
header.subjectInformation.patientID = Subj_ID;


header.sequenceParameters.TR(1) = twix_obj.hdr.MeasYaps.alTR{2};
header.sequenceParameters.TR(2) = twix_obj.hdr.MeasYaps.alTR{3};
header.sequenceParameters.flipAngle_deg(1) = twix_obj.hdr.MeasYaps.adFlipAngleDegree{2};
header.sequenceParameters.flipAngle_deg(2) = twix_obj.hdr.MeasYaps.adFlipAngleDegree{3};
TE = zeros(1,Tot_Echoes);
for i = 1:size(Gas,3)
    TE(i) = twix_obj.hdr.MeasYaps.alTE{i}+4;
end
for i = 1:size(Dis,3)
    TE(i+size(Gas,3)) = twix_obj.hdr.MeasYaps.alTE{i+8}; %Dissolved TE's appear to start at index 9?
end
TE = TE/1000;

header.sequenceParameters.TE = TE;%te90 + 40e-6;

% The Encoding (Required) - Hardcode to standards
header.encoding.trajectory = 'radial';
header.encoding.encodedSpace.fieldOfView_mm.x = twix_obj.hdr.Config.ReadFoV;
header.encoding.encodedSpace.fieldOfView_mm.y = twix_obj.hdr.Config.ReadFoV;
header.encoding.encodedSpace.fieldOfView_mm.z = twix_obj.hdr.Config.ReadFoV;
header.encoding.encodedSpace.matrixSize.x = twix_obj.hdr.MeasYaps.sKSpace.lBaseResolution;
header.encoding.encodedSpace.matrixSize.y = twix_obj.hdr.MeasYaps.sKSpace.lBaseResolution;
header.encoding.encodedSpace.matrixSize.z = twix_obj.hdr.MeasYaps.sKSpace.lBaseResolution;
% Recon Space
% (in this case same as encoding space)
header.encoding.reconSpace = header.encoding.encodedSpace;
% Encoding Limits
header.encoding.encodingLimits.kspace_encoding_step_0.minimum = 0;
header.encoding.encodingLimits.kspace_encoding_step_0.maximum = size(K,2)-1;
header.encoding.encodingLimits.kspace_encoding_step_0.center = 0;
header.encoding.encodingLimits.kspace_encoding_step_1.minimum = 0;
header.encoding.encodingLimits.kspace_encoding_step_1.maximum = size(Gas,2)-1;
header.encoding.encodingLimits.kspace_encoding_step_1.center = floor(size(K,1)/2); %Center acq for a radial scan isn't really defined. This is fine for now.
header.encoding.encodingLimits.repetition.minimum = 0;
header.encoding.encodingLimits.repetition.maximum = 0;
header.encoding.encodingLimits.repetition.center = 0;

% Custom trajectory parameters

%up(1).name = "dwellTime";
%up(1).value = Dwell * 1e6; %Convert from s to us
a(1).name = "rampTime";
a(1).value = twix_obj.hdr.MeasYaps.sWipMemBlock.alFree{8};
up(1).name = "129Xe_center_frequency";
up(1).value = twix_obj.hdr.Dicom.lFrequency;
up(2).name = "129Xe_dissolved_offset_frequency";
up(2).value = twix_obj.hdr.MeasYaps.sWipMemBlock.alFree{5};

header.encoding.trajectoryDescription.identifier = "Custom Trajectory Info";
header.encoding.trajectoryDescription.userParameterDouble = a;
header.encoding.trajectoryDescription.userParameterLong = [];
header.userParameters.userParameterLong = up;

%% Serialize and write to the data set
xmlstring = ismrmrd.xml.serialize(header);
dset.writexml(xmlstring);

dset.close();

