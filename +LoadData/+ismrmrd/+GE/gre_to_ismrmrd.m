function gre_to_ismrmrd(gre_file,mrdfile,MainInput)
% %% Handle Arguments
% if nargin == 0
%     [gre_file,mypath] = uigetfile('*.h5','Select GRE Raw Data file');
%     gre_file = fullfile(mypath,gre_file);
%     Subj_ID = inputdlg('Subject ID','Input Subject ID',[1 50]); % Or maybe better to prompt user to type in Subject ID?
%     Subj_ID = Subj_ID{1};
%     mrdfile = fullfile(mypath,[Subj_ID '_gre.h5']);
% elseif nargin == 1
%     [mypath,~,~] = fileparts(gre_file); 
%     Subj_ID = inputdlg('Subject ID','Input Subject ID',[1 50]); % Or maybe better to prompt user to type in Subject ID?
%     Subj_ID = Subj_ID{1};
%     mrdfile = fullfile(mypath,[Subj_ID '_gre.h5']);
% elseif nargin == 2
%     [mypath,~,~] = fileparts(gre_file); 
%     mrdfile = fullfile(mypath,[Subj_ID '_gre.h5']);
% end
%% 
if exist(mrdfile,'file')
    delete(mrdfile);
  %  error(['File ' mrdfile ' already exists.  Please remove first'])
end

dset = LoadData.ismrmrd.Dataset(mrdfile);

%% PJN - From here on, need GE-specific code

% ADH - Again, we are going to need this path for each site.  It won't change from scan to scan so hard to say whether
% we prompt for it every time or set it once for each site

% The path here should be the path to the local ventilation waveform file
directory = 'C:\XIPline\GE\waveforms\xe_calibration';
% Define file patterns to search for
filePatterns = {'*.mat', '*TR.fdl'};
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

wfn_file = filePaths{1};
% wfn = '/user/iibi/anhahn/code/mns-xe/xe_fgre/cart2D_129xe_fov384x360_mtx96x90_nexc90_kdt48_gmax39_smax147_dur5p22.mat';
wf = load(wfn_file);

[data,h] = LoadData.ismrmrd.GE.Functions.read_p(gre_file);

[nexc,~,nphases,nechoes,nslices,ncoils] = size(data);
nread = sum(wf.ind(1,:));

NAcq = nexc * nphases * nechoes * nslices * ncoils;
acqblock = LoadData.ismrmrd.Acquisition(NAcq);

acqblock.head.version(:) = 1;
acqblock.head.number_of_samples(:) = nread;
acqblock.head.active_channels(:) = 1;

for i = 1:NAcq
    acqblock.head.scan_counter(i) = i-1;

    Lin = mod(i-1,nexc);
    acqblock.head.idx.kspace_encode_step_1(i) = Lin;
    acqblock.head.idx.kspace_encode_step_2(i) = 0;
    acqblock.head.idx.average(i) = 0;
    acqblock.head.idx.phase(i) = 0;

    Sli = mod(floor((i-1)/(nexc*nphases*nechoes)), nslices);
    acqblock.head.idx.slice(i) = Sli;
    acqblock.head.idx.set(i) = 0;
    acqblock.head.idx.segment(i) = 0;
    acqblock.head.idx.repetition(i) = 0;
    acqblock.head.idx.contrast(i) = 0;

    % Set the flags
    acqblock.head.flagClearAll(i);
    if Lin == 0
        acqblock.head.flagSet('ACQ_FIRST_IN_SLICE', i);
    end
    if i == 1 
        acqblock.head.flagSet('ACQ_FIRST_IN_ENCODE_STEP1', i);
        acqblock.head.flagSet('ACQ_FIRST_IN_REPETITION', i);
    end
    if Lin == (nexc - 1)
        acqblock.head.flagSet('ACQ_LAST_IN_ENCODE_STEP1', i);
        acqblock.head.flagSet('ACQ_LAST_IN_SLICE', i);
    end
    if i == NAcq
        acqblock.head.flagSet('ACQ_LAST_IN_REPETITION', i);
    end

    Coil = ceil(NAcq/(nexc*nphases*nechoes*nslices));
    acqblock.data{i} = squeeze(data(Lin+1,wf.ind(Lin+1,:),1,1,Sli+1,Coil));
end

dset.appendAcquisition(acqblock);

header = [];

% Experimental Conditions (Required)
header.experimentalConditions.H1resonanceFrequency_Hz = 128000000; % 3T

% Acquisition System Information (Optional)
header.acquisitionSystemInformation.receiverChannels = 1;
header.acquisitionSystemInformation.systemFieldStrength_T = 3.0;
header.acquisitionSystemInformation.receiverChannels = 1;
header.acquisitionSystemInformation.systemVendor = 'GE';
header.acquisitionSystemInformation.systemModel = 'Premiere';   % Not sure if I can pull this out of the GE data header, will have to look into it more
header.acquisitionSystemInformation.institutionName = h.exam.ex_suid;

header.studyInformation.studyDate = ['20' h.rdb_hdr.scan_date(end-1:end) '-' h.rdb_hdr.scan_date(1:2) '-1'];


header.subjectInformation.patientID = Subj_ID;

% The Encoding (Required)
header.encoding.trajectory = 'Cartesian';
header.encoding.encodedSpace.fieldOfView_mm.x = h.rdb_hdr.fov * nexc / nread;
header.encoding.encodedSpace.fieldOfView_mm.y = h.rdb_hdr.fov;
header.encoding.encodedSpace.fieldOfView_mm.z = h.image.slthick;
header.encoding.encodedSpace.matrixSize.x = nexc;
header.encoding.encodedSpace.matrixSize.y = nread;
header.encoding.encodedSpace.matrixSize.z = 1;
% Recon Space
header.encoding.reconSpace.fieldOfView_mm.x = h.rdb_hdr.fov * nexc / nread;
header.encoding.reconSpace.fieldOfView_mm.y = h.rdb_hdr.fov;
header.encoding.reconSpace.fieldOfView_mm.z = h.image.slthick;
header.encoding.reconSpace.matrixSize.x = nexc;
header.encoding.reconSpace.matrixSize.y = nread;
header.encoding.reconSpace.matrixSize.z = 1;%NSlice;
% Encoding Limits
header.encoding.encodingLimits.kspace_encoding_step_0.minimum = 0;
header.encoding.encodingLimits.kspace_encoding_step_0.maximum = nexc-1;
header.encoding.encodingLimits.kspace_encoding_step_0.center = floor((nexc-1)/2);
header.encoding.encodingLimits.kspace_encoding_step_1.minimum = 0;
header.encoding.encodingLimits.kspace_encoding_step_1.maximum = nread-1;
header.encoding.encodingLimits.kspace_encoding_step_1.center = floor((nread-1)/2);
header.encoding.encodingLimits.slice.minimum = 0;
header.encoding.encodingLimits.slice.maximum = nslices-1;
header.encoding.encodingLimits.slice.center = floor((nslices-1)/2);

header.encoding.encodingLimits.contrast.minimum = 0;
header.encoding.encodingLimits.contrast.maximum = 0;
header.encoding.encodingLimits.contrast.center = 0;

header.encoding.encodingLimits.repetition.minimum = 0;
header.encoding.encodingLimits.repetition.maximum = 0;
header.encoding.encodingLimits.repetition.center = 0;

%% Serialize and write to the data set
xmlstring = LoadData.ismrmrd.xml.serialize(header);
dset.writexml(xmlstring);

dset.close();
