function diff_to_ismrmrd(gre_file,Subj_ID,mrdfile)
%% Handle Arguments
if nargin == 0
    [gre_file,mypath] = uigetfile('*.dat','Select Diffusion Raw Data file');
    gre_file = fullfile(mypath,gre_file);
    Subj_ID = inputdlg('Subject ID','Input Subject ID',[1 50]); % Or maybe better to prompt user to type in Subject ID?
    Subj_ID = Subj_ID{1};
    mrdfile = fullfile(mypath,[Subj_ID '_diff.h5']);
elseif nargin == 1
    [mypath,~,~] = fileparts(gre_file); 
    Subj_ID = inputdlg('Subject ID','Input Subject ID',[1 50]); % Or maybe better to prompt user to type in Subject ID?
    Subj_ID = Subj_ID{1};
    mrdfile = fullfile(mypath,[Subj_ID '_diff.h5']);
elseif nargin == 2
    [mypath,~,~] = fileparts(gre_file); 
    mrdfile = fullfile(mypath,[Subj_ID '_diff.h5']);
end
%% 
if exist(mrdfile,'file')
    delete(mrdfile);
  %  error(['File ' mrdfile ' already exists.  Please remove first'])
end

dset = ismrmrd.Dataset(mrdfile);

%Read in gre data
GRE_twix = DataImport.mapVBVD(gre_file,'ignoreSeg');
if length(GRE_twix) > 1
    GRE_twix = GRE_twix{end};
end
% Get various scan parameters
scanDate = GRE_twix.hdr.Phoenix.tReferenceImage0; 
scanDate = strsplit(scanDate,'.');
scanDate = scanDate{end};
scanDateStr = [scanDate(1:4),'-',scanDate(5:6),'-',scanDate(7:8)];
%scanDateStr = string(datetime(str2num(scanDate(1:4)),str2num(scanDate(5:6)),1));

%Read data
data = double(GRE_twix.image());

% data in format read (x coil) x phase encode x Rep? x slice 
%nX = size(data,1);
%nY = size(data,2);

acqblock = ismrmrd.Acquisition(GRE_twix.image.NAcq);

acqblock.head.version(:) = 1;
acqblock.head.number_of_samples(:) = GRE_twix.image.NCol;
acqblock.head.active_channels(:) = GRE_twix.image.NCha;

for i = 1:(GRE_twix.image.NAcq)
    acqblock.head.scan_counter(i) = i-1;
    acqblock.head.idx.kspace_encode_step_1(i) = GRE_twix.image.Lin(i)-1;
    acqblock.head.idx.kspace_encode_step_2(i) = GRE_twix.image.Par(i)-1;
    acqblock.head.idx.average(i) = GRE_twix.image.Ave(i)-1;
    acqblock.head.idx.phase(i) = GRE_twix.image.Phs(i)-1;
    acqblock.head.idx.slice(i) = GRE_twix.image.Sli(i)-1;
    acqblock.head.idx.set(i) = GRE_twix.image.Set(i)-1;
    acqblock.head.idx.segment(i) = GRE_twix.image.Seg(i)-1;
    acqblock.head.idx.repetition(i) = GRE_twix.image.Rep(i)-1;
    acqblock.head.idx.contrast(i) = GRE_twix.image.Eco(i)-1;
    if GRE_twix.image.Eco(i) == 2
        acqblock.head.user_float(1,i) = GRE_twix.hdr.MeasYaps.sWipMemBlock.adFree{12};
    else
        acqblock.head.user_float(1,i) = 0;
    end

    % Set the flags
    acqblock.head.flagClearAll(i);
    % if GRE_twix.image.Lin(i) == 1
    %     acqblock.head.flagSet('ACQ_FIRST_IN_SLICE', i);
    % end
    if i == 1 
        acqblock.head.flagSet('ACQ_FIRST_IN_ENCODE_STEP1', i);
        acqblock.head.flagSet('ACQ_FIRST_IN_REPETITION', i);
        acqblock.head.flagSet('ACQ_FIRST_IN_SLICE', i);
    elseif GRE_twix.image.Sli(i) > GRE_twix.image.Sli(i-1)
        acqblock.head.flagSet('ACQ_FIRST_IN_SLICE', i);
    end
    if i < GRE_twix.image.NAcq && GRE_twix.image.Sli(i) < GRE_twix.image.Sli(i+1)
        acqblock.head.flagSet('ACQ_LAST_IN_SLICE', i);
    elseif i == GRE_twix.image.NAcq
        acqblock.head.flagSet('ACQ_LAST_IN_ENCODE_STEP1', i);
    end
    % if GRE_twix.image.Lin(i) == GRE_twix.image.NLin
    %     acqblock.head.flagSet('ACQ_LAST_IN_ENCODE_STEP1', i);
    %     acqblock.head.flagSet('ACQ_LAST_IN_SLICE', i);
    % end
    if i == GRE_twix.image.NAcq
        acqblock.head.flagSet('ACQ_LAST_IN_REPETITION', i);
    end

    acqblock.data{i} = squeeze(data(:,:,GRE_twix.image.Lin(i),:,GRE_twix.image.Sli(i),:,:, GRE_twix.image.Eco(i)));
end

dset.appendAcquisition(acqblock);

header = [];

% Experimental Conditions (Required)
header.experimentalConditions.H1resonanceFrequency_Hz = 128000000; % 3T

%twix = DataImport.mapVBVD(ute_file);
% Acquisition System Information (Optional)
header.acquisitionSystemInformation.receiverChannels = size(data,2);
header.acquisitionSystemInformation.systemFieldStrength_T = GRE_twix.hdr.Dicom.flMagneticFieldStrength;
header.acquisitionSystemInformation.systemVendor = GRE_twix.hdr.Dicom.Manufacturer;
header.acquisitionSystemInformation.systemModel = GRE_twix.hdr.Dicom.ManufacturersModelName;
header.acquisitionSystemInformation.institutionName = GRE_twix.hdr.Dicom.InstitutionName;


header.studyInformation.studyDate = scanDateStr;
header.subjectInformation.patientID = Subj_ID;


NSlice = GRE_twix.hdr.Config.NSlc;
% The Encoding (Required)
header.encoding.trajectory = 'Cartesian';
header.encoding.encodedSpace.fieldOfView_mm.x = GRE_twix.hdr.MeasYaps.sSliceArray.asSlice{1,1}.dReadoutFOV;
header.encoding.encodedSpace.fieldOfView_mm.y = GRE_twix.hdr.MeasYaps.sSliceArray.asSlice{1,1}.dPhaseFOV;
header.encoding.encodedSpace.fieldOfView_mm.z = GRE_twix.hdr.MeasYaps.sSliceArray.asSlice{1,1}.dThickness;
header.encoding.encodedSpace.matrixSize.x = GRE_twix.image.NCol;
header.encoding.encodedSpace.matrixSize.y = GRE_twix.image.NLin;
header.encoding.encodedSpace.matrixSize.z = 1;%GRE_twix.image.NSli;
% Recon Space
header.encoding.reconSpace.fieldOfView_mm.x = GRE_twix.hdr.MeasYaps.sSliceArray.asSlice{1,1}.dReadoutFOV;
header.encoding.reconSpace.fieldOfView_mm.y = GRE_twix.hdr.MeasYaps.sSliceArray.asSlice{1,1}.dPhaseFOV;
header.encoding.reconSpace.fieldOfView_mm.z = GRE_twix.hdr.MeasYaps.sSliceArray.asSlice{1,1}.dThickness;
header.encoding.reconSpace.matrixSize.x = GRE_twix.hdr.Config.NImageCols;
header.encoding.reconSpace.matrixSize.y = GRE_twix.hdr.Config.NImageLins;
header.encoding.reconSpace.matrixSize.z = 1;%NSlice;
% Encoding Limits
header.encoding.encodingLimits.kspace_encoding_step_0.minimum = 0;
header.encoding.encodingLimits.kspace_encoding_step_0.maximum = GRE_twix.image.NCol-1;
header.encoding.encodingLimits.kspace_encoding_step_0.center = GRE_twix.image.centerCol(1)-1;
header.encoding.encodingLimits.kspace_encoding_step_1.minimum = 0;
header.encoding.encodingLimits.kspace_encoding_step_1.maximum = GRE_twix.image.NLin-1;
header.encoding.encodingLimits.kspace_encoding_step_1.center = GRE_twix.image.centerLin(1)-1;
header.encoding.encodingLimits.slice.minimum = 0;
header.encoding.encodingLimits.slice.maximum = GRE_twix.image.NSli-1;
header.encoding.encodingLimits.slice.center = floor((GRE_twix.image.NSli-1)/2);

header.encoding.encodingLimits.contrast.minimum = 0;
header.encoding.encodingLimits.contrast.maximum = GRE_twix.image.NEco-1;
header.encoding.encodingLimits.contrast.center = floor((GRE_twix.image.NEco-1)/2);

header.encoding.encodingLimits.repetition.minimum = 0;
header.encoding.encodingLimits.repetition.maximum = 0;
header.encoding.encodingLimits.repetition.center = 0;

%% Serialize and write to the data set
xmlstring = ismrmrd.xml.serialize(header);
dset.writexml(xmlstring);

dset.close();
