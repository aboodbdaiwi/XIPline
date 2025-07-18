function success = GasExchange_ProcessSingleSubject(DataLocation,NewMask,NewImages,NewReg)
%GasExchange_ProcessSingleSubject - Reconstruct and analyse gas exchange MRI from the CPIR Complete Patch
%   Will process the resulting gas exchange mri data from the CPIR complete
%    patch and export all data and analysis to same folder - for V3 and 
%    XeCTC versions of the data
%   Requires a folder with the Dissolved_Xe raw/lab/sin files,
%    Dissolved_Xe_Mask raw/lab/sin files, and two (and only two!) sets of 
%    data/list files for the dissolved and proton images
%   If no folder is provided, user will be prompted to select the folder
%
%   Syntax:  GasExchange_ProcessSingleSubject(DataLocation)
%
%   Inputs:
%      DataLocation - (optional) string containing the path of the folder
%      NewMask - (optional) 1 forces new masks to be made, 0 imports masks if present
%      NewImages - (optional) 1 forces new images to be made, 0 imports images if present
%      NewReg - (optional) 1 forces registration to run again, 0 skips it
%
%   Outputs:
%      success - 1 if ran successfully, 0 if not 
%
%   Example: 
%      GasExchange_ProcessSingleSubject()
%      GasExchange_ProcessSingleSubject('C:\Users\mwillmering\Documents\Subject1')
%      GasExchange_ProcessSingleSubject('C:\Users\mwillmering\Documents\Subject1',0,0,0)
%
%   See also: GasExchange_ProcessBatch
% 
%   Package: https://github.com/cchmc-cpir/CCHMC-Gas-Exchange-Processing-Package
%
%   Author: Matthew Willmering
%   Work email: matthew.willmering@cchmc.org
%   Personal email: matt.willmering@gmail.com
%   Website: https://cpir.cchmc.org/

clearvars -except DataLocation NewMask NewImages NewReg; close all;
opengl hardware
success = 1;
disp('Gas Exchange Processing Beginning...')
subject_tic = tic;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Import Data, variables, Trajectories                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Define Fixed Variables, Subject to Change with Sequence
%General information
try
    ReconVersion = GasExchange.get_git_hashobject;%use git to find hash of commit used
    ReconVersion = ReconVersion(1:40);%remove enter at end
catch
    ReconVersion = 'TEST';%If not connected to git, can't determine hash so state test
end
FunctionDirectory = which('GasExchange_ProcessSingleSubject');
idcs = strfind(FunctionDirectory,filesep);%determine location of file separators
FunctionDirectory = FunctionDirectory(1:idcs(end)-1);%remove file
NumPlotSlices = 7;%must be odd to work as expected


%% Determine Folder Location
%Determine folder to use
disp('Looking For Necessary Files...')
if exist('DataLocation','var') == 0
  DataLocation = uigetdir('C:\', 'Select Folder containing Gas Exchange Data');
end
if exist('NewMask','var') == 0
  NewMask = 1; %1 forces new masks to be made, 0 imports masks if present
end
if exist('NewImages','var') == 0
  NewImages = 1; %1 forces new images to be made, 0 imports images if present
end
if exist('NewReg','var') == 0
  NewReg = 1; %1 forces registration to run
end

%% Check Inputs are Valid
if(and(~isfile([DataLocation,'\GasExchangeWorkspace.mat']),(NewImages==0 || NewReg==0)))
    disp('Error. Existing Workspace Not Found...')
    NewImages = 1; %force new images
    NewReg = 1; %force new registration
    disp('Will Create New Images and Registration.')
end

if(and(~isfile([DataLocation,'\ProtonMask.nii.gz']),NewMask==0))
    disp('Error. Existing Mask Not Found...')
    NewMask = 1; %force new mask
    disp('Will Create New Mask.')
end

%check folder is valid
if DataLocation == 0 %exit if no folder selected
    disp('Error. Must select a folder')
    success = 0;
    return;
end
if isfolder(char(DataLocation)) == 0  %exit if input is not a folder
   disp('Error. Input must be a folder')
   success = 0;
   return;
end


%% Import Previous Workspace if skipping Images/Mask
if(NewImages == 0 || NewReg == 0) %if using values from previous run, import that workspace
    disp('Importing Previous Workspace for images and/or mask...')
    %Dont want to overwrite already declared varaibles
    varlist = who; %Find the variables that already exist
    varlist =strjoin(varlist','$|'); %Join into string, separating vars by '|'
    load([DataLocation,'\GasExchangeWorkspace.mat'],'-regexp',['^(?!' varlist ')\w']);%#ok<*LOAD> %import all except already declared
    clear varlist;
    disp('Importing Previous Workspace for images and/or mask Completed.')
end


%% Determine files to use
HSinFile = dir([DataLocation,'\*Dissolved_Xe_Mask_20191008.sin']);
XeSinFile = dir([DataLocation,'\*Dissolved_Xe_20191008.sin']);
ScanVersion = 'CCHMC_V3';
if (size(HSinFile,1) < 1 && size(XeSinFile,1) < 1) %didn't find that sequence naming version; check for new sequence names
    HSinFile = dir([DataLocation,'\*CPIR_Mask_Gas_Exchange.sin']);
    XeSinFile = dir([DataLocation,'\*CPIR_Gas_Exchange.sin']);
    ScanVersion = 'CCHMC_V3';
end
if (size(HSinFile,1) < 1 && size(XeSinFile,1) < 1) %didn't find that sequence version
    HSinFile = dir([DataLocation,'\*CTC_Anatomical.sin']);
    XeSinFile = dir([DataLocation,'\*CTC_GasExchange_Xe.sin']);
    ScanVersion = 'Xe_CTC';
end
if (size(HSinFile,1) < 1 && size(XeSinFile,1) < 1) %didn't find that sequence naming version; check for new sequence names
    HSinFile = dir([DataLocation,'\*XeCTC_Mask_Gas_Exchange.sin']);
    XeSinFile = dir([DataLocation,'\*XeCTC_Gas_Exchange.sin']);
    ScanVersion = 'Xe_CTC';
end
if (size(HSinFile,1) < 1 && size(XeSinFile,1) == 1) %didn't find that sequence naming version; check for CTC xe and duke mask(which is same as CTC)
    HSinFile = dir([DataLocation,'\*DukeIPF_Mask_Gas_Exchange.sin']);
    XeSinFile = dir([DataLocation,'\*XeCTC_Gas_Exchange.sin']);
    ScanVersion = 'Xe_CTC';
end
if (size(HSinFile,1) < 1 && size(XeSinFile,1) < 1) %didn't find that sequence naming version; check for duke version/Intermediate dixon
    HSinFile = dir([DataLocation,'\*DukeIPF_Mask_Gas_Exchange.sin']);
    XeSinFile = dir([DataLocation,'\*DukeIPF_Gas_Exchange.sin']);
    ScanVersion = 'DukeIPF';
end

DataFiles = dir([DataLocation,'\*.data']);
if size(DataFiles,1) ~= 2
    disp('Error. Does not have exactly 2 data files')
    success = 0;
    return;
end
if DataFiles(1).bytes < DataFiles(2).bytes
    XeDataFile = DataFiles(1);
    HDataFile = DataFiles(2);
    clear DataFiles
else
    XeDataFile = DataFiles(2);
    HDataFile = DataFiles(1);
    clear DataFiles
end
if size(HSinFile,1) ~=1
    disp('Error. Could not determine Dissolved Proton Sin file')
    success = 0;
    return;
elseif size(XeSinFile,1) ~=1
    disp('Error. Could not determine Dissolved Xe Sin file')
    success = 0;
    return;
elseif size(XeDataFile,1) ~=1
    disp('Error. Could not determine Dissolved Xe Data file')
    success = 0;
    return;
elseif size(HDataFile,1) ~=1
    disp('Error. Could not determine Dissolved Proton Data file')
    success = 0;
    return;
end
disp('Looking For Necessary Files Completed.')


%% Import Acquisition Information
disp('Importing Acquisition Information...')
HSin = GasExchange.loadSIN([HSinFile.folder,'\',HSinFile.name]);
XeSin = GasExchange.loadSIN([XeSinFile.folder,'\',XeSinFile.name]);

extraOvs = false;
OvsFactor = 1.0;
if strcmp(ScanVersion,'CCHMC_V3')
    dwell_ms = 1.1049/58; %in ms - for 58 points with 1.1049ms readout   
    if (exist('XeSin','var') && isfield(XeSin, 'sample_time_interval')) %R59 passes it so will update here
        dwell_ms = XeSin.sample_time_interval.vals(1) / 1000;
    end
    if (exist('XeSin','var') && XeSin.max_encoding_numbers.vals(1) ~= 57) %R59 oversampling; RC thinks its cartesian
        extraOvs = true;
        OvsFactor = (XeSin.max_encoding_numbers.vals(1)+1)/58;
    end
    pw = 0.65; %pulse width in ms
    DisFlipAngle = 20; %Dissolved flip angle
    GasFlipAngle = 0.5; %Gas flip angle
    freq_jump = 7143; %change in freqeuncy for dissolved and off-res
elseif strcmp(ScanVersion,'Xe_CTC')
    dwell_ms = 0.64/68; %in ms - for 64 points with 0.64ms readout   
    if (exist('XeSin','var') && isfield(XeSin, 'sample_time_interval')) %R59 passes it so will update here
        dwell_ms = XeSin.sample_time_interval.vals(1) / 1000;
    end
    if (exist('XeSin','var') && XeSin.non_cart_max_encoding_nrs.vals(1) ~= 63) %R59 oversampling
        extraOvs = true;
        OvsFactor = (XeSin.non_cart_max_encoding_nrs.vals(1)+1)/64;
    end
    pw = 0.65; %pulse width in ms
    DisFlipAngle = 20; %Dissolved flip angle
    GasFlipAngle = 0.5; %Gas flip angle
    freq_jump = 7702; %change in freqeuncy for dissolved and off-res
elseif strcmp(ScanVersion,'DukeIPF')
    dwell_ms = 0.64/68; %in ms - for 64 points with 0.64ms readout   
    if (exist('XeSin','var') && isfield(XeSin, 'sample_time_interval')) %R59 passes it so will update here
        dwell_ms = XeSin.sample_time_interval.vals(1) / 1000;
    end
    if (exist('XeSin','var') && XeSin.non_cart_max_encoding_nrs.vals(1) ~= 63) %R59 oversampling
        extraOvs = true;
        OvsFactor = (XeSin.non_cart_max_encoding_nrs.vals(1)+1)/64;
    end
    pw = 0.67; %pulse width in ms
    DisFlipAngle = 15; %Dissolved flip angle
    GasFlipAngle = 0.5; %Gas flip angle
    freq_jump = 7340; %change in freqeuncy for dissolved and off-res
end
dwell_s = dwell_ms / 1000; %convert to s


[~, folderName, ~] = fileparts(DataLocation);
folderName = strsplit(folderName);
Subject = folderName{1,2:end};
Notes = '';
if(size(folderName,2)>2)
    Notes = cat(2,folderName{1,3:end});
end
ScanDate = [folderName{1,1}(1:4),'-',folderName{1,1}(5:6),'-',folderName{1,1}(7:8)];
if strcmp(ScanVersion,'Xe_CTC') || strcmp(ScanVersion,'DukeIPF')
    XeTR = XeSin.repetition_times.vals(1)*2/1000;%account for 2 acqs, convert to s
else
    XeTR = XeSin.repetition_times.vals(1)*3/1000;%account for 3 acqs, convert to s
end
TE90 = XeSin.echo_times.vals(1);
if strcmp(ScanVersion,'Xe_CTC') || strcmp(ScanVersion,'DukeIPF')
    %CTC protocl always defines dissolved TE as mid pulse to acq
    ActTE90 = TE90; %TE90 from center of pulse
    Scanner = 'CheckRecords';
elseif(TE90>0.5 * pw)
    %R5.6.0 dMN defines dissolved TE as mid pulse to acq already
    ActTE90 = TE90; %TE90 from center of pulse
    Scanner = '3T-T1';
else
    %R5.3.1 defines dissolved TE as end pulse to acq
    ActTE90 = 0.5 * pw + TE90; %TE90 from center of pulse
    Scanner = '3T-R';
end
Xe_AcqMatrix = XeSin.scan_resolutions.vals(1);
H_AcqMatrix = HSin.scan_resolutions.vals(1);
if strcmp(ScanVersion,'Xe_CTC') || strcmp(ScanVersion,'DukeIPF') %recon 2x interpolated
    Xe_RecMatrix = 2*Xe_AcqMatrix;
    H_RecMatrix = 2*H_AcqMatrix;
else %CCHMC V3, recon 2x interpolated
    Xe_RecMatrix = 2*Xe_AcqMatrix;
    H_RecMatrix = H_AcqMatrix; %already acquired at 2x resolution
end
disp('Importing Acquisition Information Completed.')


%% Import Data
disp('Importing Data...')
%Load in FIDs
[XeData,XeInfo] = GasExchange.loadLISTDATA([XeDataFile.folder,'\',XeDataFile.name]);
XeData = squeeze(XeData);
[HData,HInfo] = GasExchange.loadLISTDATA([HDataFile.folder,'\',HDataFile.name]);
HData = squeeze(HData);

if strcmp(ScanVersion,'Xe_CTC') || strcmp(ScanVersion,'DukeIPF')
    %Dissolved k-space
    DissolvedKSpace = XeData(:,:,2);
    %Gas k-space
    GasKSpace = XeData(:,:,1);
    %Proton k-space
    HKSpace = HData;

    if(extraOvs)
        DissolvedKSpace = movmean(DissolvedKSpace,OvsFactor);
        DissolvedKSpace = downsample(DissolvedKSpace,OvsFactor);
        
        GasKSpace = movmean(GasKSpace,OvsFactor);
        GasKSpace = downsample(GasKSpace,OvsFactor);

        dwell_s = dwell_s * OvsFactor;
    end

    %Get sizes of dimensions
    Xe_nprof = size(DissolvedKSpace,2);
    Xe_interleaves = 1;
    H_nsamp = size(HKSpace,1);
    H_nprof = size(HKSpace,2);
    H_coils = size(HKSpace,3);
    H_interleaves = 1;
else
    %Split Xe data into components; formatted as [read, proj, interleave, dyn, mix]
    %Dissolved Spectra
    PreDissolvedFID = XeData(:,1,1,1,2);
    PostDissolvedFID = XeData(:,2,1,1,2);

    %%Gas Spectra - Not Needed
    %%PreGasFID = XeData(:,1,1,2,2);
    PostGasFID = XeData(:,2,1,2,2);

    %%Off-Res Spectra - Not needed
    %%PreOffResFID = XeData(:,1,1,3,2);
    %%PostOffResFID = XeData(:,2,1,3,2);

    %Dissolved k-space
    DissolvedKSpaceInit = XeData(:,:,:,1,1);
    %Gas k-space
    GasKSpaceInit = XeData(:,:,:,2,1);
    %%Off-Res k-space - Not needed
    %%OffResKSpace = XeData(:,:,:,3,1);
    %Proton k-space
    HKSpaceInit = HData;

    if(extraOvs)
        PreDissolvedFID = movmean(PreDissolvedFID,round(OvsFactor));
        PreDissolvedFID =  downsample(PreDissolvedFID,round(OvsFactor));

        PostDissolvedFID = movmean(PostDissolvedFID,round(OvsFactor));
        PostDissolvedFID =  downsample(PostDissolvedFID,round(OvsFactor));

        PostGasFID = movmean(PostGasFID,round(OvsFactor));
        PostGasFID =  downsample(PostGasFID,round(OvsFactor));

        DissolvedKSpaceInit = movmean(DissolvedKSpaceInit,round(OvsFactor));
        DissolvedKSpaceInit = downsample(DissolvedKSpaceInit,round(OvsFactor));
        
        GasKSpaceInit = movmean(GasKSpaceInit,round(OvsFactor));
        GasKSpaceInit = downsample(GasKSpaceInit,round(OvsFactor));

        dwell_s = dwell_s * round(OvsFactor);
    end

    %Correct Atten. 
    %Since R59, different levels of attn are used for the two mixes, scale here using last gas kspace and first post gas spec)
    scaleFac = abs(PostGasFID(1,1))/abs(GasKSpaceInit(1,end,end)); %These should be almost identical at k0
    DissolvedKSpaceInit = DissolvedKSpaceInit * scaleFac;
    GasKSpaceInit = GasKSpaceInit * scaleFac;

    %Get sizes of dimensions
    XeSpec_nsamp = size(DissolvedKSpaceInit,1);
    XeImg_nsamp = (XeSin.max_encoding_numbers.vals(1)+1) / OvsFactor;
    Xe_nprof = size(DissolvedKSpaceInit,2);
    Xe_interleaves = size(DissolvedKSpaceInit,3);
    H_nsamp = size(HKSpaceInit,1);
    H_nprof = size(HKSpaceInit,2);
    H_interleaves = size(HKSpaceInit,3);
    H_coils = size(HKSpaceInit,4);

    %Trim Kspaces to correct number of readout points
    DissolvedKSpaceInit = DissolvedKSpaceInit(1:XeImg_nsamp,:,:);
    GasKSpaceInit = GasKSpaceInit(1:XeImg_nsamp,:,:);

    % Put KSpaces in order and logical format
    XeOrder = GasExchange.GetProjectionOrder(XeInfo);
    HOrder = GasExchange.GetProjectionOrder(HInfo);
    DissolvedKSpace = GasExchange.SortUTEData_AcqOrder(GasExchange.reshapeUTEData(DissolvedKSpaceInit),[XeImg_nsamp, Xe_nprof, Xe_interleaves], XeOrder);
    GasKSpace = GasExchange.SortUTEData_AcqOrder(GasExchange.reshapeUTEData(GasKSpaceInit),[XeImg_nsamp, Xe_nprof, Xe_interleaves], XeOrder);
    HKSpace = zeros(H_nsamp,H_nprof*H_interleaves,H_coils);
    for i=1:H_coils
        HKSpace(:,:,i) = GasExchange.SortUTEData_AcqOrder(GasExchange.reshapeUTEData(HKSpaceInit(:,:,:,i)),[H_nsamp, H_nprof, H_interleaves], HOrder);
    end
end
disp('Importing Data Completed.')


%% Calculate Trajectories
disp('Calculating Trajectories...')
if strcmp(ScanVersion,'Xe_CTC') || strcmp(ScanVersion,'DukeIPF')
    del = 1.25;
    if (extraOvs)
        del = del * OvsFactor;
    end
    XeTraj = GasExchange.philipsradialcoords(del,2,[XeSinFile.folder,'\',XeSinFile.name]); %1.25us delay, Haltoned Spiral
    XeTraj = permute(XeTraj,[3 2 1 4]); %ro, proj, intlv, dims
    if (extraOvs)
        XeTraj = movmean(XeTraj,OvsFactor);
        XeTraj =  downsample(XeTraj,OvsFactor);
    end
    XeTraj = permute(XeTraj,[4 1 2 3]); %dims, ro, proj, intlv
    del = 1.25;
    HTraj = GasExchange.philipsradialcoords(del,2,[HSinFile.folder,'\',HSinFile.name]); %1.25us delay, Haltoned Spiral
    HTraj = permute(HTraj,[4 3 2 1]);
else
    if (strcmp(Scanner,'3T-R'))%R5.3.1
        XeTraj = GasExchange.philipsradialcoords(0.36,1,[FunctionDirectory,'\V3 Scan Info\20191008_162511_Dissolved_Xe_20191008 - 3T-R.sin']); %0.36us delay, GM, from non-spectroscopy version
    else%3T-T1; R5.6.0 dMN
        del = 0.36;
        if mod(OvsFactor,1) % strange situation for some cases that were based on 57 points readout
            del = -1;
        end
        XeTraj = GasExchange.philipsradialcoords(del,1,[FunctionDirectory,'\V3 Scan Info\20200210_133229_Dissolved_Xe_20191008 - 3T-T1.sin']); %0.36us delay, GM, from non-spectroscopy version
    end
    del = 1.65;
    if mod(OvsFactor,1) % strange situation for some cases that were based on 57 points readout
        del = -1.65;
    end
    XeTraj = permute(XeTraj,[4 3 2 1]);
    XeTraj = GasExchange.SortUTETraj_AcqOrder(GasExchange.reshapeUTETraj(XeTraj), [XeImg_nsamp, Xe_nprof, Xe_interleaves], XeOrder);%reshape and order
    HTraj = GasExchange.philipsradialcoords(del,1,[HSinFile.folder,'\',HSinFile.name]); %1.65us delay, GM
    HTraj = permute(HTraj,[4 3 2 1]);
    HTraj = GasExchange.SortUTETraj_AcqOrder(GasExchange.reshapeUTETraj(HTraj), [H_nsamp, H_nprof, H_interleaves], HOrder);%reshape and order
end
disp('Calculating Trajectories Completed.')


%% Account for Exhale/Motion in Proton Image
%Functions to fit
for coil = 1:size(HKSpace,3)
    HPhase(:,coil) = smooth(angle(HKSpace(1,:,coil)),15,'rloess');
end
HPhaseDynamics = rad2deg(mean(unwrap(HPhase,1),2));
HPulses = (1:length(HPhaseDynamics))';

%Fits to data; use only 1st third to avoid exhale/movement - assuming first 5 seconds are uncorrupted
HDynMean = mean(HPhaseDynamics(1000:round(length(HPulses)/3)));
HDynFit = HDynMean*ones(size(HPhaseDynamics(1000:round(length(HPulses)/3))));
HDynRMSE = sqrt(immse(HDynFit,HPhaseDynamics(1000:round(length(HPulses)/3))));
if (3*HDynRMSE<5) 
    HDynRMSE = 5/3; %If small phase varaition, use acceptance window of 5deg
end

%smooth to remove small excursions
HPhaseDynamics_Smoothed = medfilt1(HPhaseDynamics,15,'truncate');

%Determine Proton Data Cutoffs
HDynFit = HDynMean*ones(size(HPulses));%calculate function at all x points
keep = abs(HDynFit - HPhaseDynamics_Smoothed) < 3*HDynRMSE;%label as outlier if >3 rmse away
keep = medfilt1(single(keep),15,'truncate');%smooth to remove small excursions
keep(1:round(length(HPulses)/3)) = 1; %keep all from fit
motion_ind = find(keep==0,1,'first');%find first non keep index

%H Motion
HKSpace_Steady = HKSpace;
HTraj_Steady = HTraj;
if ~isempty(motion_ind)%only remove data if motion
    HKSpace_Steady(:,motion_ind:end,:) = [];
    HTraj_Steady(:,:,motion_ind:end) = [];
else
    motion_ind = length(HPulses);
end

%% Slow Switching Fix
[~,GasKSpace_MaxInd] = max(abs(GasKSpace));
if mode(GasKSpace_MaxInd,2) ~= 1
    Notes = 'SlowSwitchingFix'; % add SlowSwitchingFix tag to database.
    initPt = mode(GasKSpace_MaxInd,2);
    disp([newline 'IMPORTANT INFO...']);
    pause(2)
    disp('...');pause(1); disp('...');pause(1); disp('Max(GasKSpace(1,:) not at t=0 of acquisition.');
    pause(2)
    disp('...'); pause(1); disp('...');pause(1); disp('Applying SlowSwitchingFix...');
    pause(1)
    for p=1:initPt
        GasKSpace(p,:) = GasKSpace(initPt,:); 
        DissolvedKSpace(p,:) = DissolvedKSpace(initPt,:); 
        PostDissolvedFID(p) = PostDissolvedFID(initPt);
        PostGasFID(p) = PostGasFID(initPt);
        PreDissolvedFID(p) = PreDissolvedFID(initPt);
    end
    fileID_info = fopen(fullfile(DataLocation,'info.txt'),'w');
    fprintf(fileID_info,'%s    %s\n','SlowSwitchingFix','yes');
    fprintf(fileID_info,'%s    %s\n','initPt',num2str(initPt));
    if mode(GasKSpace_MaxInd,2) > 3
        fprintf(fileID_info,'%s    %s\n','GasContaminationRemoval','off');
    else
        fprintf(fileID_info,'%s    %s\n','GasContaminationRemoval','on');
    end
    fclose(fileID_info);
end

%% Import Healthy Distribution, If Available
%Healthy Cohort
ParentFolder = fullfile(DataLocation, '..');
if(exist([ParentFolder,'\HealthyCohort.mat'],'file') == 2) %if Healthy cohort data exists, use
    HealthyFile = dir([ParentFolder,'\HealthyCohort.mat']);%get file info
    HealthyData = load([HealthyFile.folder,'\',HealthyFile.name], '-regexp','^(?!Comb|Ind)...');%import all variable but figures
    VentThresh = HealthyData.thresholds.Vent;
    DissolvedThresh = HealthyData.thresholds.Dissolved;
    BarrierThresh = HealthyData.thresholds.Barrier;
    RBCThresh = HealthyData.thresholds.RBC;
    RBCBarrThresh = HealthyData.thresholds.RBCBarr;
    RBCOscThresh = HealthyData.thresholds.RBCOsc;
    HealthyCohortNum = HealthyData.CohortTag; %Healthy Cohort Data
    HealthyDistPresent = 1;
else %use thresholds from Duke 1.5T Quantitative Analysis Paper 10.1002/mp.12264,
    VentThresh = [0.51-2*0.19,0.51-1*0.19,0.51,0.51+1*0.19,0.51+2*0.19];%6
    DissolvedThresh = [0.0075-2*0.00125,0.0075-1*0.00125,0.0075,0.0075+1*0.00125,0.0075+2*0.00125];%6
    BarrierThresh = [0.0049-2*0.0015,0.0049-1*0.0015,0.0049,0.0049+1*0.0015,0.0049+2*0.0015,0.0049+3*0.0015,0.0049+4*0.0015];%8
    RBCThresh = [0.0026-2*0.0010,0.0026-1*0.0010,0.0026,0.0026+1*0.0010,0.0026+2*0.0010];%6
    RBCBarrThresh = [0.53-2*0.18,0.53-1*0.18,0.53,0.53+1*0.18,0.53+2*0.18];%6
    RBCOscThresh = [8.9596-2*10.5608,8.9596-1*10.5608,8.9596,8.9596+1*10.5608,8.9596+2*10.5608,8.9596+3*10.5608,8.9596+4*10.5608];%8
    HealthyCohortNum = 'INITIAL  GUESS'; %Healthy Cohort Data
    HealthyDistPresent = 0;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Transform Data into Images and Spectra                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Determine FOV Shift
if(NewImages == 1)
    %not using planned shifts as techs don't seem to do it consistently
    disp('Determining Image offset...')
    PixelShiftinit = [0; 0; 0]; %start with no offset
    converged = 0; %start not converged
    iter = 0;
    while (converged == 0) %run until no change in shift
        tempIm = GasExchange.Dissolved_LowResImageRecon(Xe_RecMatrix,GasKSpace,XeTraj/2,PixelShiftinit); %make image; 2x Resolution
        [tempIm, ~] = GasExchange.Vent_BiasCorrection_Coarse(tempIm);
        tempLab = bwlabeln(imbinarize(abs(tempIm)/(0.5*max(abs(tempIm(:)))))); %make mask label
        center = regionprops3(tempLab); %get centroids
        center = sortrows(center,1,'descend'); %order by size
        if (size(center,1) > 1 && table2array(center(2,1)) > 0.15 * table2array(center(1,1)))        
            CentVox(1) = (min(center.BoundingBox(:,1))+max(center.BoundingBox(:,1)+center.BoundingBox(:,4)))/2;
            CentVox(2) = (min(center.BoundingBox(:,2))+max(center.BoundingBox(:,2)+center.BoundingBox(:,5)))/2;
            CentVox(3) = (min(center.BoundingBox(:,3))+max(center.BoundingBox(:,3)+center.BoundingBox(:,6)))/2;
        else
            CentVox(1) = center.BoundingBox(1,1)+center.BoundingBox(1,4)/2;
            CentVox(2) = center.BoundingBox(1,2)+center.BoundingBox(1,5)/2;
            CentVox(3) = center.BoundingBox(1,3)+center.BoundingBox(1,6)/2;
        end
        PixelShift = 0.5*Xe_RecMatrix - CentVox'; %get pixel shift
        PixelShift = PixelShift([3, 2, 1]); %put in same order
        PixelShift(2:3) = -PixelShift(2:3); %correct for reference
        PixelShift = PixelShift + PixelShiftinit; %sum previous shifts
        iter = iter + 1;
        if (sum(abs(PixelShift - PixelShiftinit) > 1) == 0) %if all offsets changed less than 1 pixel, end
            converged = 1;
        elseif (iter == 5) %if reached max iterations
            converged = 1;
            if (sum(abs(PixelShift - PixelShiftinit) < 3) ~= 3) %if all offsets changed less than 3 pixels, keep, else set to 0
                PixelShift = [0; 0; 0]; %move back to no offset
            end
            disp('Did not converge on image offset...')
        else %continue loop
            PixelShiftinit = PixelShift; %set init for next loop
        end
    end
    disp('Determining Image Offset Completed.')
end


%% Proton Image
if (NewImages == 1)
    %Reconstruct Image
    disp('Reconstructing UTE Image...')
    ProtonImageLR = zeros(H_RecMatrix*3,H_RecMatrix*3,H_RecMatrix*3,H_coils);%3x overgridding
    ProtonImage = zeros(H_RecMatrix,H_RecMatrix,H_RecMatrix,H_coils);%3x overgridding
    for i = 1:H_coils
        if strcmp(ScanVersion,'Xe_CTC') || strcmp(ScanVersion,'DukeIPF') %recon 2x interpolated
            ProtonImageLR(:,:,:,i) = GasExchange.Dissolved_HImageRecon_CTC(H_RecMatrix,HKSpace_Steady(:,:,i),HTraj_Steady/2,PixelShift);
            ProtonImage(:,:,:,i) = GasExchange.Dissolved_HighResImageRecon(H_RecMatrix,HKSpace_Steady(:,:,i),HTraj_Steady/2,PixelShift);
        else
            ProtonImageLR(:,:,:,i) = GasExchange.Dissolved_HImageRecon(H_RecMatrix,HKSpace_Steady(:,:,i),HTraj_Steady,PixelShift);
            ProtonImage(:,:,:,i) = GasExchange.Dissolved_HighResImageRecon(H_RecMatrix,HKSpace_Steady(:,:,i),HTraj_Steady,PixelShift);            
        end
    end
    ProtonImageLR = GasExchange.SOS_Coil_Combine(ProtonImageLR);%Combine Coils
    ProtonImage = GasExchange.SOS_Coil_Combine(ProtonImage);%Combine Coils
    disp('Reconstructing UTE Image Competed.')

    %Correct Bias
    disp('Correcting Proton Bias...')
    ProtonImageLR = medfilt3(ProtonImageLR, [7 7 7]);%remove undersampling artifacts
    ProtonImageLR = imgaussfilt3(ProtonImageLR, 0.5);%reduce noise
    [ProtonImageLR, ~] = GasExchange.Dissolved_ProtonBiasCorrection(ProtonImageLR);
    [ProtonImage, ~] = GasExchange.Dissolved_ProtonBiasCorrection(ProtonImage);
    
    % Determine Levels
    ProtonMax = prctile(abs(ProtonImage(:)),99.99);
   
    disp('Proton Bias Corrected.')  
end
    
%View
ProtonMontage = figure('Name','Proton Image');set(ProtonMontage,'WindowState','minimized');
% montage(abs(ProtonImage(H_RecMatrix:2*H_RecMatrix,H_RecMatrix:2*H_RecMatrix,H_RecMatrix:2*H_RecMatrix)),'DisplayRange',[0 ProtonMax])%unregistered for these
montage(abs(ProtonImage),'DisplayRange',[0 ProtonMax])%unregistered for these


%% Spectra
disp('Fitting Spectrum...')
if strcmp(ScanVersion,'Xe_CTC') || strcmp(ScanVersion,'DukeIPF')
    try
        calFile = dir(fullfile(DataLocation,'cal','*.data'));
        [GasExResults, ~] = GasExchange.XeCTC_Calibration(fullfile(calFile.folder,calFile.name));
    catch
        error('Error. Could not analyze calibration Data file')
        success = 0;
    end
    time = GasExResults.DisFit.t;
    AppendedDissolvedNMRFit = GasExResults.DisFit;
    AppendedDissolvedFit = (time(2)-time(1))*fftshift(fft(AppendedDissolvedNMRFit.calcComponentTimeDomainSignal(time),[],1),1);
else
    %Set parameters
    time = dwell_s*(0:XeSpec_nsamp-1);

    %Dissolved Phase fitting
    %initial guesses
    area_guess = [0.27*abs(PostDissolvedFID(1)) abs(PostDissolvedFID(1)) 0.1*abs(PostDissolvedFID(1))]; 
    area_lowerBounds = [0 0 0]; % Always positive
    area_upperBounds = 1E10*[1 1 1]; % Make it finite, but huge
    freq_guess = [534 -126 -7084];
    freq_lowerBounds = [500 -2000 -12000];
    freq_upperBounds = [2000 0 -4000];
    fwhm_guesses = [300 273 44];
    fwhm_lowerBounds = [0 0 0];
    fwhm_upperBounds = inf*[1 1 1];
    fwhmG_guess = [0 275 0];
    fwhmG_lowerBounds = [0 0 0];
    fwhmG_upperBounds = [0 inf 0];
    phase_guesses = [-95 151 9];
    phase_lowerBounds = -inf*[1 1 1];
    phase_upperBounds = inf*[1 1 1];

    %fit prepended dissolved for better guesses
    PrependedDissolvedNMRFit = GasExchange_Spectro.NMR_TimeFit_v(PreDissolvedFID,time,area_guess,freq_guess,fwhm_guesses,fwhmG_guess,phase_guesses,[],[]);
    PrependedDissolvedNMRFit.fitTimeDomainSignal();%initial fit
    PrependedDissolvedNMRFit.setBounds(area_lowerBounds,area_upperBounds,...
        freq_lowerBounds,freq_upperBounds,...
        fwhm_lowerBounds,fwhm_upperBounds,...
        fwhmG_lowerBounds,fwhmG_upperBounds,...
        phase_lowerBounds,phase_upperBounds);%add bounds
    PrependedDissolvedNMRFit.fitTimeDomainSignal();%refit

    %fit appended dissolved
    AppendedDissolvedNMRFit = GasExchange_Spectro.NMR_TimeFit_v(PostDissolvedFID,time,area_guess,PrependedDissolvedNMRFit.freq,PrependedDissolvedNMRFit.fwhm,PrependedDissolvedNMRFit.fwhmG,PrependedDissolvedNMRFit.phase,[],[]);
    AppendedDissolvedNMRFit.fitTimeDomainSignal();%initial fit
    AppendedDissolvedNMRFit.setBounds(area_lowerBounds,area_upperBounds,...
        freq_lowerBounds,freq_upperBounds,...
        fwhm_lowerBounds,fwhm_upperBounds,...
        fwhmG_lowerBounds,fwhmG_upperBounds,...
        phase_lowerBounds,phase_upperBounds);%add bounds
    AppendedDissolvedNMRFit.fitTimeDomainSignal();%refit
    AppendedDissolvedFit = dwell_s*fftshift(fft(AppendedDissolvedNMRFit.calcComponentTimeDomainSignal(time),[],1),1);
end

%View
DissolvedNMR = figure('Name','Dissolved Phase Spectrum');set(DissolvedNMR,'WindowState','minimized');
set(DissolvedNMR,'color','white','Units','inches','Position',[0.25 0.25 16 9])
hold on
%data
plot(AppendedDissolvedNMRFit.f, abs(AppendedDissolvedNMRFit.spectralDomainSignal),'b')%mag
plot(AppendedDissolvedNMRFit.f, real(AppendedDissolvedNMRFit.spectralDomainSignal),'r')%real
plot(AppendedDissolvedNMRFit.f, imag(AppendedDissolvedNMRFit.spectralDomainSignal),'Color',[0,0.6,0.2]')%imag
%fits
plot(AppendedDissolvedNMRFit.f,abs(sum(AppendedDissolvedFit,2)),'Color',[0 0 1 0.33],'LineWidth',3)%mag
plot(AppendedDissolvedNMRFit.f,real(sum(AppendedDissolvedFit,2)),'Color',[1 0 0 0.33],'LineWidth',3)%real
plot(AppendedDissolvedNMRFit.f,imag(sum(AppendedDissolvedFit,2)),'Color',[0,0.6,0.2,0.33]','LineWidth',3)%imag
%components
area(AppendedDissolvedNMRFit.f,real(AppendedDissolvedFit(:,1)),'FaceColor','r','FaceAlpha',0.33,'LineStyle','none')%rbc
area(AppendedDissolvedNMRFit.f,real(AppendedDissolvedFit(:,2)),'FaceColor','b','FaceAlpha',0.33,'LineStyle','none')%barrier
area(AppendedDissolvedNMRFit.f,real(AppendedDissolvedFit(:,3)),'FaceColor',[0,0.6,0.2]','FaceAlpha',0.33,'LineStyle','none')%gas
%settings
lgnd = legend({'Spectrum - Magnitude','Spectrum - Real','Spectrum - Imaginary','Fit - Magnitude','Fit - Real','Fit - Imaginary','RBC - Real','Membrane - Real','Gas - Real'},'Location','best','NumColumns',3);
set(lgnd,'color','none');
xlim([-8000 2000])
set(gca, 'XDir','reverse','FontSize',24)
title(['Dissolved Phase Spectra and Fit: RBC/Membrane = ',num2str(AppendedDissolvedNMRFit.area(1)/AppendedDissolvedNMRFit.area(2))])
xlabel('Frequency (Hz)')
ylabel('NMR Signal (a.u.)')
hold off
disp('Fitting Spectrum Completed.')

%% Remove approch to steady state for Dissolved Phase
%copy to new variable before modifying
DissolvedKSpace_SS = DissolvedKSpace;
GasKSpace_SS = GasKSpace;
XeTraj_SS = XeTraj;
SS_ind = 60;%remove 1st 60 points
DissolvedKSpace_SS(:,1:SS_ind) = [];
GasKSpace_SS(:,1:SS_ind) = [];
XeTraj_SS(:,:,1:SS_ind) = [];


%% Gas Phase Contamination Removal
if(NewImages == 1)
    disp('Removing Gas Phase Contamination...')
    if strcmp(ScanVersion,'Xe_CTC') || strcmp(ScanVersion,'DukeIPF') ...
            || (strcmp(ScanVersion,'CCHMC_V3') && mode(GasKSpace_MaxInd,2) > 3)
       %correction not possible
       if strcmp(ScanVersion,'CCHMC_V3')
            disp([newline 'IMPORTANT INFO...']);
            pause(2)
            disp('...');pause(1); disp('...');pause(1); disp('Too much slow switching hardware failure.');
            pause(2)
            disp('...'); pause(1); disp('...');pause(1); disp('Ommitting Gas Phase Contamination...');
            pause(3)
       end
       CorrectedDissKSpace_SS = DissolvedKSpace_SS;
    else
       CorrectedDissKSpace_SS = GasExchange.GasPhaseContaminationRemoval(DissolvedKSpace_SS,GasKSpace_SS,dwell_s,-freq_jump,AppendedDissolvedNMRFit.phase(3),AppendedDissolvedNMRFit.area(3),GasFlipAngle);
    end
    disp('Removing Gas Phase Contamination Completed.')
end


%% View k0 dynamics
XePulses = (1:size(CorrectedDissKSpace_SS,2))'+SS_ind;%start after steady state

%View H and Xe Dynamic data
SigDynamics = figure('Name','Signal Dynamics');set(SigDynamics,'WindowState','minimized');
set(SigDynamics,'color','white','Units','inches','Position',[0.25 0.25 16 9])
%Xe - SS and Detrending
subplot(1,2,1);
hold on
set(gca,'FontSize',18)
%Dissolved
plot(abs(DissolvedKSpace(1,:)),'Color',[0.5 0 0],'LineWidth',1);
plot(XePulses,abs(CorrectedDissKSpace_SS(1,:)),'Color',[1.0 0 0],'LineWidth',1);
%Gas
plot(abs(GasKSpace(1,:)),'Color',[0 0 1],'LineWidth',1);
%SS
h1 = patch([0 SS_ind SS_ind 0], [max(ylim) max(ylim) 0 0], [0.66 0.66 0.66],'EdgeColor','none');%Remove due to steady state
xl = xlim; yl = ylim;
uistack(h1, 'bottom')
xlim(xl);ylim(yl);
title('Xe Signal Intensity Dynamics')
box on
xlabel('Projection')
ylabel('Signal Intensity')
legend({'Steady State Discard','Dissolved','Corrected Dissolved','Gas'},'Location','best','NumColumns',2,'FontSize',12)
hold off
%H - Motion
subplot(1,2,2);
hold on
set(gca,'FontSize',18)
h1 = plot(HPhaseDynamics,'k','LineWidth',1);
xl = xlim; yl = ylim;
patch([motion_ind length(HPulses) length(HPulses) motion_ind], [max(ylim) max(ylim) min(ylim) min(ylim)], [0.5 0.5 0.5],'EdgeColor','none')%Remove due to motion
uistack(h1, 'top');
plot(HDynFit,'Color',[0 1 0 0.5],'LineWidth',3);
plot(HDynFit+3*HDynRMSE,'Color',[0 1 0 0.5])
plot(HDynFit-3*HDynRMSE,'Color',[0 1 0 0.5])
xlim(xl);ylim(yl);
title('H Signal Phase Dynamics')
box on
xlabel('Projection')
ylabel('Signal Phase (deg)')
legend({'Unused - Motion','Data','Fit'},'Location','best','FontSize',12)
hold off
sgtitle('Signal Dynamics','FontSize',22)


%% Reconstruct Xe Images
if(NewImages == 1)
    %Vent Image
    disp('Reconstructing Ventilation Image...')
    UncorrectedVentImage = GasExchange.Dissolved_HighResImageRecon(Xe_RecMatrix,GasKSpace_SS,XeTraj_SS/2,PixelShift); %2x Resolution
    disp('Reconstructing Ventilation Image Completed.')
    disp('Correcting Ventilation Bias...')
    VentMask = (abs(UncorrectedVentImage))>(prctile(abs(UncorrectedVentImage(:)),95)/3);
    VentMask = imerode(VentMask,strel('sphere',5));
    Xeregions = regionprops3(logical(VentMask), 'all'); %get information on remaining regions
    Xeregions = sortrows(Xeregions,1,'descend'); %order by size
    VentMask = zeros(size(VentMask)); %restart mask
    try
        VentMask(Xeregions.VoxelIdxList{1,1}) = 1; %select largest region only
        if (size(Xeregions,1)>1)
            if (Xeregions.Volume(2) > 0.33*Xeregions.Volume(1)) %Select 2nd largest region if similar size
                VentMask(Xeregions.VoxelIdxList{2,1}) = 1;
            end
        end
        VentMask = logical(imdilate(VentMask,strel('sphere',6)));
    catch
        %if no signal for mask to be made, create mask of one
        success = 0;%mark as unsuccessful
        VentMask = ones(size(VentMask));
    end
    [VentImage, ~] = GasExchange.Vent_BiasCorrection(UncorrectedVentImage,VentMask);
    disp('Ventilation Bias Corrected.')

    %Gas Image
    disp('Reconstructing Gas Image...')
    GasImage = GasExchange.Dissolved_LowResImageRecon(Xe_RecMatrix,GasKSpace_SS,XeTraj_SS/2,PixelShift); %2x Resolution
    disp('Reconstructing Gas Image Completed.')

    %Dissolved Image
    disp('Reconstructing Dissolved Image...')
    DissolvedImage = GasExchange.Dissolved_LowResImageRecon(Xe_RecMatrix,DissolvedKSpace_SS,XeTraj_SS/2,PixelShift); %2x Resolution
    disp('Reconstructing Dissolved Image Completed.')
    disp('Reconstructing Corrected Dissolved Image...')
    CorrDissolvedImage = GasExchange.Dissolved_LowResImageRecon(Xe_RecMatrix,CorrectedDissKSpace_SS,XeTraj_SS/2,PixelShift); %2x Resolution
    disp('Reconstructing Corrected Dissolved Image Completed.')
elseif (VentMask == ones(size(VentMask)))%Vent mask already exists since not 
    %make success = 0 since indicates no signal
    success = 0;%mark as unsuccessful
end

%View
VentMontage = figure('Name','Vent Image');set(VentMontage,'WindowState','minimized');
montage(abs(UncorrectedVentImage),'DisplayRange',[0 max(abs(UncorrectedVentImage(:)))])

CorrVentMontage = figure('Name','Corrected Vent Image');set(CorrVentMontage,'WindowState','minimized');
montage(abs(VentImage),'DisplayRange',[0 max(abs(VentImage(:)))])

GasMontage = figure('Name','Gas Image');set(GasMontage,'WindowState','minimized');
montage(abs(GasImage),'DisplayRange',[0 max(abs(GasImage(:)))])

DissolvedMontage = figure('Name','Dissolved Image');set(DissolvedMontage,'WindowState','minimized');
montage(abs(DissolvedImage),'DisplayRange',[0 max(abs(DissolvedImage(:)))])

CorrDissolvedMontage = figure('Name','Corrected Dissolved Image');set(CorrDissolvedMontage,'WindowState','minimized');
montage(abs(CorrDissolvedImage),'DisplayRange',[])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Process Data                                                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Obtain Parameters of Interest from Spectrum
SpecRBCBarrierRatio = AppendedDissolvedNMRFit.area(1)/AppendedDissolvedNMRFit.area(2); %To determine expected ratio in image, only use appended to be closer to be at steady state
%Determine actual separation angle
SpecDissolvedPhaseDifference = AppendedDissolvedNMRFit.phase(2) - AppendedDissolvedNMRFit.phase(1); %To determine accuracy of separation at k=0
if SpecDissolvedPhaseDifference >= 180  %Shift phase difference to be between +180 and -180
    SpecDissolvedPhaseDifference = SpecDissolvedPhaseDifference - 360;
elseif SpecDissolvedPhaseDifference <= -180
    SpecDissolvedPhaseDifference = SpecDissolvedPhaseDifference + 360;
end
%T2* Calculations
GasT2Star = 1/(pi * AppendedDissolvedNMRFit.fwhm(3))*1000; %in ms
BarrierT2Star = 1/(pi * max([AppendedDissolvedNMRFit.fwhm(2),AppendedDissolvedNMRFit.fwhmG(2)]))*1000; %in ms
BarrierLtoGRatio = AppendedDissolvedNMRFit.fwhm(2)/AppendedDissolvedNMRFit.fwhmG(2);
RBCT2Star = 1/(pi * AppendedDissolvedNMRFit.fwhm(1))*1000; %in ms
%Freqs
GasFreq = AppendedDissolvedNMRFit.freq(3)+freq_jump;
BarrierFreq = AppendedDissolvedNMRFit.freq(2)+freq_jump;
RBCFreq = AppendedDissolvedNMRFit.freq(1)+freq_jump;
%Contamination
SpecGasPhaseContaminationPercentage = AppendedDissolvedNMRFit.area(3)/(AppendedDissolvedNMRFit.area(1)+AppendedDissolvedNMRFit.area(2))*100;
 

%% Make Proton Mask
if (strcmp(ScanVersion,'CCHMC_V3'))
    automask = 1; % using ML trained model to create a mask
else
    automask = 0; % no ML trained model to create a mask
end

if(NewMask==1)
    if automask==1
        %=========================registration=============================
        %Set up registration
        disp('Registering Proton Image to Gas Image...')
        [optimizer,metric] = imregconfig('multimodal');
        
        InitMask = ~imbinarize(ProtonImageLR/(prctile(abs(ProtonImageLR(:)),99.9)));%initial otsu thresholding; select low intensity
        ProtonSignal = mean(ProtonImageLR(~InitMask));%determine mean proton signal for thresholding
        %Make mask of noise
        NoiseMask = InitMask;
        NoiseMask(ProtonImageLR>0.25*ProtonSignal) = 0; %remove clear proton signal to disconnect noise and lung

        %Modify images for better registration
        VentReg = medfilt3(abs(UncorrectedVentImage), [3 3 3]);
        ProtonReg = ProtonImageLR;
        ProtonReg = medfilt3(abs(ProtonReg), [5 5 5]);
        ProtonReg(ProtonImageLR>ProtonSignal) = ProtonMax;%make anything clearly proton the max

        %Make mask of noise
        NoiseMask = InitMask;
        NoiseMask(ProtonImageLR>0.25*ProtonSignal) = 0; %remove clear proton signal to disconnect noise and lung
        temp = imclearborder(NoiseMask);%remove noise
        NoiseMask = NoiseMask - temp;%leave only noise
        NoiseMask = imclose(NoiseMask,strel('sphere',7));%smooth out holes
        NoiseMask = imdilate(NoiseMask,strel('sphere',7));%expand to avoid edges
        NoiseMask = imclose(NoiseMask,strel('sphere',7));%smooth out holes
        NoiseMask = logical(NoiseMask);
        %Initial Lung Mask
        LungInitMask = InitMask - NoiseMask;
        LungInitMask(LungInitMask<0) = 0;
        LungInitMask = logical(LungInitMask);
        LungStd =  std(ProtonImageLR(LungInitMask));%determine mean lung signal for thresholding
        LungInitMask = imerode(LungInitMask,strel('sphere',5));%only for-certain lung
        LungSignal =  mean(ProtonImageLR(LungInitMask));%determine mean lung signal for thresholding

        ProtonReg(ProtonImageLR<LungSignal) = LungSignal-(3*LungStd);%make anything clearly lung the min
        ProtonReg(NoiseMask) = ProtonMax;%make anything clearly noise the max
        ProtonReg = abs(ProtonReg - ProtonMax);%Invert to hopefully help
    
        %create base tform since using 3x overgridded H image
        tformBase = affine3d([1,0,0,0;0,1,0,0;0,0,1,0;-H_RecMatrix,-H_RecMatrix,-H_RecMatrix,1]);
        
        %start with similarity transform
        optimizer.GrowthFactor = 1.025;
        optimizer.Epsilon = 1e-10;
        optimizer.InitialRadius = 0.00005;
        optimizer.MaximumIterations = 150;
        tformSimilarity = imregtform(abs(ProtonReg), abs(VentReg), 'similarity', optimizer, metric,'InitialTransformation',tformBase); %initial similarity transform
    
        %affine transform
        optimizer.GrowthFactor = 1.01;
        optimizer.Epsilon = 1e-12;
        optimizer.InitialRadius = 0.00005;
        optimizer.MaximumIterations = 500;
        tform = imregtform(abs(ProtonReg), abs(VentReg), 'affine', optimizer, metric,'InitialTransformation',tformSimilarity,'PyramidLevels',4); %affine transform using inital similarity transform
        ProtonImageLRRegistered = imwarp(abs(ProtonImageLR), tform,'OutputView',imref3d(size(abs(VentImage))));
   
        sizeFirstVariable = [H_RecMatrix, H_RecMatrix, H_RecMatrix];
        sizeSecondVariable = sizeFirstVariable * 3;
        ProtonImageHR3 = zeros(sizeSecondVariable);
        centerIndices = ceil(sizeSecondVariable / 2) - ceil(sizeFirstVariable / 2);
        ProtonImageHR3(centerIndices(1):(centerIndices(1)+sizeFirstVariable(1)-1), ...
                      centerIndices(2):(centerIndices(2)+sizeFirstVariable(2)-1), ...
                      centerIndices(3):(centerIndices(3)+sizeFirstVariable(3)-1)) = ProtonImage;
        ProtonImageRegistered = imwarp(abs(ProtonImageHR3), tform,'OutputView',imref3d(size(abs(VentImage))));
        
        disp('Registering Proton Image to Gas Image Completed.')
        %======================================================================
        
        %============================Making Lung Mask==========================
        disp('Making Lung Mask (trained ML model method)...');

        % create a folder if not exist
        automasking_folder = 'automasking_folder';
        destinationFolderPath = join([FunctionDirectory(1:3),automasking_folder]);
        cd(FunctionDirectory(1:3));
        if ~exist(automasking_folder, 'dir')
            mkdir(destinationFolderPath);
        else
            disp('Folder already exists in the destination folder.');
        end

        % List of filenames to copy
        fileNames = {
            'Segmentation3D_Module.py'
            'GasExchange_3D_HLR_100e_20250330.hdf5'
            'GasExchange_3D_Xe_100e_20250325.hdf5'
            'GasExchange_3D_XeHLR_100e_20250330.hdf5'
        };
        sourceFolderPath = [FunctionDirectory,'\+GasExchange'];
        
        % Loop over each file
        for i = 1:length(fileNames)
            fileName = fileNames{i};
            sourceFilePath = fullfile(sourceFolderPath, fileName);
            destFilePath = fullfile(destinationFolderPath, fileName);
        
            if ~exist(destFilePath, 'file')
                copyfile(sourceFilePath, destinationFolderPath);
                fprintf('Copied: %s\n', fileName);
            else
                fprintf('Already exists: %s\n', fileName);
            end
        end

        % copy model to the automasking folder

%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %Need a raw CPython installation (NOT anaconda).  These commands below
%         % shouldn't need to be called if everthing is set up properly.
        % terminate(pyenv)
        % pyenv('Version','C:\Users\MCM5BK\AppData\Local\Programs\Python\Python39\pythonw.exe'); %3.10.2
        % %Call Python 3.10.2
        % % pyenv('Version', pythonPath);
        % system('pip install numpy')
        % system('pip install keras==2.10.0'); % Specific version of Keras
        % system('pip install tensorflow==2.10.1'); % Specific version of TensorFlow
        % system('pip install nibabel')
        % system('pip install scipy') 
        % terminate(pyenv)
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        cd(FunctionDirectory)
        % Define full file path
        fileToDelete = fullfile(destinationFolderPath, 'Mat2Py_preprocessing.mat');
        
        % Check and delete
        if exist(fileToDelete, 'file') == 2
            delete(fileToDelete);
            disp('Mat2Py_preprocessing.mat deleted successfully.');
        else
            disp('Mat2Py_preprocessing.mat does not exist in the folder.');
        end

          %Write Xe and Proton images to external file for Python-based ML segmentation
        ImageInput = 'H'; % 'Xe' | 'H' | 'XeH'
        Mat2Py_Images = GasExchange.preprocess_images_for_auto_segmentation(VentImage,ProtonImageRegistered, ImageInput, destinationFolderPath);
        % figure; imslice(squeeze(Mat2Py_Images));
        cd(destinationFolderPath);
        % pathToMod = fileparts(which('+GasExchange\Segmentation3D_Module.py'));
        pathToMod = fileparts(which('Segmentation3D_Module.py'));
        if count(py.sys.path,pathToMod)==0
            insert(py.sys.path,int32(0),pathToMod)
            disp('Uploaded Segmentation3DModule')
        end
        
        % SCRIPTLOC = join([FunctionDirectory,'\+GasExchange']);
        % cd(SCRIPTLOC)
        SCRIPTLOC = destinationFolderPath;
        switch ImageInput
            case 'Xe'
                MODELLOC = join([destinationFolderPath,'\GasExchange_3D_Xe_100e_20250325.hdf5']);
            case 'H'
                MODELLOC = join([destinationFolderPath,'\GasExchange_3D_HLR_100e_20250330.hdf5']);
            case 'XeH'
                MODELLOC = join([destinationFolderPath,'\GasExchange_3D_XeHLR_100e_20250330.hdf5']);
        end        
            IMGLOC = join([destinationFolderPath,'\Mat2Py_preprocessing.mat']);

        % run python script 
        cd(SCRIPTLOC)
        command = string(strcat('python Segmentation3D_Module.py',{' '},IMGLOC,{' '},MODELLOC));
        LungMask = system(command);
        % load auto mask
        cd(destinationFolderPath);
        if exist('AutoMask.mat')
            load([destinationFolderPath,'\AutoMask.mat']);
            LungMask = squeeze(AutoMask(1,:,:,:));
            for i = 1:size(LungMask,3)
                sliceSum = sum(LungMask(:,:,i), 'all');
                if sliceSum <= 10
                    LungMask(:,:,i) = 0;
                end
            end
            disp('auto mask process Completed.')
            ProtonMaskRegistered = LungMask>0;
            %Determine Slices to Plot
                %Coronal
            Co_MaskVoxels = squeeze(sum(ProtonMaskRegistered,[1,2]));
            Co_StartIndex = find(Co_MaskVoxels,1,'first');
            Co_EndIndex = find(Co_MaskVoxels,1,'last');
            Co_MiddleIndex = round((Co_StartIndex+Co_EndIndex)/2);
            Co_Step = floor((Co_EndIndex-Co_StartIndex)/NumPlotSlices);
            Slices_Co = (Co_MiddleIndex-Co_Step*(NumPlotSlices-1)/2):Co_Step:(Co_MiddleIndex+Co_Step*(NumPlotSlices-1)/2);
                %Axial
            Ax_MaskVoxels = squeeze(sum(ProtonMaskRegistered,[2,3]));
            Ax_StartIndex = find(Ax_MaskVoxels,1,'first');
            Ax_EndIndex = find(Ax_MaskVoxels,1,'last');
            Ax_MiddleIndex = round((Ax_StartIndex+Ax_EndIndex)/2);
            Ax_Step = floor((Ax_EndIndex-Ax_StartIndex)/NumPlotSlices);
            Slices_Ax = (Ax_MiddleIndex-Ax_Step*(NumPlotSlices-1)/2):Ax_Step:(Ax_MiddleIndex+Ax_Step*(NumPlotSlices-1)/2);
            delete([destinationFolderPath,'\AutoMask.mat'])
        else
            disp('Auto mask does not exist in the destination folder.');
        end
        cd(DataLocation)
    else
        disp('Making Lung Mask...')
        %Initial mask
        InitMask = ~imbinarize(ProtonImage/(prctile(abs(ProtonImage(:)),99.9)));%initial otsu thresholding; select low intensity
        ProtonSignal = mean(ProtonImage(~InitMask));%determine mean proton signal for thresholding
        %Make mask of noise
        NoiseMask = InitMask;
        NoiseMask(ProtonImage>0.25*ProtonSignal) = 0; %remove clear proton signal to disconnect noise and lung
        temp = imclearborder(NoiseMask);%remove noise
        NoiseMask = NoiseMask - temp;%leave only noise
        NoiseMask = imclose(NoiseMask,strel('sphere',7));%smooth out holes
        NoiseMask = imdilate(NoiseMask,strel('sphere',7));%expand to avoid edges
        NoiseMask = imclose(NoiseMask,strel('sphere',7));%smooth out holes
        NoiseMask = logical(NoiseMask);
        %Initial Lung Mask
        LungInitMask = InitMask - NoiseMask;
        LungInitMask(LungInitMask<0) = 0;
        LungInitMask = logical(LungInitMask);
        LungStd =  std(ProtonImage(LungInitMask));%determine mean lung signal for thresholding
        LungInitMask = imerode(LungInitMask,strel('sphere',5));%only for-certain lung
        LungSignal =  mean(ProtonImage(LungInitMask));%determine mean lung signal for thresholding
        LungInitMask = zeros(size(LungInitMask));
        LungInitMask(ProtonImage<(LungSignal+(ProtonSignal-LungSignal)/2)) = 1;
        LungInitMask = LungInitMask - NoiseMask;
        LungInitMask(LungInitMask<0) = 0;
        Hregions = regionprops3(logical(LungInitMask), 'all'); %get information on remaining regions
        Hregions = sortrows(Hregions,1,'descend'); %order by size
        LungInitMask = zeros(size(LungInitMask)); %start final mask
        for ii=1:size(Hregions,1)
            if (Hregions.Volume(ii) > 0.33*Hregions.Volume(1)) %Select only largest regions
                LungInitMask(Hregions.VoxelIdxList{ii,1}) = 1;
            else
                break
            end
        end
        %Final touches of Lung Mask
        LungMask = LungInitMask;
        LungMask(ProtonImage>(LungSignal+5*LungStd)) = 0;%remove clear proton
        Hregions = regionprops3(logical(LungInitMask), 'all'); %get information on remaining regions
        Hregions = sortrows(Hregions,1,'descend'); %order by size
        LungMask = zeros(size(LungMask)); %start final mask
        for ii=1:size(Hregions,1)
            if (Hregions.Volume(ii) > 0.33*Hregions.Volume(1)) %Select only largest regions
                LungMask(Hregions.VoxelIdxList{ii,1}) = 1;
            else
                break
            end
        end
        LungMask = imerode(logical(LungMask),strel('sphere',2));%avoid edges
        LungMask = imclose(logical(LungMask),strel('sphere',4));%avoid edges
    end
else
    disp('Importing Lung Mask...')
    %Calculate Proton and Lung signals
    InitMask = ~imbinarize(ProtonImage/(prctile(abs(ProtonImage(:)),99.9)));%initial otsu thresholding; select low intensity
    ProtonSignal = mean(ProtonImage(~InitMask));%determine mean proton signal for thresholding
    NoiseMask = InitMask;
    NoiseMask(ProtonImage>0.25*ProtonSignal) = 0; %remove clear proton signal to disconnect noise and lung
    temp = imclearborder(NoiseMask);%remove noise
    NoiseMask = NoiseMask - temp;%leave only noise
    NoiseMask = imclose(NoiseMask,strel('sphere',7));%smooth out holes
    NoiseMask = imdilate(NoiseMask,strel('sphere',7));%expand to avoid edges
    NoiseMask = imclose(NoiseMask,strel('sphere',7));%smooth out holes
    NoiseMask = logical(NoiseMask);
    LungInitMask = InitMask - NoiseMask;
    LungInitMask(LungInitMask<0) = 0;
    LungInitMask = logical(LungInitMask);
    LungStd =  std(ProtonImage(LungInitMask));
    LungInitMask = imerode(LungInitMask,strel('sphere',5));%only for-certain lung
    LungSignal =  mean(ProtonImage(LungInitMask));%determine mean lung signal for thresholding
    %Import Mask
    LungMask = logical(fliplr(rot90(niftiread([DataLocation,'\ProtonMask.nii.gz']),-1)));
end

if LungMask == zeros(size(LungMask))
    %if no mask, create mask of ones
    success = 0;%mark as unsuccessful
    LungMask = logical(ones(size(LungMask)));
end

%View
ProtonMaskMontage = figure('Name','Lung Mask');set(ProtonMaskMontage,'WindowState','minimized');
if automask==1
    montage(LungMask,'DisplayRange',[])%unregistered for these
else
    montage(LungMask(H_RecMatrix:2*H_RecMatrix,H_RecMatrix:2*H_RecMatrix,H_RecMatrix:2*H_RecMatrix),'DisplayRange',[])%unregistered for these
end
disp('Lung Mask Completed.')


%% Register Proton to Xenon
if exist('tform','var') == 0 %if doesn't exist
    tform = affine3d([1,0,0,0;0,1,0,0;0,0,1,0;-H_RecMatrix,-H_RecMatrix,-H_RecMatrix,1]);%initialize transform if not present
end
if (success==0 && automask == 0)%if already failed, dont run full registration as it likely will only hurt with low signal
    disp('Registering Proton Image to Gas Image...')
    %create base tform since using 3x overgridded H image
    tformBase = affine3d([1,0,0,0;0,1,0,0;0,0,1,0;-H_RecMatrix,-H_RecMatrix,-H_RecMatrix,1]);
    
    %transform according to base
    ProtonImageRegistered = imwarp(abs(ProtonImage), tformBase,'OutputView',imref3d(size(abs(VentImage))));
    ProtonMaskRegistered = logical(imwarp(abs(LungMask), tformBase,'OutputView',imref3d(size(abs(VentImage)))));
    disp('Registering Proton Image to Gas Image Completed.')

    %Determine Slices to Plot
        %Coronal
    Co_MaskVoxels = squeeze(sum(ProtonMaskRegistered,[1,2]));
    Co_StartIndex = find(Co_MaskVoxels,1,'first');
    Co_EndIndex = find(Co_MaskVoxels,1,'last');
    Co_MiddleIndex = round((Co_StartIndex+Co_EndIndex)/2);
    Co_Step = floor((Co_EndIndex-Co_StartIndex)/NumPlotSlices);
    Slices_Co = (Co_MiddleIndex-Co_Step*(NumPlotSlices-1)/2):Co_Step:(Co_MiddleIndex+Co_Step*(NumPlotSlices-1)/2);
        %Axial
    Ax_MaskVoxels = squeeze(sum(ProtonMaskRegistered,[2,3]));
    Ax_StartIndex = find(Ax_MaskVoxels,1,'first');
    Ax_EndIndex = find(Ax_MaskVoxels,1,'last');
    Ax_MiddleIndex = round((Ax_StartIndex+Ax_EndIndex)/2);
    Ax_Step = floor((Ax_EndIndex-Ax_StartIndex)/NumPlotSlices);
    Slices_Ax = (Ax_MiddleIndex-Ax_Step*(NumPlotSlices-1)/2):Ax_Step:(Ax_MiddleIndex+Ax_Step*(NumPlotSlices-1)/2);
elseif(NewReg == 1 && automask == 0)
    %move to imregdemons?
    %Set up registration
    disp('Registering Proton Image to Gas Image...')
    [optimizer,metric] = imregconfig('multimodal');
    
    %Modify images for better registration
    VentReg = medfilt3(abs(UncorrectedVentImage), [3 3 3]);
    ProtonReg = ProtonImage;
    ProtonReg = medfilt3(abs(ProtonReg), [5 5 5]);
    ProtonReg(ProtonImage>ProtonSignal) = ProtonMax;%make anything clearly proton the max
    ProtonReg(ProtonImage<LungSignal) = LungSignal-(3*LungStd);%make anything clearly lung the min
    ProtonReg(NoiseMask) = ProtonMax;%make anything clearly noise the max
    ProtonReg = abs(ProtonReg - ProtonMax);%Invert to hopefully help

    %create base tform since using 3x overgridded H image
    tformBase = affine3d([1,0,0,0;0,1,0,0;0,0,1,0;-H_RecMatrix,-H_RecMatrix,-H_RecMatrix,1]);
    
    %start with similarity transform
    optimizer.GrowthFactor = 1.025;
    optimizer.Epsilon = 1e-10;
    optimizer.InitialRadius = 0.00005;
    optimizer.MaximumIterations = 150;
    tformSimilarity = imregtform(abs(ProtonReg), abs(VentReg), 'similarity', optimizer, metric,'InitialTransformation',tformBase); %initial similarity transform

    %affine transform
    optimizer.GrowthFactor = 1.01;
    optimizer.Epsilon = 1e-12;
    optimizer.InitialRadius = 0.00005;
    optimizer.MaximumIterations = 500;
    tform = imregtform(abs(ProtonReg), abs(VentReg), 'affine', optimizer, metric,'InitialTransformation',tformSimilarity,'PyramidLevels',4); %affine transform using inital similarity transform
    ProtonImageRegistered = imwarp(abs(ProtonImage), tform,'OutputView',imref3d(size(abs(VentImage))));
    ProtonMaskRegistered = logical(imwarp(abs(LungMask), tform,'OutputView',imref3d(size(abs(VentImage)))));
    disp('Registering Proton Image to Gas Image Completed.')

    %Determine Slices to Plot
        %Coronal
    Co_MaskVoxels = squeeze(sum(ProtonMaskRegistered,[1,2]));
    Co_StartIndex = find(Co_MaskVoxels,1,'first');
    Co_EndIndex = find(Co_MaskVoxels,1,'last');
    Co_MiddleIndex = round((Co_StartIndex+Co_EndIndex)/2);
    Co_Step = floor((Co_EndIndex-Co_StartIndex)/NumPlotSlices);
    Slices_Co = (Co_MiddleIndex-Co_Step*(NumPlotSlices-1)/2):Co_Step:(Co_MiddleIndex+Co_Step*(NumPlotSlices-1)/2);
        %Axial
    Ax_MaskVoxels = squeeze(sum(ProtonMaskRegistered,[2,3]));
    Ax_StartIndex = find(Ax_MaskVoxels,1,'first');
    Ax_EndIndex = find(Ax_MaskVoxels,1,'last');
    Ax_MiddleIndex = round((Ax_StartIndex+Ax_EndIndex)/2);
    Ax_Step = floor((Ax_EndIndex-Ax_StartIndex)/NumPlotSlices);
    Slices_Ax = (Ax_MiddleIndex-Ax_Step*(NumPlotSlices-1)/2):Ax_Step:(Ax_MiddleIndex+Ax_Step*(NumPlotSlices-1)/2);
elseif(automask == 0 || (automask == 1 && NewMask==0)) %if not new reg, transform according to previous transform
    %only warp according to previously determined tform
    ProtonImageRegistered = imwarp(abs(ProtonImage), tform,'OutputView',imref3d(size(abs(VentImage))));
    if(automask == 0)
        ProtonMaskRegistered = logical(imwarp(abs(LungMask), tform,'OutputView',imref3d(size(abs(VentImage)))));
    else
        ProtonMaskRegistered = LungMask;
    end
    disp('Registering Proton Image to Gas Image Completed.')

    %Determine Slices to Plot
        %Coronal
    Co_MaskVoxels = squeeze(sum(ProtonMaskRegistered,[1,2]));
    Co_StartIndex = find(Co_MaskVoxels,1,'first');
    Co_EndIndex = find(Co_MaskVoxels,1,'last');
    Co_MiddleIndex = round((Co_StartIndex+Co_EndIndex)/2);
    Co_Step = floor((Co_EndIndex-Co_StartIndex)/NumPlotSlices);
    Slices_Co = (Co_MiddleIndex-Co_Step*(NumPlotSlices-1)/2):Co_Step:(Co_MiddleIndex+Co_Step*(NumPlotSlices-1)/2);
        %Axial
    Ax_MaskVoxels = squeeze(sum(ProtonMaskRegistered,[2,3]));
    Ax_StartIndex = find(Ax_MaskVoxels,1,'first');
    Ax_EndIndex = find(Ax_MaskVoxels,1,'last');
    Ax_MiddleIndex = round((Ax_StartIndex+Ax_EndIndex)/2);
    Ax_Step = floor((Ax_EndIndex-Ax_StartIndex)/NumPlotSlices);
    Slices_Ax = (Ax_MiddleIndex-Ax_Step*(NumPlotSlices-1)/2):Ax_Step:(Ax_MiddleIndex+Ax_Step*(NumPlotSlices-1)/2);
end

%View
    Registration = figure('Name','Registration','units','normalized','outerposition',[0 0 1 4/NumPlotSlices]);set(Registration,'WindowState','minimized');
set(Registration,'color','white','Units','inches','Position',[0.25 0.25 2*NumPlotSlices 4*2+1])
tiledlayout(4,NumPlotSlices,'TileSpacing','none','Padding','compact');
for slice=1:NumPlotSlices
    nexttile
    if automask == 1
        imshowpair(abs(ProtonImage(:,:,Slices_Co(slice)))/ProtonMax,abs(VentImage(:,:,Slices_Co(slice)))/max(abs(VentImage(ProtonMaskRegistered(:)))),'Scaling','none','ColorChannels','red-cyan')
    else
        imshowpair(abs(ProtonImage(H_RecMatrix:2*H_RecMatrix,H_RecMatrix:2*H_RecMatrix,H_RecMatrix+Slices_Co(slice)))/ProtonMax,abs(VentImage(:,:,Slices_Co(slice)))/max(abs(VentImage(ProtonMaskRegistered(:)))),'Scaling','none','ColorChannels','red-cyan')
    end
    if slice == round(NumPlotSlices/2)
        title('Before Registration','FontSize',24)
    end
end
for slice=1:NumPlotSlices
    nexttile
    if automask == 1
        imshowpair(fliplr(rot90(squeeze(abs(ProtonImage(Slices_Ax(slice),:,:))),-1))/ProtonMax,fliplr(rot90(squeeze(abs(VentImage(Slices_Ax(slice),:,:))),-1))/max(abs(VentImage(ProtonMaskRegistered(:)))),'Scaling','none','ColorChannels','red-cyan')
    else
        imshowpair(fliplr(rot90(squeeze(abs(ProtonImage(H_RecMatrix+Slices_Ax(slice),H_RecMatrix:2*H_RecMatrix,H_RecMatrix:2*H_RecMatrix))),-1))/ProtonMax,fliplr(rot90(squeeze(abs(VentImage(Slices_Ax(slice),:,:))),-1))/max(abs(VentImage(ProtonMaskRegistered(:)))),'Scaling','none','ColorChannels','red-cyan')
    end    
end
for slice=1:NumPlotSlices
    nexttile
    imshowpair(abs(ProtonImageRegistered(:,:,Slices_Co(slice)))/ProtonMax,abs(VentImage(:,:,Slices_Co(slice)))/max(abs(VentImage(ProtonMaskRegistered(:)))),'Scaling','none','ColorChannels','red-cyan')
    if slice == round(NumPlotSlices/2)
        title('After Registration','FontSize',24)
    end
end
for slice=1:NumPlotSlices
    nexttile
    imshowpair(fliplr(rot90(squeeze(abs(ProtonImageRegistered(Slices_Ax(slice),:,:))),-1))/ProtonMax,fliplr(rot90(squeeze(abs(VentImage(Slices_Ax(slice),:,:))),-1))/max(abs(VentImage(ProtonMaskRegistered(:)))),'Scaling','none','ColorChannels','red-cyan')
end


%% Create Dixon Images
%Determine Keys for RBC Oscillations
cd(FunctionDirectory)
disp('Analyzing RBC Osclications...')
[RBCOsc_High_Key,RBCOsc_Low_Key,DissolvedKSpace_SS_Detrend_Normalized,RBC2Bar_struct,OscWorkFlow_Fig] = GasExchange.dissolved_RBC_Osc_Keyhole(XeTraj_SS,GasKSpace_SS,CorrectedDissKSpace_SS,XeTR,SpecRBCBarrierRatio);
disp('Analyzing RBC Osclications Completed.')

%Reconstruct RBC Oscillation Images
disp('Reconstructing RBC keyhole Images...')
RBCOsc_High_Image = GasExchange.Dissolved_RBCOscImageRecon(Xe_RecMatrix,RBCOsc_High_Key,XeTraj_SS/2,PixelShift);%2x Resolution
RBCOsc_Low_Image = GasExchange.Dissolved_RBCOscImageRecon(Xe_RecMatrix,RBCOsc_Low_Key,XeTraj_SS/2,PixelShift);%2x Resolution
RBCOsc_Normalization = GasExchange.Dissolved_RBCOscImageRecon(Xe_RecMatrix,DissolvedKSpace_SS_Detrend_Normalized,XeTraj_SS/2,PixelShift);%2x Resolution
disp('Reconstructing RBC keyhole Images Completed.')

%Dixon Separation of Dissolved Phases
disp('Dixon separation of Dissolved Phases...')
    %Standard Dissolved Image
[BarrierImage,RBCImage,Delta_angle_deg,~] = GasExchange.SinglePointDixon(CorrDissolvedImage,-SpecRBCBarrierRatio,GasImage,ProtonMaskRegistered);
[~,RBCOsc_High,~,~] = GasExchange.SinglePointDixon(RBCOsc_High_Image,-RBC2Bar_struct.High,GasImage,ProtonMaskRegistered);
[~,RBCOsc_Low,~,~] = GasExchange.SinglePointDixon(RBCOsc_Low_Image,-RBC2Bar_struct.Low,GasImage,ProtonMaskRegistered);
[~,RBCOsc_Norm,~,~] = GasExchange.SinglePointDixon(RBCOsc_Normalization,-RBC2Bar_struct.Tot,GasImage,ProtonMaskRegistered);
disp('Dixon separation of Dissolved Phases Completed.')

%Remove Values below zero
BarrierImage(BarrierImage<0) = 0;
RBCImage(RBCImage<0) = 0;
RBCOsc_High(RBCOsc_High<0) = 0;
RBCOsc_Low(RBCOsc_Low<0) = 0;
RBCOsc_Norm(RBCOsc_Norm<0) = 0;

%View
BarrierMontage = figure('Name','Barrier Image');set(BarrierMontage,'WindowState','minimized');
montage(abs(BarrierImage),'DisplayRange',[])
RBCMontage = figure('Name','RBC Image');set(RBCMontage,'WindowState','minimized');
montage(abs(RBCImage),'DisplayRange',[])


%% Quantitative Corrections
disp('Correcting Image Intensities for Comparison...')
%T2* scaling
%M0 = Mt/(e^(-t/T2*))
GasImageScaled = GasImage/(exp(-ActTE90/GasT2Star));
RBCfraction = AppendedDissolvedNMRFit.area(1)/(AppendedDissolvedNMRFit.area(2)+AppendedDissolvedNMRFit.area(1));
DissolvedT2Star = BarrierT2Star*(1-RBCfraction) + RBCT2Star*RBCfraction;
DissolvedImageCorrected = CorrDissolvedImage/(exp(-ActTE90/DissolvedT2Star));
BarrierImageCorrected = BarrierImage/(exp(-ActTE90/BarrierT2Star));
RBCImageCorrected = RBCImage/(exp(-ActTE90/RBCT2Star));

%Flip Angle Correction
GasImageScaled = GasImageScaled*sind(DisFlipAngle)/sind(GasFlipAngle); %10.1002/mp.12264
disp('Correcting Image Intensities for Comparison Completed.')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Analyze Data                                                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate Ratios
disp('Calculating Ratios...')
DissGasRatio = abs(DissolvedImageCorrected)./abs(GasImageScaled);
BarrGasRatio = abs(BarrierImageCorrected)./abs(GasImageScaled);
RBCGasRatio = abs(RBCImageCorrected)./abs(GasImageScaled);
RBCBarrRatio = abs(RBCImageCorrected)./abs(BarrierImageCorrected);

RBCOsc = (abs(RBCOsc_High) - abs(RBCOsc_Low))/mean(abs(RBCOsc_Norm(and(ProtonMaskRegistered,RBCGasRatio>RBCThresh(1)))))*100;%percent
disp('Calculating Ratios Completed.')


%% Scale Ventilation Intensities
disp('Scaling Ventilation Image...')
VentMax = prctile(abs(VentImage(ProtonMaskRegistered(:))),99);
ScaledVentImage = VentImage/VentMax;
ScaledVentImage(ScaledVentImage>1) = 1;
disp('Scaling Ventilation Image Completed.')


%% Calculate Binned Maps
disp('Binning Ratios...')
%Create Binning Maps
SixBinMap = [1 0 0; 1 0.7143 0; 0.4 0.7 0.4; 0 1 0; 0 0.57 0.71; 0 0 1]; %Used for Vent and RBC
EightBinMap = [1 0 0; 1 0.7143 0; 0.4 0.7 0.4; 0 1 0; 184/255 226/255 145/255; 243/255 205/255 213/255; 225/255 129/255 162/255; 197/255 27/255 125/255]; %Used for barrier
SixBinRBCBarMap = [197/255 27/255 125/255; 225/255 129/255 162/255; 0.4 0.7 0.4; 0 1 0; 1 0.7143 0; 1 0 0]; %Used for RBC/Barrier ratio

%Vent Binning
VentBinMap = GasExchange.BinImages(ScaledVentImage, VentThresh);
VentBinMap = VentBinMap.*ProtonMaskRegistered;%Mask to lung volume
VentBinMask = logical((VentBinMap-1).*ProtonMaskRegistered); %Mask for ventilated areas only

if ProtonMaskRegistered == ones(size(ProtonMaskRegistered))
    %if proton mask is entire FOV, dont mask
    success = 0;%mark as unsuccessful
    VentBinMask = ones(size(VentBinMask));
end

%Dissolved Binning
DissolvedBinMap = GasExchange.BinImages(DissGasRatio, DissolvedThresh);
DissolvedBinMap = DissolvedBinMap.*VentBinMask;%Mask to ventilated volume

%Barrier Binning
BarrierBinMap = GasExchange.BinImages(BarrGasRatio, BarrierThresh);
BarrierBinMap = BarrierBinMap.*VentBinMask;%Mask to ventilated volume

%RBC Binning
RBCBinMap = GasExchange.BinImages(RBCGasRatio, RBCThresh);
RBCBinMap = RBCBinMap.*VentBinMask;%Mask to ventilated volume
RBCBinMask = logical((RBCBinMap-1).*VentBinMask); %Mask to ventilated volume

if VentBinMask == ones(size(VentBinMask))
    %if vent mask is entire FOV, dont mask
    success = 0;%mark as unsuccessful
    RBCBinMask = ones(size(RBCBinMask));
end

%RBC/Barrier Binning
RBCBarrierBinMap = GasExchange.BinImages(RBCBarrRatio, RBCBarrThresh);
RBCBarrierBinMap = RBCBarrierBinMap.*VentBinMask;%Mask to ventilated volume

%RBC Oscillation Binning
RBCOscBinMap = GasExchange.BinImages(RBCOsc, RBCOscThresh);
RBCOscBinMap = RBCOscBinMap.*RBCBinMask;%Mask to RBC-transfer volume
disp('Binning Ratios Completed.')


%% Calculate SNRs
disp('Calculating SNR...')
NoiseMask = imerode(~ProtonMaskRegistered,strel('sphere',7));%avoid edges/artifacts/partialvolume
%((mean signal) - (mean noise)) / (stddev noise)
VentSNR = (mean(abs(VentImage(VentBinMask(:)))) - mean(abs(VentImage(NoiseMask(:)))) ) / std(abs(VentImage(NoiseMask(:))));
GasSNR = (mean(abs(GasImage(VentBinMask(:)))) - mean(abs(GasImage(NoiseMask(:)))) ) / std(abs(GasImage(NoiseMask(:))));
DissolvedSNR = (mean(abs(DissolvedImage(VentBinMask(:)))) - mean(abs(DissolvedImage(NoiseMask(:)))) ) / std(abs(DissolvedImage(NoiseMask(:))));
CorrDissolvedSNR = (mean(abs(CorrDissolvedImage(VentBinMask(:)))) - mean(abs(CorrDissolvedImage(NoiseMask(:)))) ) / std(abs(CorrDissolvedImage(NoiseMask(:))));
BarrierSNR = (mean(BarrierImage(VentBinMask(:))) - mean(BarrierImage(NoiseMask(:))) ) / std(BarrierImage(NoiseMask(:)));
RBCSNR = (mean(RBCImage(VentBinMask(:))) - mean(RBCImage(NoiseMask(:))) ) / std(RBCImage(NoiseMask(:)));
disp('Calculating SNR Completed.')


%% Calculate Histograms
disp('Calculating Histograms...')
%Calculate bin ranges
VentEdges = [-0.5, linspace(0,1,100) 1.5];
BarEdges = [-100, linspace(0,BarrierThresh(3)*3,100) 100];
DissEdges = [-100, linspace(0,DissolvedThresh(3)*3,100) 100];
RBCEdges = [-100, linspace(0,RBCThresh(3)*3,100) 100];
RBCBarEdges = [-100, linspace(0,RBCBarrThresh(3)*3,100) 100];
RBCOscEdges = [-1000, linspace(-100,100,100) 1000];

%Plot Histograms
if HealthyDistPresent
    VentHistFig = GasExchange.CalculateDissolvedHistogram(ScaledVentImage(ProtonMaskRegistered(:)),VentEdges,VentThresh,SixBinMap,HealthyData.VentEdges,HealthyData.HealthyVentFit);
    DissHistFig = GasExchange.CalculateDissolvedHistogram(DissGasRatio(VentBinMask(:)),DissEdges,DissolvedThresh,SixBinMap,HealthyData.DissolvedEdges,HealthyData.HealthyDissFit);
    BarHistFig = GasExchange.CalculateDissolvedHistogram(BarrGasRatio(VentBinMask(:)),BarEdges,BarrierThresh,EightBinMap,HealthyData.BarsEdges,HealthyData.HealthyBarsFit);
    RBCHistFig = GasExchange.CalculateDissolvedHistogram(RBCGasRatio(VentBinMask(:)),RBCEdges,RBCThresh,SixBinMap,HealthyData.RBCEdges,HealthyData.HealthyRBCsFit);
    RBCBarHistFig = GasExchange.CalculateDissolvedHistogram(RBCBarrRatio(VentBinMask(:)),RBCBarEdges,RBCBarrThresh,SixBinRBCBarMap,HealthyData.RBCBarEdges,HealthyData.HealthyRBCBarFit);
    RBCOscHistFig = GasExchange.CalculateDissolvedHistogram(RBCOsc(RBCBinMask(:)),RBCOscEdges,RBCOscThresh,EightBinMap,HealthyData.RBCOscEdges,HealthyData.HealthyRBCOscFit);
else
    VentHistFig = GasExchange.CalculateDissolvedHistogram(ScaledVentImage(ProtonMaskRegistered(:)),VentEdges,VentThresh,SixBinMap,[],[]);
    DissHistFig = GasExchange.CalculateDissolvedHistogram(DissGasRatio(VentBinMask(:)),DissEdges,DissolvedThresh,SixBinMap,[],[]);
    BarHistFig = GasExchange.CalculateDissolvedHistogram(BarrGasRatio(VentBinMask(:)),BarEdges,BarrierThresh,EightBinMap,[],[]);
    RBCHistFig = GasExchange.CalculateDissolvedHistogram(RBCGasRatio(VentBinMask(:)),RBCEdges,RBCThresh,SixBinMap,[],[]);
    RBCBarHistFig = GasExchange.CalculateDissolvedHistogram(RBCBarrRatio(VentBinMask(:)),RBCBarEdges,RBCBarrThresh,SixBinRBCBarMap,[],[]);
    RBCOscHistFig = GasExchange.CalculateDissolvedHistogram(RBCOsc(RBCBinMask(:)),RBCOscEdges,RBCOscThresh,EightBinMap,[],[]);
end

%Edit Names/Titles
set(VentHistFig,'Name','Ventilation Histogram');title(VentHistFig.CurrentAxes,'99 Percentile Scaled Ventilation Histogram');ylabel(VentHistFig.CurrentAxes,'Percentage of Voxels');xlabel(VentHistFig.CurrentAxes,'Scaled Ventilation (a.u.)');
set(DissHistFig,'Name','Dissolved Histogram');title(DissHistFig.CurrentAxes,'Dissolved:Gas Histogram');ylabel(DissHistFig.CurrentAxes,'Percentage of Voxels');xlabel(DissHistFig.CurrentAxes,'Dissolved/Gas Ratio');
set(BarHistFig,'Name','Membrane Histogram');title(BarHistFig.CurrentAxes,'Membrane:Gas Histogram');ylabel(BarHistFig.CurrentAxes,'Percentage of Voxels');xlabel(BarHistFig.CurrentAxes,'Membrane/Gas Ratio (Membrane-uptake)');
set(RBCHistFig,'Name','RBC Histogram');title(RBCHistFig.CurrentAxes,'RBC:Gas Histogram');ylabel(RBCHistFig.CurrentAxes,'Percentage of Voxels');xlabel(RBCHistFig.CurrentAxes,'RBC/Gas Ratio (RBC-transfer)');
set(RBCBarHistFig,'Name','RBC:Membrane Histogram');title(RBCBarHistFig.CurrentAxes,'RBC:Membrane Histogram');ylabel(RBCBarHistFig.CurrentAxes,'Percentage of Voxels');xlabel(RBCBarHistFig.CurrentAxes,'RBC/Membrane Ratio');
set(RBCOscHistFig,'Name','RBC Oscillation Histogram');title(RBCOscHistFig.CurrentAxes,'RBC Oscillation Histogram');ylabel(RBCOscHistFig.CurrentAxes,'Percentage of Voxels');xlabel(RBCOscHistFig.CurrentAxes,'RBC oscilations (% of RBC Signal)');

disp('Calculating Histograms Completed.')


%% Quantify Defects, Intensities, etc.
disp('Quantifying Images...')
%Ventilation Quantification
VentMean = mean(ScaledVentImage(ProtonMaskRegistered(:)));
VentStd = std(ScaledVentImage(ProtonMaskRegistered(:)));
VentBinPercents = zeros(6,1);
for bin = 1:6
    VentBinPercents(bin) = sum(VentBinMap(:)==bin)/sum(ProtonMaskRegistered(:)==1)*100;
end

%Dissolved/Gas Quantification
DissolvedMean = mean(DissGasRatio(VentBinMask(:)));
DissolvedStd = std(DissGasRatio(VentBinMask(:)));
DissolvedBinPercents = zeros(6,1);
for bin = 1:6
    DissolvedBinPercents(bin) = sum(DissolvedBinMap(:)==bin)/sum(VentBinMask(:)==1)*100;
end


%Barrier(Membrane)/Gas Quantification
BarrierUptakeMean = mean(BarrGasRatio(VentBinMask(:)));
BarrierUptakeStd = std(BarrGasRatio(VentBinMask(:)));
BarrierUptakeBinPercents = zeros(8,1);
for bin = 1:8
    BarrierUptakeBinPercents(bin) = sum(BarrierBinMap(:)==bin)/sum(VentBinMask(:)==1)*100;
end

%RBC/Gas Quantification
RBCTransferMean = mean(RBCGasRatio(VentBinMask(:)));
RBCTransferStd = std(RBCGasRatio(VentBinMask(:)));
RBCTransferBinPercents = zeros(6,1);
for bin = 1:6
    RBCTransferBinPercents(bin) = sum(RBCBinMap(:)==bin)/sum(VentBinMask(:)==1)*100;
end

%RBC/Barrier Quantification
RBCBarrierMean = nanmean(RBCBarrRatio(VentBinMask(:)&~isinf(RBCBarrRatio(:))));
RBCBarrierStd = nanstd(RBCBarrRatio(VentBinMask(:)&~isinf(RBCBarrRatio(:))));
RBCBarrierBinPercents = zeros(6,1);
for bin = 1:6
    RBCBarrierBinPercents(bin) = sum(RBCBarrierBinMap(:)==bin)/sum(VentBinMask(:)==1)*100;
end

%RBC Oscillation Quantification
%RBCOscBin1Percent
RBCOscMean = mean(RBCOsc(RBCBinMask(:)));
RBCOscStd = std(RBCOsc(RBCBinMask(:)));
RBCOscBinPercents = zeros(8,1);
for bin = 1:8
    RBCOscBinPercents(bin) = sum(RBCOscBinMap(:)==bin)/sum(RBCBinMask(:)==1)*100;
end
disp('Quantifying Images Completed.')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Export Data and Results                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Summary Figures
%Vent
SumVentFig = figure('Name','Ventitlation Binned','units','normalized','outerposition',[0 0 1 4/NumPlotSlices]);set(SumVentFig,'WindowState','minimized');
set(SumVentFig,'color','white','Units','inches','Position',[0.25 0.25 2*NumPlotSlices 2*2])
tiledlayout(2,NumPlotSlices,'TileSpacing','none','Padding','compact');
for slice=1:NumPlotSlices %Plot Gas Phase Mask Bins slices
    nexttile
    [~,~] = GasExchange_imoverlay.imoverlay(squeeze(abs(ProtonImageRegistered(:,:,Slices_Co(slice),1))),squeeze(VentBinMap(:,:,Slices_Co(slice),1)),[1,6],[0,0.99*ProtonMax],SixBinMap,1,gca);
    colormap(gca,SixBinMap)
    if slice == 1
        colorbar(gca,'north','Ticks',[])
    end
    if slice == round(NumPlotSlices/2)
        title('Ventilation Binned','FontSize',24)
    end
end
for slice=1:NumPlotSlices %Plot Gas Phase Mask Bins slices
    nexttile
    [~,~] = GasExchange_imoverlay.imoverlay(fliplr(rot90(squeeze(abs(ProtonImageRegistered(Slices_Ax(slice),:,:,1))),-1)),fliplr(rot90(squeeze(VentBinMap(Slices_Ax(slice),:,:,1)),-1)),[1,6],[0,0.99*ProtonMax],SixBinMap,1,gca);
    colormap(gca,SixBinMap)
    if slice == 1
        colorbar(gca,'north','Ticks',[])
    end
end
set(SumVentFig,'WindowState','minimized');

%Dissolved
SumDissFig = figure('Name','Dissolved Binned','units','normalized','outerposition',[0 0 1 4/NumPlotSlices]);set(SumDissFig,'WindowState','minimized');
set(SumDissFig,'color','white','Units','inches','Position',[0.25 0.25 2*NumPlotSlices 2*2])
tiledlayout(2,NumPlotSlices,'TileSpacing','none','Padding','compact');
for slice=1:NumPlotSlices %Plot Dissolved Phase Mask Bins slices
    nexttile
    [~,~] = GasExchange_imoverlay.imoverlay(squeeze(abs(ProtonImageRegistered(:,:,Slices_Co(slice),1))),squeeze(DissolvedBinMap(:,:,Slices_Co(slice),1)),[1,6],[0,0.99*ProtonMax],SixBinMap,1,gca);
    colormap(gca,SixBinMap)
    if slice == 1
        colorbar(gca,'north','Ticks',[])
    end
    if slice == round(NumPlotSlices/2)
        title('Dissolved Binned','FontSize',24)
    end
end
for slice=1:NumPlotSlices %Plot Dissolved Phase Mask Bins slices
    nexttile
    [~,~] = GasExchange_imoverlay.imoverlay(fliplr(rot90(squeeze(abs(ProtonImageRegistered(Slices_Ax(slice),:,:,1))),-1)),fliplr(rot90(squeeze(DissolvedBinMap(Slices_Ax(slice),:,:,1)),-1)),[1,6],[0,0.99*ProtonMax],SixBinMap,1,gca);
    colormap(gca,SixBinMap)
    if slice == 1
        colorbar(gca,'north','Ticks',[])
    end
end
set(SumDissFig,'WindowState','minimized');

%Barrier
SumBarrFig = figure('Name','Membrane:Gas Binned','units','normalized','outerposition',[0 0 1 4/NumPlotSlices]);set(SumBarrFig,'WindowState','minimized');
set(SumBarrFig,'color','white','Units','inches','Position',[0.25 0.25 2*NumPlotSlices 2*2])
tiledlayout(2,NumPlotSlices,'TileSpacing','none','Padding','compact');
for slice=1:NumPlotSlices %Plot Barrier Phase Mask Bins slices
    nexttile
    [~,~] = GasExchange_imoverlay.imoverlay(squeeze(abs(ProtonImageRegistered(:,:,Slices_Co(slice),1))),squeeze(BarrierBinMap(:,:,Slices_Co(slice),1)),[1,8],[0,0.99*ProtonMax],EightBinMap,1,gca);
    colormap(gca,EightBinMap)
    if slice == 1
        colorbar(gca,'north','Ticks',[])
    end
    if slice == round(NumPlotSlices/2)
        title('Membrane:Gas Binned','FontSize',24)
    end
end
for slice=1:NumPlotSlices %Plot Barrier Phase Mask Bins slices
    nexttile
    [~,~] = GasExchange_imoverlay.imoverlay(fliplr(rot90(squeeze(abs(ProtonImageRegistered(Slices_Ax(slice),:,:,1))),-1)),fliplr(rot90(squeeze(BarrierBinMap(Slices_Ax(slice),:,:,1)),-1)),[1,8],[0,0.99*ProtonMax],EightBinMap,1,gca);
    colormap(gca,EightBinMap)
    if slice == 1
        colorbar(gca,'north','Ticks',[])
    end
end
set(SumBarrFig,'WindowState','minimized');

%RBC
SumRBCFig = figure('Name','RBC:Gas Binned','units','normalized','outerposition',[0 0 1 4/NumPlotSlices]);set(SumRBCFig,'WindowState','minimized');
set(SumRBCFig,'color','white','Units','inches','Position',[0.25 0.25 2*NumPlotSlices 2*2])
tiledlayout(2,NumPlotSlices,'TileSpacing','none','Padding','compact');
for slice=1:NumPlotSlices %Plot RBC Phase Mask Bins slices
    nexttile
    [~,~] = GasExchange_imoverlay.imoverlay(squeeze(abs(ProtonImageRegistered(:,:,Slices_Co(slice),1))),squeeze(RBCBinMap(:,:,Slices_Co(slice),1)),[1,6],[0,0.99*ProtonMax],SixBinMap,1,gca);
    colormap(gca,SixBinMap)
    if slice == 1
        colorbar(gca,'north','Ticks',[])
    end
    if slice == round(NumPlotSlices/2)
        title('RBC:Gas Binned','FontSize',24)
    end
end
for slice=1:NumPlotSlices %Plot RBC Phase Mask Bins slices
    nexttile
    [~,~] = GasExchange_imoverlay.imoverlay(fliplr(rot90(squeeze(abs(ProtonImageRegistered(Slices_Ax(slice),:,:,1))),-1)),fliplr(rot90(squeeze(RBCBinMap(Slices_Ax(slice),:,:,1)),-1)),[1,6],[0,0.99*ProtonMax],SixBinMap,1,gca);
    colormap(gca,SixBinMap)
    if slice == 1
        colorbar(gca,'north','Ticks',[])
    end
end
set(SumRBCFig,'WindowState','minimized');

%RBC:Barrier
SumRBCBarFig = figure('Name','RBC:Membrane Ratio Binned','units','normalized','outerposition',[0 0 1 4/NumPlotSlices]);set(SumRBCBarFig,'WindowState','minimized');
set(SumRBCBarFig,'color','white','Units','inches','Position',[0.25 0.25 2*NumPlotSlices 2*2])
tiledlayout(2,NumPlotSlices,'TileSpacing','none','Padding','compact');
for slice=1:NumPlotSlices %Plot RBC Barrier Mask Bins slices
    nexttile
    [~,~] = GasExchange_imoverlay.imoverlay(squeeze(abs(ProtonImageRegistered(:,:,Slices_Co(slice),1))),squeeze(RBCBarrierBinMap(:,:,Slices_Co(slice),1)),[1,6],[0,0.99*ProtonMax],SixBinRBCBarMap,1,gca);
    colormap(gca,SixBinRBCBarMap)
    if slice == 1
        colorbar(gca,'north','Ticks',[])
    end
    if slice == round(NumPlotSlices/2)
        title('RBC:Membrane Binned','FontSize',24)
    end
end
for slice=1:NumPlotSlices %Plot RBC Barrier Mask Bins slices
    nexttile
    [~,~] = GasExchange_imoverlay.imoverlay(fliplr(rot90(squeeze(abs(ProtonImageRegistered(Slices_Ax(slice),:,:,1))),-1)),fliplr(rot90(squeeze(RBCBarrierBinMap(Slices_Ax(slice),:,:,1)),-1)),[1,6],[0,0.99*ProtonMax],SixBinRBCBarMap,1,gca);
    colormap(gca,SixBinRBCBarMap)
    if slice == 1
        colorbar(gca,'north','Ticks',[])
    end
end
set(SumRBCBarFig,'WindowState','minimized');

%RBC Oscillations
SumRBCOscFig = figure('Name','RBC Oscillation Percentage Binned','units','normalized','outerposition',[0 0 1 4/NumPlotSlices]);set(SumRBCOscFig,'WindowState','minimized');
set(SumRBCOscFig,'color','white','Units','inches','Position',[0.25 0.25 2*NumPlotSlices 2*2])
tiledlayout(2,NumPlotSlices,'TileSpacing','none','Padding','compact');
for slice=1:NumPlotSlices %Plot RBC Osc Mask Bins slices
    nexttile
    [~,~] = GasExchange_imoverlay.imoverlay(squeeze(abs(ProtonImageRegistered(:,:,Slices_Co(slice),1))),squeeze(RBCOscBinMap(:,:,Slices_Co(slice),1)),[1,8],[0,0.99*ProtonMax],EightBinMap,1,gca);
    colormap(gca,EightBinMap)
    if slice == 1
        colorbar(gca,'north','Ticks',[])
    end
    if slice == round(NumPlotSlices/2)
        title('RBC Oscillation % Binned','FontSize',24)
    end
end
for slice=1:NumPlotSlices %Plot RBC Osc Mask Bins slices
    nexttile
    [~,~] = GasExchange_imoverlay.imoverlay(fliplr(rot90(squeeze(abs(ProtonImageRegistered(Slices_Ax(slice),:,:,1))),-1)),fliplr(rot90(squeeze(RBCOscBinMap(Slices_Ax(slice),:,:,1)),-1)),[1,8],[0,0.99*ProtonMax],EightBinMap,1,gca);
    colormap(gca,EightBinMap)
    if slice == 1
        colorbar(gca,'north','Ticks',[])
    end
end
set(SumRBCOscFig,'WindowState','minimized');
pause(0.001)


%% Save Images
disp('Saving Niftis...')
%Niftis (Intensity Images)
niftiwrite(abs(fliplr(rot90(ProtonImage,-1))),[DataLocation,'\ProtonImage.nii'],'Compressed',true);
info = niftiinfo([DataLocation,'\ProtonImage.nii.gz']);
info.Description = strcat('Package Version: ', ReconVersion);
niftiwrite(abs(fliplr(rot90(ProtonImage,-1))),[DataLocation,'\ProtonImage.nii'],info,'Compressed',true);

niftiwrite(abs(fliplr(rot90(ProtonImageRegistered,-1))),[DataLocation,'\ProtonImageRegistered.nii'],'Compressed',true);
info = niftiinfo([DataLocation,'\ProtonImageRegistered.nii.gz']);
info.Description = strcat('Package Version: ', ReconVersion);
niftiwrite(abs(fliplr(rot90(ProtonImageRegistered,-1))),[DataLocation,'\ProtonImageRegistered.nii'],info,'Compressed',true);

niftiwrite(abs(fliplr(rot90(LungMask,-1))),[DataLocation,'\ProtonMask.nii'],'Compressed',true);
info = niftiinfo([DataLocation,'\ProtonMask.nii.gz']);
info.Description = strcat('Package Version: ', ReconVersion);
niftiwrite(abs(fliplr(rot90(LungMask,-1))),[DataLocation,'\ProtonMask.nii'],info,'Compressed',true);

niftiwrite(abs(fliplr(rot90(ProtonMaskRegistered,-1))),[DataLocation,'\ProtonMaskRegistered.nii'],'Compressed',true);%for obtaining other variables later
info = niftiinfo([DataLocation,'\ProtonMaskRegistered.nii.gz']);
info.Description = strcat('Package Version: ', ReconVersion);
niftiwrite(abs(fliplr(rot90(ProtonMaskRegistered,-1))),[DataLocation,'\ProtonMaskRegistered.nii'],info,'Compressed',true);

niftiwrite(abs(fliplr(rot90(VentImage,-1))),[DataLocation,'\VentImage.nii'],'Compressed',true);
info = niftiinfo([DataLocation,'\VentImage.nii.gz']);
info.Description = strcat('Package Version: ', ReconVersion);
niftiwrite(abs(fliplr(rot90(VentImage,-1))),[DataLocation,'\VentImage.nii'],info,'Compressed',true);

niftiwrite(abs(fliplr(rot90(GasImage,-1))),[DataLocation,'\GasImage.nii'],'Compressed',true);
info = niftiinfo([DataLocation,'\GasImage.nii.gz']);
info.Description = strcat('Package Version: ', ReconVersion);
niftiwrite(abs(fliplr(rot90(GasImage,-1))),[DataLocation,'\GasImage.nii'],info,'Compressed',true);

niftiwrite(abs(fliplr(rot90(CorrDissolvedImage,-1))),[DataLocation,'\DissolvedImage.nii'],'Compressed',true);
info = niftiinfo([DataLocation,'\DissolvedImage.nii.gz']);
info.Description = strcat('Package Version: ', ReconVersion);
niftiwrite(abs(fliplr(rot90(CorrDissolvedImage,-1))),[DataLocation,'\DissolvedImage.nii'],info,'Compressed',true);

niftiwrite(abs(fliplr(rot90(BarrierImage,-1))),[DataLocation,'\BarrierImage.nii'],'Compressed',true);
info = niftiinfo([DataLocation,'\BarrierImage.nii.gz']);
info.Description = strcat('Package Version: ', ReconVersion);
niftiwrite(abs(fliplr(rot90(BarrierImage,-1))),[DataLocation,'\BarrierImage.nii'],info,'Compressed',true);

niftiwrite(abs(fliplr(rot90(RBCImage,-1))),[DataLocation,'\RBCImage.nii'],'Compressed',true);
info = niftiinfo([DataLocation,'\RBCImage.nii.gz']);
info.Description = strcat('Package Version: ', ReconVersion);
niftiwrite(abs(fliplr(rot90(RBCImage,-1))),[DataLocation,'\RBCImage.nii'],info,'Compressed',true);

disp('Saving Niftis Completed.')

%Tiffs (Binned Images)
tiff = figure('MenuBar','none','ToolBar','none','DockControls','off','Resize','off','WindowState','minimized');%figure for tiffs
ax1 = axes('Parent',tiff);ax2 = axes('Parent',tiff);%make axis for both images
set(ax1,'Visible','off');set(ax2,'Visible','off');%turn off axis
set(ax1,'units','inches');set(ax2,'units','inches');%make axis units inches
set(ax1,'position',[0 0 2 2]);set(ax2,'position',[0 0 2 2]);%make axis same as image
set(gcf,'units','inches'); % set the figure units to pixels
set(gcf,'position',[1 1 2 2])% set the position of the figure to axes
disp('Saving Vent Tiff...')

%Vent Binned
for slice=1:H_RecMatrix %repeat for rest of slices
    Xe = imshow(VentBinMap(:,:,slice),[1,6],'Parent',ax2);%plot Xe image
    set(Xe, 'AlphaData', ProtonMaskRegistered(:,:,slice));%mask Xe image
    colormap(ax2,SixBinMap);%change colors
    imshow(abs(ProtonImageRegistered(:,:,slice)),[0,0.99*ProtonMax],'Parent',ax1);%plot H
    X = print('-RGBImage',['-r',num2str(H_RecMatrix/2)]);%2 inches
    if (slice == 1)
        imwrite(X,[DataLocation,'\BinnedVent.tif'],'Description',strcat('Package Version: ', ReconVersion,'; Cohort: ', HealthyCohortNum));%write new/ overwrite tiff
    else
        imwrite(X,[DataLocation,'\BinnedVent.tif'],'WriteMode','append','Description',strcat('Package Version: ', ReconVersion,'; Cohort: ', HealthyCohortNum));%append tiff
    end
end
disp('Saving Vent Tiff Completed. Saving Dissolved Tiff...')

%Dissolved Binned
for slice=1:H_RecMatrix %repeat for rest of slices
    Xe = imshow(DissolvedBinMap(:,:,slice),[1,6],'Parent',ax2);%plot Xe image
    set(Xe, 'AlphaData', VentBinMask(:,:,slice));%mask Xe image
    colormap(ax2,SixBinMap);%change colors
    imshow(abs(ProtonImageRegistered(:,:,slice)),[0,0.99*ProtonMax],'Parent',ax1);%plot H
    X = print('-RGBImage',['-r',num2str(H_RecMatrix/2)]);%2 inches
    if (slice == 1)
        imwrite(X,[DataLocation,'\BinnedDissolved.tif'],'Description',strcat('Package Version: ', ReconVersion,'; Cohort: ', HealthyCohortNum));%write new/ overwrite tiff
    else
        imwrite(X,[DataLocation,'\BinnedDissolved.tif'],'WriteMode','append','Description',strcat('Package Version: ', ReconVersion,'; Cohort: ', HealthyCohortNum));%append tiff
    end
end
disp('Saving Dissolved Tiff Completed. Saving Barrier Tiff...')

%Barrier Binned
for slice=1:H_RecMatrix %repeat for rest of slices
    Xe = imshow(BarrierBinMap(:,:,slice),[1,8],'Parent',ax2);%plot Xe image
    set(Xe, 'AlphaData', VentBinMask(:,:,slice));%mask Xe image
    colormap(ax2,EightBinMap);%change colors
    imshow(abs(ProtonImageRegistered(:,:,slice)),[0,0.99*ProtonMax],'Parent',ax1);%plot H
    X = print('-RGBImage',['-r',num2str(H_RecMatrix/2)]);%2 inches
    if (slice == 1)
        imwrite(X,[DataLocation,'\BinnedBarrierUptake.tif'],'Description',strcat('Package Version: ', ReconVersion,'; Cohort: ', HealthyCohortNum));%write new/ overwrite tiff
    else
        imwrite(X,[DataLocation,'\BinnedBarrierUptake.tif'],'WriteMode','append','Description',strcat('Package Version: ', ReconVersion,'; Cohort: ', HealthyCohortNum));%append tiff
    end
end
disp('Saving Barrier Tiff Completed. Saving RBC Tiff...')

%RBC Binned
for slice=1:H_RecMatrix %repeat for rest of slices
    Xe = imshow(RBCBinMap(:,:,slice),[1,6],'Parent',ax2);%plot Xe image
    set(Xe, 'AlphaData', VentBinMask(:,:,slice));%mask Xe image
    colormap(ax2,SixBinMap);%change colors
    imshow(abs(ProtonImageRegistered(:,:,slice)),[0,0.99*ProtonMax],'Parent',ax1);%plot H
    X = print('-RGBImage',['-r',num2str(H_RecMatrix/2)]);%2 inches
    if (slice == 1)
        imwrite(X,[DataLocation,'\BinnedRBCTransfer.tif'],'Description',strcat('Package Version: ', ReconVersion,'; Cohort: ', HealthyCohortNum));%write new/ overwrite tiff
    else
        imwrite(X,[DataLocation,'\BinnedRBCTransfer.tif'],'WriteMode','append','Description',strcat('Package Version: ', ReconVersion,'; Cohort: ', HealthyCohortNum));%append tiff
    end
end
disp('Saving RBC Tiff Completed. Saving RBC:Barrier Tiff...')

%RBC:Barrier Binned
for slice=1:H_RecMatrix %repeat for rest of slices
    Xe = imshow(RBCBarrierBinMap(:,:,slice),[1,6],'Parent',ax2);%plot Xe image
    set(Xe, 'AlphaData', VentBinMask(:,:,slice));%mask Xe image
    colormap(ax2,SixBinRBCBarMap);%change colors
    imshow(abs(ProtonImageRegistered(:,:,slice)),[0,0.99*ProtonMax],'Parent',ax1);%plot H
    X = print('-RGBImage',['-r',num2str(H_RecMatrix/2)]);%2 inches
    if (slice == 1)
        imwrite(X,[DataLocation,'\BinnedRBCBarrier.tif'],'Description',strcat('Package Version: ', ReconVersion,'; Cohort: ', HealthyCohortNum));%write new/ overwrite tiff
    else
        imwrite(X,[DataLocation,'\BinnedRBCBarrier.tif'],'WriteMode','append','Description',strcat('Package Version: ', ReconVersion,'; Cohort: ', HealthyCohortNum));%append tiff
    end
end
disp('Saving RBC:Barrier Tiff Completed. Saving RBC Oscillation Tiff...')

%RBC Oscillation Binned
for slice=1:H_RecMatrix %repeat for rest of slices
    Xe = imshow(RBCOscBinMap(:,:,slice),[1,8],'Parent',ax2);%plot Xe image
    set(Xe, 'AlphaData', RBCBinMask(:,:,slice));%mask Xe image
    colormap(ax2,EightBinMap);%change colors
    imshow(abs(ProtonImageRegistered(:,:,slice)),[0,0.99*ProtonMax],'Parent',ax1);%plot H
    X = print('-RGBImage',['-r',num2str(H_RecMatrix/2)]);%2 inches
    if (slice == 1)
        imwrite(X,[DataLocation,'\BinnedRBCOscillation.tif'],'Description',strcat('Package Version: ', ReconVersion,'; Cohort: ', HealthyCohortNum));%write new/ overwrite tiff
    else
        imwrite(X,[DataLocation,'\BinnedRBCOscillation.tif'],'WriteMode','append','Description',strcat('Package Version: ', ReconVersion,'; Cohort: ', HealthyCohortNum));%append tiff
    end
end
disp('Saving RBC Oscillation Tiff Completed.')
close(tiff)


%% Save Report
disp('Exporting Results to Powerpoint Summary...')
%Start new presentation
isOpen  = GasExchange_exportToPPTX.exportToPPTX();
if ~isempty(isOpen) %If PowerPoint already started, then close first and then open a new one
    GasExchange_exportToPPTX.exportToPPTX('close');
end
ReportTitle = [Subject,'_',ScanDate,'_',Notes,'_GasExchangeSummary'];
GasExchange_exportToPPTX.exportToPPTX('new','Dimensions',[16 9], ...
    'Title',ReportTitle, ...
    'Author','CPIR @ CCHMC');

%Add slides
GasExchange_exportToPPTX.exportToPPTX('addslide'); %Title Slide
GasExchange_exportToPPTX.exportToPPTX('addtext',sprintf(Subject), ...
    'Position',[0 0 16 4.5], ...
    'VerticalAlignment','bottom', ...
    'HorizontalAlignment','center', ...
    'FontSize',72, ...
    'FontWeight','bold');
GasExchange_exportToPPTX.exportToPPTX('addtext',sprintf(['Scan Date: ',ScanDate,'\nProcessing Date: ',datestr(date,29)]), ...
    'Position',[0 4.5 16 4.5], ...
    'VerticalAlignment','top', ...
    'HorizontalAlignment','center', ...
    'FontSize',36);

%-------------------------------------------------------------CBM WIP
[snr_kspace_beginning,snr_kspace_end,snr_spectrum_end] = GasExchange_DissolvedSNR(CorrectedDissKSpace_SS,PostDissolvedFID);
Fig1_inputs = struct('AppendedDissolvedNMRFit',AppendedDissolvedNMRFit,...
    'AppendedDissolvedFit',AppendedDissolvedFit,...
    'DissolvedKSpace',DissolvedKSpace,...
    'CorrectedDissKSpace_SS',CorrectedDissKSpace_SS,...
    'GasKSpace',GasKSpace,...
    'HPhaseDynamics',HPhaseDynamics,...
    'SS_ind',SS_ind, ...
    'HPulses',HPulses, ...
    'motion_ind',motion_ind,...
    'HDynFit',HDynFit, ...
    'HDynRMSE',HDynRMSE, ...
    'snr_kspace_beginning',snr_kspace_beginning, ...
    'snr_kspace_end',snr_kspace_end,...
    'snr_spectrum_end',snr_spectrum_end...
    );
DissolvedNMR_SigDynamics_fhandle = GasExchange_DissolvedNMR_SigDynamics_fig(Fig1_inputs);  %Signal Dynamics and Dissolved NMR Spectrum
%----------------------------------------------------------------------------------------
GasExchange_exportToPPTX.exportToPPTX('addslide');
GasExchange_exportToPPTX.exportToPPTX('addpicture',DissolvedNMR_SigDynamics_fhandle,'Scale','maxfixed');
GasExchange_exportToPPTX.exportToPPTX('addtext',sprintf('Signal Dynamics and Dissolved-Phase Spectrum'),'Position',[0 0 5 3]);


GasExchange_exportToPPTX.exportToPPTX('addslide'); %Vent Summary
GasExchange_exportToPPTX.exportToPPTX('addpicture',SumVentFig,'Position',[(16-2*NumPlotSlices)/2 0 2*NumPlotSlices 2*2]);
GasExchange_exportToPPTX.exportToPPTX('addpicture',VentHistFig,'Position',[0 4.5 8 4.5]);
if (HealthyDistPresent)
    [~,pm,~,zm] = ztest(VentMean,HealthyData.HealthyMeans.MeanVent,HealthyData.HealthyMeans.StdVent);
    GasExchange_exportToPPTX.exportToPPTX('addtext',sprintf(['Mean: ',num2str(VentMean,'%1.3f'),' (Healthy = ',num2str(HealthyData.HealthyMeans.MeanVent,'%1.3f'),'±',num2str(HealthyData.HealthyMeans.StdVent,'%1.3f'),'); Z = ',num2str(zm,'%1.3f'),', P = ',num2str(pm,'%1.3f')]),'Position',[8 4.5 8 4.5],'HorizontalAlignment','center');
    p=zeros(6,1);
    z=zeros(6,1);
    for bin = 1:6
        [~,p(bin),~,z(bin)] = ztest(VentBinPercents(bin),HealthyData.BinPercentMeans.Vent(bin),HealthyData.BinPercentStds.Vent(bin));
    end
    BinTable   = { ...
        '',{'Defect','BackgroundColor',SixBinMap(1,:),'FontWeight','bold'},{'Low','BackgroundColor',SixBinMap(2,:),'FontWeight','bold'},{'Healthy','BackgroundColor',SixBinMap(3,:),'FontWeight','bold'},{'Healthy','BackgroundColor',SixBinMap(4,:),'FontWeight','bold'},{'High','BackgroundColor',SixBinMap(5,:),'FontWeight','bold','Color','w'},{'High','BackgroundColor',SixBinMap(6,:),'FontWeight','bold','Color','w'};...
        'Subject','Bin 1','Bin 2','Bin 3','Bin 4','Bin 5','Bin 6';...
        Subject,num2str(VentBinPercents(1),'%1.2f'),num2str(VentBinPercents(2),'%1.2f'),num2str(VentBinPercents(3),'%1.2f'),num2str(VentBinPercents(4),'%1.2f'),num2str(VentBinPercents(5),'%1.2f'),num2str(VentBinPercents(6),'%1.2f');...
        'Healthy Cohort',[num2str(HealthyData.BinPercentMeans.Vent(1),'%1.2f'),'±',num2str(HealthyData.BinPercentStds.Vent(1),'%1.2f')],...
                         [num2str(HealthyData.BinPercentMeans.Vent(2),'%1.2f'),'±',num2str(HealthyData.BinPercentStds.Vent(2),'%1.2f')],...
                         [num2str(HealthyData.BinPercentMeans.Vent(3),'%1.2f'),'±',num2str(HealthyData.BinPercentStds.Vent(3),'%1.2f')],...
                         [num2str(HealthyData.BinPercentMeans.Vent(4),'%1.2f'),'±',num2str(HealthyData.BinPercentStds.Vent(4),'%1.2f')],...
                         [num2str(HealthyData.BinPercentMeans.Vent(5),'%1.2f'),'±',num2str(HealthyData.BinPercentStds.Vent(5),'%1.2f')],...
                         [num2str(HealthyData.BinPercentMeans.Vent(6),'%1.2f'),'±',num2str(HealthyData.BinPercentStds.Vent(6),'%1.2f')],;...
        'Z',num2str(z(1),'%1.2f'),num2str(z(2),'%1.2f'),num2str(z(3),'%1.2f'),num2str(z(4),'%1.2f'),num2str(z(5),'%1.2f'),num2str(z(6),'%1.2f');...
        'P',num2str(p(1),'%1.2f'),num2str(p(2),'%1.2f'),num2str(p(3),'%1.2f'),num2str(p(4),'%1.2f'),num2str(p(5),'%1.2f'),num2str(p(6),'%1.2f');...
        };
else
    GasExchange_exportToPPTX.exportToPPTX('addtext',sprintf(['Mean: ',num2str(VentMean,'%1.3f')]),'Position',[8 4.5 8 4.5],'HorizontalAlignment','center');
    BinTable   = { ...
    '',{'Defect','BackgroundColor',SixBinMap(1,:),'FontWeight','bold'},{'Low','BackgroundColor',SixBinMap(2,:),'FontWeight','bold'},{'Healthy','BackgroundColor',SixBinMap(3,:),'FontWeight','bold'},{'Healthy','BackgroundColor',SixBinMap(4,:),'FontWeight','bold'},{'High','BackgroundColor',SixBinMap(5,:),'FontWeight','bold','Color','w'},{'High','BackgroundColor',SixBinMap(6,:),'FontWeight','bold','Color','w'};...
    'Subject','Bin 1','Bin 2','Bin 3','Bin 4','Bin 5','Bin 6';...
    Subject,num2str(VentBinPercents(1),'%1.2f'),num2str(VentBinPercents(2),'%1.2f'),num2str(VentBinPercents(3),'%1.2f'),num2str(VentBinPercents(4),'%1.2f'),num2str(VentBinPercents(5),'%1.2f'),num2str(VentBinPercents(6),'%1.2f');...
    };
end
GasExchange_exportToPPTX.exportToPPTX('addtable',BinTable,'Position',[8 5 7.75 2],'Vert','middle','Horiz','center','FontSize',12);
GasExchange_exportToPPTX.exportToPPTX('addtext',sprintf('Ventilation Summary'),'Position',[0 0 5 3]);


GasExchange_exportToPPTX.exportToPPTX('addslide'); %Barrier Summary
GasExchange_exportToPPTX.exportToPPTX('addpicture',SumBarrFig,'Position',[(16-2*NumPlotSlices)/2 0 2*NumPlotSlices 2*2]);
GasExchange_exportToPPTX.exportToPPTX('addpicture',BarHistFig,'Position',[0 4.5 8 4.5]);
if (HealthyDistPresent)
    [~,pm,~,zm] = ztest(BarrierUptakeMean,HealthyData.HealthyMeans.MeanBarrier,HealthyData.HealthyMeans.StdBarrier);
    GasExchange_exportToPPTX.exportToPPTX('addtext',sprintf(['Mean: ',num2str(BarrierUptakeMean,'%1.5f'),' (Healthy = ',num2str(HealthyData.HealthyMeans.MeanBarrier,'%1.5f'),'±',num2str(HealthyData.HealthyMeans.StdBarrier,'%1.5f'),'); Z = ',num2str(zm,'%1.3f'),', P = ',num2str(pm,'%1.3f')]),'Position',[8 4.5 8 4.5],'HorizontalAlignment','center');
    p=zeros(8,1);
    z=zeros(8,1);
    for bin = 1:8
        [~,p(bin),~,z(bin)] = ztest(BarrierUptakeBinPercents(bin),HealthyData.BinPercentMeans.Barrier(bin),HealthyData.BinPercentStds.Barrier(bin));
    end
    BinTable   = { ...
    '',{'Defect','BackgroundColor',EightBinMap(1,:),'FontWeight','bold'},{'Low','BackgroundColor',EightBinMap(2,:),'FontWeight','bold'},{'Healthy','BackgroundColor',EightBinMap(3,:),'FontWeight','bold'},{'Healthy','BackgroundColor',EightBinMap(4,:),'FontWeight','bold'},{'Elevated','BackgroundColor',EightBinMap(5,:),'FontWeight','bold'},{'Elevated','BackgroundColor',EightBinMap(6,:),'FontWeight','bold'},{'High','BackgroundColor',EightBinMap(7,:),'FontWeight','bold','Color','w'},{'High','BackgroundColor',EightBinMap(8,:),'FontWeight','bold','Color','w'};...
    'Subject','Bin 1','Bin 2','Bin 3','Bin 4','Bin 5','Bin 6','Bin 7','Bin 8';...
    Subject,num2str(BarrierUptakeBinPercents(1),'%1.2f'),num2str(BarrierUptakeBinPercents(2),'%1.2f'),num2str(BarrierUptakeBinPercents(3),'%1.2f'),num2str(BarrierUptakeBinPercents(4),'%1.2f'),num2str(BarrierUptakeBinPercents(5),'%1.2f'),num2str(BarrierUptakeBinPercents(6),'%1.2f'),num2str(BarrierUptakeBinPercents(7),'%1.2f'),num2str(BarrierUptakeBinPercents(8),'%1.2f');...
    'Healthy Cohort',[num2str(HealthyData.BinPercentMeans.Barrier(1),'%1.2f'),'±',num2str(HealthyData.BinPercentStds.Barrier(1),'%1.2f')],...
                     [num2str(HealthyData.BinPercentMeans.Barrier(2),'%1.2f'),'±',num2str(HealthyData.BinPercentStds.Barrier(2),'%1.2f')],...
                     [num2str(HealthyData.BinPercentMeans.Barrier(3),'%1.2f'),'±',num2str(HealthyData.BinPercentStds.Barrier(3),'%1.2f')],...
                     [num2str(HealthyData.BinPercentMeans.Barrier(4),'%1.2f'),'±',num2str(HealthyData.BinPercentStds.Barrier(4),'%1.2f')],...
                     [num2str(HealthyData.BinPercentMeans.Barrier(5),'%1.2f'),'±',num2str(HealthyData.BinPercentStds.Barrier(5),'%1.2f')],...
                     [num2str(HealthyData.BinPercentMeans.Barrier(6),'%1.2f'),'±',num2str(HealthyData.BinPercentStds.Barrier(6),'%1.2f')],...
                     [num2str(HealthyData.BinPercentMeans.Barrier(7),'%1.2f'),'±',num2str(HealthyData.BinPercentStds.Barrier(7),'%1.2f')],...
                     [num2str(HealthyData.BinPercentMeans.Barrier(8),'%1.2f'),'±',num2str(HealthyData.BinPercentStds.Barrier(8),'%1.2f')],;...
    'Z',num2str(z(1),'%1.2f'),num2str(z(2),'%1.2f'),num2str(z(3),'%1.2f'),num2str(z(4),'%1.2f'),num2str(z(5),'%1.2f'),num2str(z(6),'%1.2f'),num2str(z(7),'%1.2f'),num2str(z(8),'%1.2f');...
    'P',num2str(p(1),'%1.2f'),num2str(p(2),'%1.2f'),num2str(p(3),'%1.2f'),num2str(p(4),'%1.2f'),num2str(p(5),'%1.2f'),num2str(p(6),'%1.2f'),num2str(p(7),'%1.2f'),num2str(p(8),'%1.2f');...
    };
else
    GasExchange_exportToPPTX.exportToPPTX('addtext',sprintf(['Mean: ',num2str(BarrierUptakeMean,'%1.5f')]),'Position',[8 4.5 8 4.5],'HorizontalAlignment','center');
        BinTable   = { ...
    '',{'Defect','BackgroundColor',EightBinMap(1,:),'FontWeight','bold'},{'Low','BackgroundColor',EightBinMap(2,:),'FontWeight','bold'},{'Healthy','BackgroundColor',EightBinMap(3,:),'FontWeight','bold'},{'Healthy','BackgroundColor',EightBinMap(4,:),'FontWeight','bold'},{'Elevated','BackgroundColor',EightBinMap(5,:),'FontWeight','bold'},{'Elevated','BackgroundColor',EightBinMap(6,:),'FontWeight','bold'},{'High','BackgroundColor',EightBinMap(7,:),'FontWeight','bold','Color','w'},{'High','BackgroundColor',EightBinMap(8,:),'FontWeight','bold','Color','w'};...
    'Subject','Bin 1','Bin 2','Bin 3','Bin 4','Bin 5','Bin 6','Bin 7','Bin 8';...
    Subject,num2str(BarrierUptakeBinPercents(1),'%1.2f'),num2str(BarrierUptakeBinPercents(2),'%1.2f'),num2str(BarrierUptakeBinPercents(3),'%1.2f'),num2str(BarrierUptakeBinPercents(4),'%1.2f'),num2str(BarrierUptakeBinPercents(5),'%1.2f'),num2str(BarrierUptakeBinPercents(6),'%1.2f'),num2str(BarrierUptakeBinPercents(7),'%1.2f'),num2str(BarrierUptakeBinPercents(8),'%1.2f');...
    };
end
GasExchange_exportToPPTX.exportToPPTX('addtable',BinTable,'Position',[8 5 7.75 2],'Vert','middle','Horiz','center','FontSize',12);
GasExchange_exportToPPTX.exportToPPTX('addtext',sprintf('Membrane-uptake Summary'),'Position',[0 0 5 3]);

GasExchange_exportToPPTX.exportToPPTX('addslide'); %RBC Summary
GasExchange_exportToPPTX.exportToPPTX('addpicture',SumRBCFig,'Position',[(16-2*NumPlotSlices)/2 0 2*NumPlotSlices 2*2]);
GasExchange_exportToPPTX.exportToPPTX('addpicture',RBCHistFig,'Position',[0 4.5 8 4.5]);
if (HealthyDistPresent)
    [~,pm,~,zm] = ztest(RBCTransferMean,HealthyData.HealthyMeans.MeanRBC,HealthyData.HealthyMeans.StdRBC);
    GasExchange_exportToPPTX.exportToPPTX('addtext',sprintf(['Mean: ',num2str(RBCTransferMean,'%1.3e'),' (Healthy = ',num2str(HealthyData.HealthyMeans.MeanRBC,'%1.3e'),'±',num2str(HealthyData.HealthyMeans.StdRBC,'%1.3e'),'); Z = ',num2str(zm,'%1.3f'),', P = ',num2str(pm,'%1.3f')]),'Position',[8 4.5 8 4.5],'HorizontalAlignment','center');
    p=zeros(6,1);
    z=zeros(6,1);
    for bin = 1:6
        [~,p(bin),~,z(bin)] = ztest(RBCTransferBinPercents(bin),HealthyData.BinPercentMeans.RBC(bin),HealthyData.BinPercentStds.RBC(bin));
    end
    BinTable   = { ...
    '',{'Defect','BackgroundColor',SixBinMap(1,:),'FontWeight','bold'},{'Low','BackgroundColor',SixBinMap(2,:),'FontWeight','bold'},{'Healthy','BackgroundColor',SixBinMap(3,:),'FontWeight','bold'},{'Healthy','BackgroundColor',SixBinMap(4,:),'FontWeight','bold'},{'High','BackgroundColor',SixBinMap(5,:),'FontWeight','bold','Color','w'},{'High','BackgroundColor',SixBinMap(6,:),'FontWeight','bold','Color','w'};...
    'Subject','Bin 1','Bin 2','Bin 3','Bin 4','Bin 5','Bin 6';...
    Subject,num2str(RBCTransferBinPercents(1),'%1.2f'),num2str(RBCTransferBinPercents(2),'%1.2f'),num2str(RBCTransferBinPercents(3),'%1.2f'),num2str(RBCTransferBinPercents(4),'%1.2f'),num2str(RBCTransferBinPercents(5),'%1.2f'),num2str(RBCTransferBinPercents(6),'%1.2f');...
    'Healthy Cohort',[num2str(HealthyData.BinPercentMeans.RBC(1),'%1.2f'),'±',num2str(HealthyData.BinPercentStds.RBC(1),'%1.2f')],...
                     [num2str(HealthyData.BinPercentMeans.RBC(2),'%1.2f'),'±',num2str(HealthyData.BinPercentStds.RBC(2),'%1.2f')],...
                     [num2str(HealthyData.BinPercentMeans.RBC(3),'%1.2f'),'±',num2str(HealthyData.BinPercentStds.RBC(3),'%1.2f')],...
                     [num2str(HealthyData.BinPercentMeans.RBC(4),'%1.2f'),'±',num2str(HealthyData.BinPercentStds.RBC(4),'%1.2f')],...
                     [num2str(HealthyData.BinPercentMeans.RBC(5),'%1.2f'),'±',num2str(HealthyData.BinPercentStds.RBC(5),'%1.2f')],...
                     [num2str(HealthyData.BinPercentMeans.RBC(6),'%1.2f'),'±',num2str(HealthyData.BinPercentStds.RBC(6),'%1.2f')],;...
    'Z',num2str(z(1),'%1.2f'),num2str(z(2),'%1.2f'),num2str(z(3),'%1.2f'),num2str(z(4),'%1.2f'),num2str(z(5),'%1.2f'),num2str(z(6),'%1.2f');...
    'P',num2str(p(1),'%1.2f'),num2str(p(2),'%1.2f'),num2str(p(3),'%1.2f'),num2str(p(4),'%1.2f'),num2str(p(5),'%1.2f'),num2str(p(6),'%1.2f');...
    };
else
    GasExchange_exportToPPTX.exportToPPTX('addtext',sprintf(['Mean: ',num2str(RBCTransferMean,'%1.3e')]),'Position',[8 4.5 8 4.5],'HorizontalAlignment','center');
    BinTable   = { ...
    '',{'Defect','BackgroundColor',SixBinMap(1,:),'FontWeight','bold'},{'Low','BackgroundColor',SixBinMap(2,:),'FontWeight','bold'},{'Healthy','BackgroundColor',SixBinMap(3,:),'FontWeight','bold'},{'Healthy','BackgroundColor',SixBinMap(4,:),'FontWeight','bold'},{'High','BackgroundColor',SixBinMap(5,:),'FontWeight','bold','Color','w'},{'High','BackgroundColor',SixBinMap(6,:),'FontWeight','bold','Color','w'};...
    'Subject','Bin 1','Bin 2','Bin 3','Bin 4','Bin 5','Bin 6';...
    Subject,num2str(RBCTransferBinPercents(1),'%1.2f'),num2str(RBCTransferBinPercents(2),'%1.2f'),num2str(RBCTransferBinPercents(3),'%1.2f'),num2str(RBCTransferBinPercents(4),'%1.2f'),num2str(RBCTransferBinPercents(5),'%1.2f'),num2str(RBCTransferBinPercents(6),'%1.2f');...
    };
end
GasExchange_exportToPPTX.exportToPPTX('addtable',BinTable,'Position',[8 5 7.75 2],'Vert','middle','Horiz','center','FontSize',12);
GasExchange_exportToPPTX.exportToPPTX('addtext',sprintf('RBC-Transfer Summary'),'Position',[0 0 5 3]);

GasExchange_exportToPPTX.exportToPPTX('addslide'); %RBC:Barrier Summary
GasExchange_exportToPPTX.exportToPPTX('addpicture',SumRBCBarFig,'Position',[(16-2*NumPlotSlices)/2 0 2*NumPlotSlices 2*2]);
GasExchange_exportToPPTX.exportToPPTX('addpicture',RBCBarHistFig,'Position',[0 4.5 8 4.5]);
if (HealthyDistPresent)
    [~,pm,~,zm] = ztest(RBCBarrierMean,HealthyData.HealthyMeans.MeanRBCBar,HealthyData.HealthyMeans.StdRBCBar);
    GasExchange_exportToPPTX.exportToPPTX('addtext',sprintf(['Mean: ',num2str(RBCBarrierMean,'%1.3f'),' (Healthy = ',num2str(HealthyData.HealthyMeans.MeanRBCBar,'%1.3f'),'±',num2str(HealthyData.HealthyMeans.StdRBCBar,'%1.3f'),'); Z = ',num2str(zm,'%1.3f'),', P = ',num2str(pm,'%1.3f')]),'Position',[8 4.5 8 4.5],'HorizontalAlignment','center');
    p=zeros(6,1);
    z=zeros(6,1);
    for bin = 1:6
        [~,p(bin),~,z(bin)] = ztest(RBCBarrierBinPercents(bin),HealthyData.BinPercentMeans.RBCBar(bin),HealthyData.BinPercentStds.RBCBar(bin));
    end
    BinTable   = { ...
    '',{'Low RBC:Membrane','BackgroundColor',SixBinRBCBarMap(1,:),'FontWeight','bold','Color','w'},{'Low RBC:Membrane','BackgroundColor',SixBinRBCBarMap(2,:),'FontWeight','bold','Color','w'},{'Healthy','BackgroundColor',SixBinRBCBarMap(3,:),'FontWeight','bold'},{'Healthy','BackgroundColor',SixBinRBCBarMap(4,:),'FontWeight','bold'},{'Elevated RBC:Membrane','BackgroundColor',SixBinRBCBarMap(5,:),'FontWeight','bold'},{'High RBC:Membrane','BackgroundColor',SixBinRBCBarMap(6,:),'FontWeight','bold'};...
    'Subject','Bin 1','Bin 2','Bin 3','Bin 4','Bin 5','Bin 6';...
    Subject,num2str(RBCBarrierBinPercents(1),'%1.2f'),num2str(RBCBarrierBinPercents(2),'%1.2f'),num2str(RBCBarrierBinPercents(3),'%1.2f'),num2str(RBCBarrierBinPercents(4),'%1.2f'),num2str(RBCBarrierBinPercents(5),'%1.2f'),num2str(RBCBarrierBinPercents(6),'%1.2f');...
    'Healthy Cohort',[num2str(HealthyData.BinPercentMeans.RBCBar(1),'%1.2f'),'±',num2str(HealthyData.BinPercentStds.RBCBar(1),'%1.2f')],...
                     [num2str(HealthyData.BinPercentMeans.RBCBar(2),'%1.2f'),'±',num2str(HealthyData.BinPercentStds.RBCBar(2),'%1.2f')],...
                     [num2str(HealthyData.BinPercentMeans.RBCBar(3),'%1.2f'),'±',num2str(HealthyData.BinPercentStds.RBCBar(3),'%1.2f')],...
                     [num2str(HealthyData.BinPercentMeans.RBCBar(4),'%1.2f'),'±',num2str(HealthyData.BinPercentStds.RBCBar(4),'%1.2f')],...
                     [num2str(HealthyData.BinPercentMeans.RBCBar(5),'%1.2f'),'±',num2str(HealthyData.BinPercentStds.RBCBar(5),'%1.2f')],...
                     [num2str(HealthyData.BinPercentMeans.RBCBar(6),'%1.2f'),'±',num2str(HealthyData.BinPercentStds.RBCBar(6),'%1.2f')],;...
    'Z',num2str(z(1),'%1.2f'),num2str(z(2),'%1.2f'),num2str(z(3),'%1.2f'),num2str(z(4),'%1.2f'),num2str(z(5),'%1.2f'),num2str(z(6),'%1.2f');...
    'P',num2str(p(1),'%1.2f'),num2str(p(2),'%1.2f'),num2str(p(3),'%1.2f'),num2str(p(4),'%1.2f'),num2str(p(5),'%1.2f'),num2str(p(6),'%1.2f');...
    };
else
    GasExchange_exportToPPTX.exportToPPTX('addtext',sprintf(['Mean: ',num2str(RBCBarrierMean,'%1.3f')]),'Position',[8 4.5 8 4.5],'HorizontalAlignment','center');
    BinTable   = { ...
    '',{'Defect','BackgroundColor',SixBinRBCBarMap(1,:),'FontWeight','bold','Color','w'},{'Low','BackgroundColor',SixBinRBCBarMap(2,:),'FontWeight','bold','Color','w'},{'Healthy','BackgroundColor',SixBinRBCBarMap(3,:),'FontWeight','bold'},{'Healthy','BackgroundColor',SixBinRBCBarMap(4,:),'FontWeight','bold'},{'High','BackgroundColor',SixBinRBCBarMap(5,:),'FontWeight','bold'},{'High','BackgroundColor',SixBinRBCBarMap(6,:),'FontWeight','bold'};...
    'Subject','Bin 1','Bin 2','Bin 3','Bin 4','Bin 5','Bin 6';...
    Subject,num2str(RBCBarrierBinPercents(1),'%1.2f'),num2str(RBCBarrierBinPercents(2),'%1.2f'),num2str(RBCBarrierBinPercents(3),'%1.2f'),num2str(RBCBarrierBinPercents(4),'%1.2f'),num2str(RBCBarrierBinPercents(5),'%1.2f'),num2str(RBCBarrierBinPercents(6),'%1.2f');...
    };
end
GasExchange_exportToPPTX.exportToPPTX('addtable',BinTable,'Position',[8 5 7.75 2],'Vert','middle','Horiz','center','FontSize',12);
GasExchange_exportToPPTX.exportToPPTX('addtext',sprintf('RBC:Membrane Summary'),'Position',[0 0 5 3]);

GasExchange_exportToPPTX.exportToPPTX('addslide'); %Dissolved Summary
GasExchange_exportToPPTX.exportToPPTX('addpicture',SumDissFig,'Position',[(16-2*NumPlotSlices)/2 0 2*NumPlotSlices 2*2]);
GasExchange_exportToPPTX.exportToPPTX('addpicture',DissHistFig,'Position',[0 4.5 8 4.5]);
if (HealthyDistPresent)
    [~,pm,~,zm] = ztest(DissolvedMean,HealthyData.HealthyMeans.MeanDissolved,HealthyData.HealthyMeans.StdDissolved);
    GasExchange_exportToPPTX.exportToPPTX('addtext',sprintf(['Mean: ',num2str(DissolvedMean,'%1.5f'),' (Healthy = ',num2str(HealthyData.HealthyMeans.MeanDissolved,'%1.5f'),'±',num2str(HealthyData.HealthyMeans.StdDissolved,'%1.5f'),'); Z = ',num2str(zm,'%1.3f'),', P = ',num2str(pm,'%1.3f')]),'Position',[8 4.5 8 4.5],'HorizontalAlignment','center');
    p=zeros(6,1);
    z=zeros(6,1);
    for bin = 1:6
        [~,p(bin),~,z(bin)] = ztest(DissolvedBinPercents(bin),HealthyData.BinPercentMeans.Dissolved(bin),HealthyData.BinPercentStds.Dissolved(bin));
    end
    BinTable   = { ...
    '',{'Defect','BackgroundColor',SixBinMap(1,:),'FontWeight','bold'},{'Low','BackgroundColor',SixBinMap(2,:),'FontWeight','bold'},{'Healthy','BackgroundColor',SixBinMap(3,:),'FontWeight','bold'},{'Healthy','BackgroundColor',SixBinMap(4,:),'FontWeight','bold'},{'High','BackgroundColor',SixBinMap(5,:),'FontWeight','bold','Color','w'},{'High','BackgroundColor',SixBinMap(6,:),'FontWeight','bold','Color','w'};...
    'Subject','Bin 1','Bin 2','Bin 3','Bin 4','Bin 5','Bin 6';...
    Subject,num2str(DissolvedBinPercents(1),'%1.2f'),num2str(DissolvedBinPercents(2),'%1.2f'),num2str(DissolvedBinPercents(3),'%1.2f'),num2str(DissolvedBinPercents(4),'%1.2f'),num2str(DissolvedBinPercents(5),'%1.2f'),num2str(DissolvedBinPercents(6),'%1.2f');...
    'Healthy Cohort',[num2str(HealthyData.BinPercentMeans.Dissolved(1),'%1.2f'),'±',num2str(HealthyData.BinPercentStds.Dissolved(1),'%1.2f')],...
                     [num2str(HealthyData.BinPercentMeans.Dissolved(2),'%1.2f'),'±',num2str(HealthyData.BinPercentStds.Dissolved(2),'%1.2f')],...
                     [num2str(HealthyData.BinPercentMeans.Dissolved(3),'%1.2f'),'±',num2str(HealthyData.BinPercentStds.Dissolved(3),'%1.2f')],...
                     [num2str(HealthyData.BinPercentMeans.Dissolved(4),'%1.2f'),'±',num2str(HealthyData.BinPercentStds.Dissolved(4),'%1.2f')],...
                     [num2str(HealthyData.BinPercentMeans.Dissolved(5),'%1.2f'),'±',num2str(HealthyData.BinPercentStds.Dissolved(5),'%1.2f')],...
                     [num2str(HealthyData.BinPercentMeans.Dissolved(6),'%1.2f'),'±',num2str(HealthyData.BinPercentStds.Dissolved(6),'%1.2f')],;...
    'Z',num2str(z(1),'%1.2f'),num2str(z(2),'%1.2f'),num2str(z(3),'%1.2f'),num2str(z(4),'%1.2f'),num2str(z(5),'%1.2f'),num2str(z(6),'%1.2f');...
    'P',num2str(p(1),'%1.2f'),num2str(p(2),'%1.2f'),num2str(p(3),'%1.2f'),num2str(p(4),'%1.2f'),num2str(p(5),'%1.2f'),num2str(p(6),'%1.2f');...
    };
else
    GasExchange_exportToPPTX.exportToPPTX('addtext',sprintf(['Mean: ',num2str(DissolvedMean,'%1.5f')]),'Position',[8 4.5 8 4.5],'HorizontalAlignment','center');
    BinTable   = { ...
    '',{'Defect','BackgroundColor',SixBinMap(1,:),'FontWeight','bold'},{'Low','BackgroundColor',SixBinMap(2,:),'FontWeight','bold'},{'Healthy','BackgroundColor',SixBinMap(3,:),'FontWeight','bold'},{'Healthy','BackgroundColor',SixBinMap(4,:),'FontWeight','bold'},{'High','BackgroundColor',SixBinMap(5,:),'FontWeight','bold','Color','w'},{'High','BackgroundColor',SixBinMap(6,:),'FontWeight','bold','Color','w'};...
    'Subject','Bin 1','Bin 2','Bin 3','Bin 4','Bin 5','Bin 6';...
    Subject,num2str(DissolvedBinPercents(1),'%1.2f'),num2str(DissolvedBinPercents(2),'%1.2f'),num2str(DissolvedBinPercents(3),'%1.2f'),num2str(DissolvedBinPercents(4),'%1.2f'),num2str(DissolvedBinPercents(5),'%1.2f'),num2str(DissolvedBinPercents(6),'%1.2f');...
    };
end
GasExchange_exportToPPTX.exportToPPTX('addtable',BinTable,'Position',[8 5 7.75 2],'Vert','middle','Horiz','center','FontSize',12);
GasExchange_exportToPPTX.exportToPPTX('addtext',sprintf('Dissolved Summary'),'Position',[0 0 5 3]);

GasExchange_exportToPPTX.exportToPPTX('addslide'); %Registration
GasExchange_exportToPPTX.exportToPPTX('addpicture',Registration);
GasExchange_exportToPPTX.exportToPPTX('addtext',sprintf('Registration'),'Position',[0 0 5 3]);

GasExchange_exportToPPTX.exportToPPTX('addslide'); %RBC Osc Summary
GasExchange_exportToPPTX.exportToPPTX('addpicture',SumRBCOscFig,'Position',[(16-2*NumPlotSlices)/2 0 2*NumPlotSlices 2*2]);
GasExchange_exportToPPTX.exportToPPTX('addpicture',RBCOscHistFig,'Position',[0 4.5 8 4.5]);
if (HealthyDistPresent)
    [~,pm,~,zm] = ztest(RBCOscMean,HealthyData.HealthyMeans.MeanRBCOsc,HealthyData.HealthyMeans.StdRBCOsc);
    GasExchange_exportToPPTX.exportToPPTX('addtext',sprintf(['Mean: ',num2str(RBCOscMean,'%2.1f'),' (Healthy = ',num2str(HealthyData.HealthyMeans.MeanRBCOsc,'%2.1f'),'±',num2str(HealthyData.HealthyMeans.StdRBCOsc,'%2.1f'),'); Z = ',num2str(zm,'%1.3f'),', P = ',num2str(pm,'%1.3f')]),'Position',[8 4.5 8 4.5],'HorizontalAlignment','center');
    p=zeros(8,1);
    z=zeros(8,1);
    for bin = 1:8
        [~,p(bin),~,z(bin)] = ztest(RBCOscBinPercents(bin),HealthyData.BinPercentMeans.RBCOsc(bin),HealthyData.BinPercentStds.RBCOsc(bin));
    end
    BinTable   = { ...
    '',{'Defect','BackgroundColor',EightBinMap(1,:),'FontWeight','bold'},{'Low','BackgroundColor',EightBinMap(2,:),'FontWeight','bold'},{'Healthy','BackgroundColor',EightBinMap(3,:),'FontWeight','bold'},{'Healthy','BackgroundColor',EightBinMap(4,:),'FontWeight','bold'},{'Elevated','BackgroundColor',EightBinMap(5,:),'FontWeight','bold'},{'Elevated','BackgroundColor',EightBinMap(6,:),'FontWeight','bold'},{'High','BackgroundColor',EightBinMap(7,:),'FontWeight','bold','Color','w'},{'High','BackgroundColor',EightBinMap(8,:),'FontWeight','bold','Color','w'};...
    'Subject','Bin 1','Bin 2','Bin 3','Bin 4','Bin 5','Bin 6','Bin 7','Bin 8';...
    Subject,num2str(RBCOscBinPercents(1),'%1.2f'),num2str(RBCOscBinPercents(2),'%1.2f'),num2str(RBCOscBinPercents(3),'%1.2f'),num2str(RBCOscBinPercents(4),'%1.2f'),num2str(RBCOscBinPercents(5),'%1.2f'),num2str(RBCOscBinPercents(6),'%1.2f'),num2str(RBCOscBinPercents(7),'%1.2f'),num2str(RBCOscBinPercents(8),'%1.2f');...
    'Healthy Cohort',[num2str(HealthyData.BinPercentMeans.RBCOsc(1),'%1.2f'),'±',num2str(HealthyData.BinPercentStds.RBCOsc(1),'%1.2f')],...
                 [num2str(HealthyData.BinPercentMeans.RBCOsc(2),'%1.2f'),'±',num2str(HealthyData.BinPercentStds.RBCOsc(2),'%1.2f')],...
                 [num2str(HealthyData.BinPercentMeans.RBCOsc(3),'%1.2f'),'±',num2str(HealthyData.BinPercentStds.RBCOsc(3),'%1.2f')],...
                 [num2str(HealthyData.BinPercentMeans.RBCOsc(4),'%1.2f'),'±',num2str(HealthyData.BinPercentStds.RBCOsc(4),'%1.2f')],...
                 [num2str(HealthyData.BinPercentMeans.RBCOsc(5),'%1.2f'),'±',num2str(HealthyData.BinPercentStds.RBCOsc(5),'%1.2f')],...
                 [num2str(HealthyData.BinPercentMeans.RBCOsc(6),'%1.2f'),'±',num2str(HealthyData.BinPercentStds.RBCOsc(6),'%1.2f')],...
                 [num2str(HealthyData.BinPercentMeans.RBCOsc(7),'%1.2f'),'±',num2str(HealthyData.BinPercentStds.RBCOsc(7),'%1.2f')],...
                 [num2str(HealthyData.BinPercentMeans.RBCOsc(8),'%1.2f'),'±',num2str(HealthyData.BinPercentStds.RBCOsc(8),'%1.2f')],;...
    'Z',num2str(z(1),'%1.2f'),num2str(z(2),'%1.2f'),num2str(z(3),'%1.2f'),num2str(z(4),'%1.2f'),num2str(z(5),'%1.2f'),num2str(z(6),'%1.2f'),num2str(z(7),'%1.2f'),num2str(z(8),'%1.2f');...
    'P',num2str(p(1),'%1.2f'),num2str(p(2),'%1.2f'),num2str(p(3),'%1.2f'),num2str(p(4),'%1.2f'),num2str(p(5),'%1.2f'),num2str(p(6),'%1.2f'),num2str(p(7),'%1.2f'),num2str(p(8),'%1.2f');...
    };
else
    GasExchange_exportToPPTX.exportToPPTX('addtext',sprintf(['Mean: ',num2str(RBCOscMean,'%2.1f')]),'Position',[8 4.5 8 4.5],'HorizontalAlignment','center');
    BinTable   = { ...
    '',{'Defect','BackgroundColor',EightBinMap(1,:),'FontWeight','bold'},{'Low','BackgroundColor',EightBinMap(2,:),'FontWeight','bold'},{'Healthy','BackgroundColor',EightBinMap(3,:),'FontWeight','bold'},{'Healthy','BackgroundColor',EightBinMap(4,:),'FontWeight','bold'},{'Elevated','BackgroundColor',EightBinMap(5,:),'FontWeight','bold'},{'Elevated','BackgroundColor',EightBinMap(6,:),'FontWeight','bold'},{'High','BackgroundColor',EightBinMap(7,:),'FontWeight','bold','Color','w'},{'High','BackgroundColor',EightBinMap(8,:),'FontWeight','bold','Color','w'};...
    'Subject','Bin 1','Bin 2','Bin 3','Bin 4','Bin 5','Bin 6','Bin 7','Bin 8';...
    Subject,num2str(RBCOscBinPercents(1),'%1.2f'),num2str(RBCOscBinPercents(2),'%1.2f'),num2str(RBCOscBinPercents(3),'%1.2f'),num2str(RBCOscBinPercents(4),'%1.2f'),num2str(RBCOscBinPercents(5),'%1.2f'),num2str(RBCOscBinPercents(6),'%1.2f'),num2str(RBCOscBinPercents(7),'%1.2f'),num2str(RBCOscBinPercents(8),'%1.2f');...
    };
end
GasExchange_exportToPPTX.exportToPPTX('addtable',BinTable,'Position',[8 5 7.75 2],'Vert','middle','Horiz','center','FontSize',12);
GasExchange_exportToPPTX.exportToPPTX('addtext',sprintf('RBC Oscillation Summary'),'Position',[0 0 5 3]);

GasExchange_exportToPPTX.exportToPPTX('addslide'); %RBC Osc Workflow
GasExchange_exportToPPTX.exportToPPTX('addpicture',OscWorkFlow_Fig);
GasExchange_exportToPPTX.exportToPPTX('addtext',sprintf('RBC Oscillation Workflow'),'Position',[0 0 5 3]);


GasExchange_exportToPPTX.exportToPPTX('addslide'); %Proton Image
GasExchange_exportToPPTX.exportToPPTX('addpicture',ProtonMontage);
GasExchange_exportToPPTX.exportToPPTX('addtext',sprintf('Proton Image'),'Position',[0 0 5 3]);

GasExchange_exportToPPTX.exportToPPTX('addslide'); %Proton Mask
GasExchange_exportToPPTX.exportToPPTX('addpicture',ProtonMaskMontage);
GasExchange_exportToPPTX.exportToPPTX('addtext',sprintf('Proton Mask'),'Position',[0 0 5 3]);

GasExchange_exportToPPTX.exportToPPTX('addslide'); %Uncorrected Vent Image
GasExchange_exportToPPTX.exportToPPTX('addpicture',VentMontage);
GasExchange_exportToPPTX.exportToPPTX('addtext',sprintf(['Uncorrected Ventilation Image-\nSNR: ', num2str(VentSNR)]),'Position',[0 0 5 3]);

GasExchange_exportToPPTX.exportToPPTX('addslide'); %Vent Image
GasExchange_exportToPPTX.exportToPPTX('addpicture',CorrVentMontage);
GasExchange_exportToPPTX.exportToPPTX('addtext',sprintf(['Corrected Ventilation Image-\nUncorrected SNR: ', num2str(VentSNR)]),'Position',[0 0 5 3]);


GasExchange_exportToPPTX.exportToPPTX('addslide'); %Gas Image
GasExchange_exportToPPTX.exportToPPTX('addpicture',GasMontage);
GasExchange_exportToPPTX.exportToPPTX('addtext',sprintf(['Gas Image-\nSNR: ', num2str(GasSNR)]),'Position',[0 0 5 3]);

GasExchange_exportToPPTX.exportToPPTX('addslide'); %Dissolved Image
GasExchange_exportToPPTX.exportToPPTX('addpicture',DissolvedMontage);
GasExchange_exportToPPTX.exportToPPTX('addtext',sprintf(['Dissolved Image-\nSNR: ', num2str(DissolvedSNR)]),'Position',[0 0 5 3]);

GasExchange_exportToPPTX.exportToPPTX('addslide'); %Corrected Dissolved Image
GasExchange_exportToPPTX.exportToPPTX('addpicture',CorrDissolvedMontage);
GasExchange_exportToPPTX.exportToPPTX('addtext',sprintf(['Corrected Dissolved Image-\nSNR: ', num2str(CorrDissolvedSNR)]),'Position',[0 0 5 3]);

GasExchange_exportToPPTX.exportToPPTX('addslide'); %Barrier Image
GasExchange_exportToPPTX.exportToPPTX('addpicture',BarrierMontage);
GasExchange_exportToPPTX.exportToPPTX('addtext',sprintf(['Membrane Image-\nSNR: ', num2str(BarrierSNR)]),'Position',[0 0 5 3]);

GasExchange_exportToPPTX.exportToPPTX('addslide'); %RBC Image
GasExchange_exportToPPTX.exportToPPTX('addpicture',RBCMontage);
GasExchange_exportToPPTX.exportToPPTX('addtext',sprintf(['RBC Image-\nSNR: ', num2str(RBCSNR)]),'Position',[0 0 5 3]);

%Save presentation and close presentation -- overwrite file if it already exists
GasExchange_exportToPPTX.exportToPPTX('save',fullfile(DataLocation, ReportTitle));
if ~exist(fullfile([ParentFolder,'\All Summary Powerpoints']), 'dir')%create powerpowint summary folder if not present
   mkdir(fullfile([ParentFolder,'\All Summary Powerpoints']))
end
GasExchange_exportToPPTX.exportToPPTX('save',fullfile([ParentFolder,'\All Summary Powerpoints'], ReportTitle));
GasExchange_exportToPPTX.exportToPPTX('close');
disp('Exporting Results to Powerpoint Summary Completed.')


%% Save To All Subject Summary
disp('Saving Values of Interest to All Subject Summary Data...')
SubjectMatch = [];
    try
        load([ParentFolder,'\AllSubjectSummary.mat'],'AllSubjectSummary');
        SubjectMatch = find(strcmpi(AllSubjectSummary.Subject,Subject) &... 
                    strcmpi(AllSubjectSummary.Scan_Date,ScanDate) &...
                    strcmpi(AllSubjectSummary.Notes,Notes));
    catch
        disp('Did not find existing AllSubjectSummary table; creating new one...')
        headers = {'Subject', 'Scan_Date', 'Notes', 'Scanner'...%Subject Info
                    'TE90', 'Actual_TE90', 'Flip_Angle_Factor',...%Acquisition Info
                    'Package_Version', 'Process_Date', 'Healthy_Cohort',...%Reconstruction Info
                    'Gas_Frequency', 'Barrier_Frequency', 'RBC_Frequency',...%Spectro Fitting Results - Freq
                    'Gas_T2star', 'Barrier_T2star', 'Barrier_LtoG_ratio', 'RBC_T2star',...%Spectro Fitting Results - T2*
                    'Gas_Phase', 'Barrier_Phase', 'RBC_Phase',...%Spectro Fitting Results - Phase
                    'Dixon_Phase_Separation', 'RBC_Barrier_Ratio', 'Gas_Contamination', 'Dixon_Phase_Shift',...%Reconstruction Parameters
                    'Ventilation_SNR', 'Gas_SNR', 'Dissolved_SNR','CorrDissolved_SNR', 'Barrier_SNR', 'RBC_SNR',...%SNRs
                    'Ventilation_Mean', 'Ventilation_Std_Dev',...%Quantitative Distributions - Vent
                    'Dissolved_Mean', 'Dissolved_Std_Dev',...%Quantitative Distributions - Dissolved
                    'Barrier_Uptake_Mean', 'Barrier_Uptake_Std_Dev',...%Quantitative Distributions - Barrier
                    'RBC_Transfer_Mean', 'RBC_Transfer_Std_Dev',...%Quantitative Distributions - RBC
                    'RBC2Barrier_Mean', 'RBC2Barrier_Std_Dev',...%Quantitative Distributions - RBC:Barrier
                    'RBC_Osc_K0%', 'RBC_Osc_Mean', 'RBC_Osc_Std_Dev',...%Quantitative Distributions - RBC Oscillations
                    'Ventilation_Bin1_Percent', 'Ventilation_Bin2_Percent', 'Ventilation_Bin3_Percent', 'Ventilation_Bin4_Percent', 'Ventilation_Bin5_Percent', 'Ventilation_Bin6_Percent',...%Binning Results - Vent
                    'Dissolved_Bin1_Percent', 'Dissolved_Bin2_Percent', 'Dissolved_Bin3_Percent', 'Dissolved_Bin4_Percent', 'Dissolved_Bin5_Percent', 'Dissolved_Bin6_Percent',...%Binning Results - Dissolved
                    'Barrier_Uptake_Bin1_Percent', 'Barrier_Uptake_Bin2_Percent', 'Barrier_Uptake_Bin3_Percent', 'Barrier_Uptake_Bin4_Percent', 'Barrier_Uptake_Bin5_Percent', 'Barrier_Uptake_Bin6_Percent', 'Barrier_Uptake_Bin7_Percent', 'Barrier_Uptake_Bin8_Percent',...%Binning Results - Barrier
                    'RBC_Transfer_Bin1_Percent', 'RBC_Transfer_Bin2_Percent', 'RBC_Transfer_Bin3_Percent', 'RBC_Transfer_Bin4_Percent', 'RBC_Transfer_Bin5_Percent', 'RBC_Transfer_Bin6_Percent',...%Binning Results - RBC
                    'RBC2Barrier_Bin1_Percent', 'RBC2Barrier_Bin2_Percent', 'RBC2Barrier_Bin3_Percent', 'RBC2Barrier_Bin4_Percent', 'RBC2Barrier_Bin5_Percent', 'RBC2Barrier_Bin6_Percent',...%Binning Results - RBC:Barrier
                    'RBC_Osc_Bin1_Percent', 'RBC_Osc_Bin2_Percent', 'RBC_Osc_Bin3_Percent', 'RBC_Osc_Bin4_Percent', 'RBC_Osc_Bin5_Percent', 'RBC_Osc_Bin6_Percent', 'RBC_Osc_Bin7_Percent', 'RBC_Osc_Bin8_Percent',...%Binning Results - RBC Oscillations
                    };
        AllSubjectSummary = cell2table(cell(0,size(headers,2)));
        AllSubjectSummary.Properties.VariableNames = headers;
    end
    NewData = {Subject, ScanDate, Notes, Scanner,...%Subject Info
                TE90, ActTE90, XeSin.flip_angles.vals(1,1)/DisFlipAngle,...%Acquisition Info
                ReconVersion, datestr(date,29), HealthyCohortNum,...%Reconstruction Info
                GasFreq, BarrierFreq, RBCFreq,...%Spectro Fitting Results - Freq
                GasT2Star, BarrierT2Star, BarrierLtoGRatio, RBCT2Star,...%Spectro Fitting Results - T2*
                AppendedDissolvedNMRFit.phase(3), AppendedDissolvedNMRFit.phase(2), AppendedDissolvedNMRFit.phase(1),...%Spectro Fitting Results - Phase
                SpecDissolvedPhaseDifference, SpecRBCBarrierRatio, SpecGasPhaseContaminationPercentage, Delta_angle_deg,...%Reconstruction Parameters
                VentSNR, GasSNR, DissolvedSNR, CorrDissolvedSNR, BarrierSNR, RBCSNR,...%SNRs
                VentMean, VentStd,...%Quantitative Distributions - Vent
                DissolvedMean, DissolvedStd,...%Quantitative Distributions - Dissolved
                BarrierUptakeMean, BarrierUptakeStd,...%Quantitative Distributions - Barrier
                RBCTransferMean, RBCTransferStd,...%Quantitative Distributions - RBC
                RBCBarrierMean, RBCBarrierStd,...%Quantitative Distributions - RBC:Barrier
                RBC2Bar_struct.RBCk0Osc,RBCOscMean, RBCOscStd,...%Quantitative Distributions - RBC Oscillations
                VentBinPercents(1), VentBinPercents(2), VentBinPercents(3), VentBinPercents(4), VentBinPercents(5), VentBinPercents(6),...%Binning Results - Vent
                DissolvedBinPercents(1), DissolvedBinPercents(2), DissolvedBinPercents(3), DissolvedBinPercents(4), DissolvedBinPercents(5), DissolvedBinPercents(6),...%Binning Results - Dissolved
                BarrierUptakeBinPercents(1), BarrierUptakeBinPercents(2), BarrierUptakeBinPercents(3), BarrierUptakeBinPercents(4), BarrierUptakeBinPercents(5), BarrierUptakeBinPercents(6), BarrierUptakeBinPercents(7), BarrierUptakeBinPercents(8),...%Binning Results - Barrier
                RBCTransferBinPercents(1), RBCTransferBinPercents(2), RBCTransferBinPercents(3), RBCTransferBinPercents(4), RBCTransferBinPercents(5), RBCTransferBinPercents(6),...%Binning Results - RBC
                RBCBarrierBinPercents(1), RBCBarrierBinPercents(2), RBCBarrierBinPercents(3), RBCBarrierBinPercents(4), RBCBarrierBinPercents(5), RBCBarrierBinPercents(6),...%Binning Results - RBC:Barrier
                RBCOscBinPercents(1), RBCOscBinPercents(2), RBCOscBinPercents(3), RBCOscBinPercents(4), RBCOscBinPercents(5), RBCOscBinPercents(6), RBCOscBinPercents(7), RBCOscBinPercents(8), ...%Binning Results - RBC Oscillations
                };
    if (isempty(SubjectMatch))%if no match
        AllSubjectSummary = [AllSubjectSummary;NewData];%append
    else
        AllSubjectSummary(SubjectMatch,:) = NewData;%overwrite
    end
    save([ParentFolder,'\AllSubjectSummary.mat'],'AllSubjectSummary')
    try
        writetable(AllSubjectSummary,[ParentFolder,'\AllSubjectSummary.xlsx'],'Sheet',1)
    catch
        writetable(AllSubjectSummary,[ParentFolder,'\AllSubjectSummary (backup copy).xlsx'],'Sheet',1)
    end


disp('Saving Values of Interest to All Subject Summary Data Completed.')


%% Export Mat File for Cohort/Additional Analysis
disp('Exporting Workspace...')
close all; % no need to save open figures as they will reopen automatically. can replot if needed
save([DataLocation, '\GasExchangeWorkspace.mat']);
disp('Exporting Workspace Completed.')

%% Try making short summary report for QIA meetings
try
    GasExchange_GenerateQIAslides(DataLocation)                                   
catch
    disp("unable to make short summary report for QIA meetings")
end

%% Completed
disp(['Gas Exchange Imaging Process Completed for ', Subject, ': Scanned on ', ScanDate])
close all; % no need to save open figures as they will reopen automatically. can replot if needed
toc(subject_tic);
end
