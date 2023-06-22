
function [UncorrectedVentImage,...
          VentImage,...
          GasImage,...
          DissolvedImage,...
          CorrDissolvedImage,...
          AppendedDissolvedNMRFit,...
          RBC2Bar_struct,...
          RBCOsc_High_Image,...
          RBCOsc_Low_Image,...
          RBCOsc_Normalization,...
          ActTE90,...
          DisFlipAngle,...
          PixelShift,...
          freq_jump,...
          DissolvedNMR,...
          SigDynamics] = LoadData_Gas_GasExchange_Philips_Sin(GasDataLocation,Institute)
%Will process the resulting gas exchange mri data
%   for V3 and XeCTC versions of the data
%   Requires a folder with the Dissolved_Xe raw/lab/sin files,
%   Dissolved_Xe_Mask raw/lab/sin files, and two (and only two!) sets of 
%   data/list files for the dissolved and proton images
%   If no folder is provided, user will be prompted to select the folder
%
%   Syntax:  LoadData_Gas_GasExchange_Philips_Sin(GasDataLocation)
%
%   Inputs:
%      GasDataLocation - string containing the path of the folder
%
%   Outputs:
%       UncorrectedVentImage,...
%       VentImage,...
%       GasImage,...
%       DissolvedImage,...
%       CorrDissolvedImage,...
%       AppendedDissolvedNMRFit,...
%       RBC2Bar_struct,...
%       RBCOsc_High_Image,...
%       RBCOsc_Low_Image,...
%       RBCOsc_Normalization,...
%       ActTE90,...
%       DisFlipAngle,...
%       PixelShift,...
%       freq_jump,...
%       DissolvedNMR,...
%       SigDynamics
%
%   Example: 
%   GasExchange_ProcessSingleSubject('C:\Users\mwillmering\Documents\Subject1')
%
% 
%   Package: https://github.com/cchmc-cpir/CCHMC-Gas-Exchange-Processing-Package
%
%   Author: Matthew Willmering
%   Work email: matthew.willmering@cchmc.org
%   Personal email: matt.willmering@gmail.com
%   Website: https://cpir.cchmc.org/
%
%   Co-Author: Abdullah Bdaiwi
%   Work email: abdullah.bdaiwi@cchmc.org
%   Personal email: abdaiwi89@gmail.com
%   Website: https://cpir.cchmc.org/
%% Define Fixed Variables, Subject to Change with Sequence
%General information
NewImages = 1; %1 forces new images to be made
FunctionDirectory = which('LoadData.LoadData_Gas_GasExchange_Philips_Sin');
idcs = strfind(FunctionDirectory,filesep);%determine location of file separators
FunctionDirectory = FunctionDirectory(1:idcs(end)-1);%remove file

cd(GasDataLocation)
mkdir([GasDataLocation '\Gax Exchange Analysis']);
outputpath = [GasDataLocation '\Gax Exchange Analysis'];
%% Determine files to use
if strcmp(Institute,'XeCTC') == 1    
    ScanVersion = 'Xe_CTC';
%     XeSinFile = dir([GasDataLocation,'\*.sin']);    
elseif strcmp(Institute,'CCHMC') == 1                 
    ScanVersion = 'CCHMC_V3';
%     XeSinFile = dir([GasDataLocation,'\*.sin']);    
end
XeSinFile = dir([GasDataLocation,'\*.sin']);    

% if size(XeSinFile,1) < 1 %didn't find that sequence version
%     XeSinFile = dir([GasDataLocation,'\*CTC_GasExchange_Xe.sin']);
%     ScanVersion = 'Xe_CTC';
% end
DataFiles = dir([GasDataLocation,'\*.data']);
XeDataFile = DataFiles(1);
disp('Looking For Necessary Files Completed.')

%% Import Acquisition Information
if strcmp(ScanVersion,'CCHMC_V3')
    dwell_ms = 1.1049/58; %in ms - for 58 points with 1.1049ms readout
    dwell_s = dwell_ms / 1000; %convert to s
    pw = 0.65; %pulse width in ms
    DisFlipAngle = 20; %Dissolved flip angle
    GasFlipAngle = 0.5; %Gas flip angle
    freq_jump = 7143; %change in freqeuncy for dissolved and off-res
elseif strcmp(ScanVersion,'Xe_CTC')
    dwell_ms = 0.64/68; %in ms - for 64 points with 0.64ms readout
    dwell_s = dwell_ms / 1000; %convert to s
    pw = 0.65; %pulse width in ms
    DisFlipAngle = 20; %Dissolved flip angle
    GasFlipAngle = 0.5; %Gas flip angle
    freq_jump = 7702; %change in freqeuncy for dissolved and off-res
end

disp('Importing Acquisition Information...')
XeSin = GasExchangeFunctions.loadSIN([XeSinFile.folder,'\',XeSinFile.name]);

% ScanDate = [folderName{1,1}(1:4),'-',folderName{1,1}(5:6),'-',folderName{1,1}(7:8)];
if strcmp(ScanVersion,'Xe_CTC')
    XeTR = XeSin.repetition_times.vals(1)*2/1000;%account for 2 acqs, convert to s
else
    XeTR = XeSin.repetition_times.vals(1)*3/1000;%account for 3 acqs, convert to s
end
TE90 = XeSin.echo_times.vals(1);
if strcmp(ScanVersion,'Xe_CTC')
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
if strcmp(ScanVersion,'Xe_CTC') %recon 2x interpolated
    Xe_RecMatrix = 2*Xe_AcqMatrix;
else %CCHMC V3, recon 2x interpolated
    Xe_RecMatrix = 2*Xe_AcqMatrix;
end
disp('Importing Acquisition Information Completed.')

%% Import Data
disp('Importing Data...')
%Load in FIDs
[XeData,XeInfo] = GasExchangeFunctions.loadLISTDATA([XeDataFile.folder,'\',XeDataFile.name]);
XeData = squeeze(XeData);

if strcmp(ScanVersion,'Xe_CTC')
    %Dissolved k-space
    DissolvedKSpace = XeData(:,:,2);
    %Gas k-space
    GasKSpace = XeData(:,:,1);

    %Get sizes of dimensions
    Xe_nprof = size(DissolvedKSpace,2);
    Xe_interleaves = 1;
else
    %Split Xe data into components; formatted as [read, proj, interleave, dyn, mix]
    %Dissolved Spectra
    PreDissolvedFID = XeData(:,1,1,1,2);
    PostDissolvedFID = XeData(:,2,1,1,2);

    %%Gas Spectra - Not Needed
    %%PreGasFID = XeData(:,1,1,2,2);
    %%PostGasFID = XeData(:,2,1,2,2);

    %%Off-Res Spectra - Not needed
    %%PreOffResFID = XeData(:,1,1,3,2);
    %%PostOffResFID = XeData(:,2,1,3,2);

    %Dissolved k-space
    DissolvedKSpaceInit = XeData(:,:,:,1,1);
    %Gas k-space
    GasKSpaceInit = XeData(:,:,:,2,1);
    %%Off-Res k-space - Not needed
    %%OffResKSpace = XeData(:,:,:,3,1);

    %Get sizes of dimensions
    XeSpec_nsamp = size(DissolvedKSpaceInit,1);
    XeImg_nsamp = XeSin.max_encoding_numbers.vals(1)+1;
    Xe_nprof = size(DissolvedKSpaceInit,2);
    Xe_interleaves = size(DissolvedKSpaceInit,3);

    %Trim Kspaces to correct number of readout points
    DissolvedKSpaceInit = DissolvedKSpaceInit(1:XeImg_nsamp,:,:);
    GasKSpaceInit = GasKSpaceInit(1:XeImg_nsamp,:,:);
    
    % Put KSpaces in order and logical format
    XeOrder = GasExchangeFunctions.GetProjectionOrder(XeInfo);
    DissolvedKSpace = GasExchangeFunctions.SortUTEData_AcqOrder(GasExchangeFunctions.reshapeUTEData(DissolvedKSpaceInit),[XeImg_nsamp, Xe_nprof, Xe_interleaves], XeOrder);
    GasKSpace = GasExchangeFunctions.SortUTEData_AcqOrder(GasExchangeFunctions.reshapeUTEData(GasKSpaceInit),[XeImg_nsamp, Xe_nprof, Xe_interleaves], XeOrder);
end
disp('Importing Data Completed.')

%% Calculate Trajectories
disp('Calculating Trajectories...')
if strcmp(ScanVersion,'Xe_CTC')
    XeTraj = GasExchangeFunctions.philipsradialcoords(1.25,2,[XeSinFile.folder,'\',XeSinFile.name]); %1.25us delay, Haltoned Spiral
    XeTraj = permute(XeTraj,[4 3 2 1]);
else
    if (strcmp(Scanner,'3T-R'))%R5.3.1
        XeTraj = GasExchangeFunctions.philipsradialcoords(0.36,1,[FunctionDirectory,'\V3 Scan Info\20191008_162511_Dissolved_Xe_20191008 - 3T-R.sin']); %0.36us delay, GM, from non-spectroscopy version
    else%3T-T1; R5.6.0 dMN
        XeTraj = GasExchangeFunctions.philipsradialcoords(0.36,1,[FunctionDirectory,'\V3 Scan Info\20200210_133229_Dissolved_Xe_20191008 - 3T-T1.sin']); %0.36us delay, GM, from non-spectroscopy version
    end
    XeTraj = permute(XeTraj,[4 3 2 1]);
    XeTraj = GasExchangeFunctions.SortUTETraj_AcqOrder(GasExchangeFunctions.reshapeUTETraj(XeTraj), [XeImg_nsamp, Xe_nprof, Xe_interleaves], XeOrder);%reshape and order
end
disp('Calculating Trajectories Completed.')

%% Determine FOV Shift
%not using planned shifts as techs don't seem to do it consistently
disp('Determining Image offset...')
PixelShiftinit = [0; 0; 0]; %start with no offset
converged = 0; %start not converged
iter = 0;
while (converged == 0) %run until no change in shift
    tempIm = GasExchangeFunctions.Dissolved_LowResImageRecon(Xe_RecMatrix,GasKSpace,XeTraj/2,PixelShiftinit); %make image; 2x Resolution
    [tempIm, ~] = GasExchangeFunctions.Vent_BiasCorrection_Coarse(tempIm);
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


%% Spectra
disp('Fitting Spectrum...')
if strcmp(ScanVersion,'Xe_CTC')
    try
%         calFile = dir(fullfile(GasDataLocation,'*.data'));
        calFile = dir(fullfile(GasDataLocation,'cal','*.data'));
        [GasExResults, ~] = GasExchangeFunctions.XeCTC_Calibration(fullfile(calFile.folder,calFile.name));
    catch
        error('Error. Could not analyze calibration Data file')
    end
    time = GasExResults.DisFit.t;
    AppendedDissolvedNMRFit = GasExResults.DisFit;
    AppendedDissolvedFit = dwell_s*fftshift(fft(AppendedDissolvedNMRFit.calcComponentTimeDomainSignal(time),[],1),1);
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
    PrependedDissolvedNMRFit = GasExchangeFunctions.GasExchange_Spectro.NMR_TimeFit_v(PreDissolvedFID,time,area_guess,freq_guess,fwhm_guesses,fwhmG_guess,phase_guesses,[],[]);
    PrependedDissolvedNMRFit.fitTimeDomainSignal();%initial fit
    PrependedDissolvedNMRFit.setBounds(area_lowerBounds,area_upperBounds,...
        freq_lowerBounds,freq_upperBounds,...
        fwhm_lowerBounds,fwhm_upperBounds,...
        fwhmG_lowerBounds,fwhmG_upperBounds,...
        phase_lowerBounds,phase_upperBounds);%add bounds
    PrependedDissolvedNMRFit.fitTimeDomainSignal();%refit

    %fit appended dissolved
    AppendedDissolvedNMRFit = GasExchangeFunctions.GasExchange_Spectro.NMR_TimeFit_v(PostDissolvedFID,time,area_guess,PrependedDissolvedNMRFit.freq,PrependedDissolvedNMRFit.fwhm,PrependedDissolvedNMRFit.fwhmG,PrependedDissolvedNMRFit.phase,[],[]);
    AppendedDissolvedNMRFit.fitTimeDomainSignal();%initial fit
    AppendedDissolvedNMRFit.setBounds(area_lowerBounds,area_upperBounds,...
        freq_lowerBounds,freq_upperBounds,...
        fwhm_lowerBounds,fwhm_upperBounds,...
        fwhmG_lowerBounds,fwhmG_upperBounds,...
        phase_lowerBounds,phase_upperBounds);%add bounds
    AppendedDissolvedNMRFit.fitTimeDomainSignal();%refit
    AppendedDissolvedFit = dwell_s*fftshift(fft(AppendedDissolvedNMRFit.calcComponentTimeDomainSignal(time),[],1),1);
end
SpecRBCBarrierRatio = AppendedDissolvedNMRFit.area(1)/AppendedDissolvedNMRFit.area(2); %To determine expected ratio in image, only use appended to be closer to be at steady state

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
lgnd = legend({'Spectrum - Magnitude','Spectrum - Real','Spectrum - Imaginary','Fit - Magnitude','Fit - Real','Fit - Imaginary','RBC - Real','Barrier - Real','Gas - Real'},'Location','best','NumColumns',3);
set(lgnd,'color','none');
xlim([-8000 2000])
set(gca, 'XDir','reverse','FontSize',24)
title(['Dissolved Phase Spectra and Fit: RBC/Barrier = ',num2str(AppendedDissolvedNMRFit.area(1)/AppendedDissolvedNMRFit.area(2))])
xlabel('Frequency (Hz)')
ylabel('NMR Signal (a.u.)')
hold off
disp('Fitting Spectrum Completed.')
cd(outputpath)
savefig('DissolvedNMR.fig')
close(gcf)
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
    if strcmp(ScanVersion,'Xe_CTC')
        %correction not possible
        CorrectedDissKSpace_SS = DissolvedKSpace_SS;
    else
        CorrectedDissKSpace_SS = GasExchangeFunctions.GasPhaseContaminationRemoval(DissolvedKSpace_SS,GasKSpace_SS,dwell_s,-freq_jump,AppendedDissolvedNMRFit.phase(3),AppendedDissolvedNMRFit.area(3),GasFlipAngle);
    end
    disp('Removing Gas Phase Contamination Completed.')
end

%% View k0 dynamics
XePulses = (1:size(CorrectedDissKSpace_SS,2))'+SS_ind;%start after steady state
%View Xe Dynamic data
SigDynamics = figure('Name','Signal Dynamics');set(SigDynamics,'WindowState','minimized');
set(SigDynamics,'color','white','Units','inches','Position',[0.25 0.25 16 9])
%Xe - SS and Detrending
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
% sgtitle('Signal Dynamics','FontSize',22)
cd(outputpath)
savefig('SigDynamics.fig')
close(gcf)

%% Reconstruct Xe Images
%Vent Image
disp('Reconstructing Ventilation Image...')
UncorrectedVentImage = GasExchangeFunctions.Dissolved_HighResImageRecon(Xe_RecMatrix,GasKSpace_SS,XeTraj_SS/2,PixelShift); %2x Resolution
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
    VentMask = ones(size(VentMask));
end
[VentImage, ~] = GasExchangeFunctions.Vent_BiasCorrection(UncorrectedVentImage,VentMask);
disp('Ventilation Bias Corrected.')

cd(outputpath)
%Gas Image
disp('Reconstructing Gas Image...')
GasImage = GasExchangeFunctions.Dissolved_LowResImageRecon(Xe_RecMatrix,GasKSpace_SS,XeTraj_SS/2,PixelShift); %2x Resolution
disp('Reconstructing Gas Image Completed.')

%Dissolved Image
disp('Reconstructing Dissolved Image...')
DissolvedImage = GasExchangeFunctions.Dissolved_LowResImageRecon(Xe_RecMatrix,DissolvedKSpace_SS,XeTraj_SS/2,PixelShift); %2x Resolution
disp('Reconstructing Dissolved Image Completed.')
disp('Reconstructing Corrected Dissolved Image...')
CorrDissolvedImage = GasExchangeFunctions.Dissolved_LowResImageRecon(Xe_RecMatrix,CorrectedDissKSpace_SS,XeTraj_SS/2,PixelShift); %2x Resolution
disp('Reconstructing Corrected Dissolved Image Completed.')

%View
VentMontage = figure('Name','Vent Image');set(VentMontage,'WindowState','minimized');
montage(abs(UncorrectedVentImage),'DisplayRange',[0 max(abs(UncorrectedVentImage(:)))])
savefig('VentMontage.fig')
close(gcf)

CorrVentMontage = figure('Name','Corrected Vent Image');set(CorrVentMontage,'WindowState','minimized');
montage(abs(VentImage),'DisplayRange',[0 max(abs(VentImage(:)))])
savefig('CorrVentMontage.fig')
close(gcf)

GasMontage = figure('Name','Gas Image');set(GasMontage,'WindowState','minimized');
montage(abs(GasImage),'DisplayRange',[0 max(abs(GasImage(:)))])
savefig('GasMontage.fig')
close(gcf)

DissolvedMontage = figure('Name','Dissolved Image');set(DissolvedMontage,'WindowState','minimized');
montage(abs(DissolvedImage),'DisplayRange',[0 max(abs(DissolvedImage(:)))])
savefig('DissolvedMontage.fig')
close(gcf)

CorrDissolvedMontage = figure('Name','Corrected Dissolved Image');set(CorrDissolvedMontage,'WindowState','minimized');
montage(abs(CorrDissolvedImage),'DisplayRange',[])
savefig('CorrDissolvedMontage.fig')
close(gcf)


%% Determine Keys for RBC Oscillations
disp('Analyzing RBC Osclications...')
[RBCOsc_High_Key,RBCOsc_Low_Key,DissolvedKSpace_SS_Detrend_Normalized,RBC2Bar_struct,OscWorkFlow_Fig] = GasExchangeFunctions.dissolved_RBC_Osc_Keyhole(XeTraj_SS,GasKSpace_SS,CorrectedDissKSpace_SS,XeTR,SpecRBCBarrierRatio);
disp('Analyzing RBC Osclications Completed.')
saveas(OscWorkFlow_Fig,'OscWorkFlow_Fig.png')

%Reconstruct RBC Oscillation Images
disp('Reconstructing RBC keyhole Images...')
RBCOsc_High_Image = GasExchangeFunctions.Dissolved_RBCOscImageRecon(Xe_RecMatrix,RBCOsc_High_Key,XeTraj_SS/2,PixelShift);%2x Resolution
RBCOsc_Low_Image = GasExchangeFunctions.Dissolved_RBCOscImageRecon(Xe_RecMatrix,RBCOsc_Low_Key,XeTraj_SS/2,PixelShift);%2x Resolution
RBCOsc_Normalization = GasExchangeFunctions.Dissolved_RBCOscImageRecon(Xe_RecMatrix,DissolvedKSpace_SS_Detrend_Normalized,XeTraj_SS/2,PixelShift);%2x Resolution
disp('Reconstructing RBC keyhole Images Completed.')



end
