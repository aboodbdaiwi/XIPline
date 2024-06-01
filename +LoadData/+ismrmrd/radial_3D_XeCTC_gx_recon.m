
function [GasExchange] = radial_3D_XeCTC_gx_recon(MainInput,GasExchange)
%% Working with an existing ISMRMRD data set

% This is a simple example of how to reconstruct images from data
% acquired on a fully sampled cartesian grid
%
% Capabilities:
%   2D/3D
%   multiple slices/slabs
%   multiple contrasts, repetitions
%   
% Limitations:
%   only works with a single encoded space
%   fully sampled k-space (no partial fourier or undersampling)
%   multiple repetitions
%   doesn't handle averages, phases, segments and sets
%   ignores noise scans (no pre-whitening)
% 


%   Co-Author: Abdullah Bdaiwi
%   Work email: abdullah.bdaiwi@cchmc.org
%   Personal email: abdaiwi89@gmail.com
%   Website: https://www.cincinnatichildrens.org/research/divisions/c/cpir

%
clc

cd(MainInput.XeDataLocation)
if exist(MainInput.XeFileName, 'file')
    dset = LoadData.ismrmrd.Dataset(MainInput.XeFileName, 'dataset');
else
    error(['File ' MainInput.XeFileName ' does not exist.  Please generate it.'])
end
GasDataLocation = MainInput.XeDataLocation;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Read some fields from the XML header %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We need to check if optional fields exists before trying to read them

hdr = LoadData.ismrmrd.xml.deserialize(dset.readxml);
%% Read all the data
clc
% Reading can be done one acquisition (or chunk) at a time, 
% but this is much faster for data sets that fit into RAM.
D = dset.readAcquisition();
%% Encoding and reconstruction information
% Matrix size
enc_Nx = hdr.encoding.encodedSpace.matrixSize.x;
enc_Ny = hdr.encoding.encodedSpace.matrixSize.y;
if strcmp(MainInput.Scanner,'Siemens')
    enc_Ny = 5000/2;
end
enc_Nz = hdr.encoding.encodedSpace.matrixSize.z;
rec_Nx = hdr.encoding.reconSpace.matrixSize.x;
rec_Ny = hdr.encoding.reconSpace.matrixSize.y;
rec_Nz = hdr.encoding.reconSpace.matrixSize.z;

% Field of View
enc_FOVx = hdr.encoding.encodedSpace.fieldOfView_mm.x;
enc_FOVy = hdr.encoding.encodedSpace.fieldOfView_mm.y;
enc_FOVz = hdr.encoding.encodedSpace.fieldOfView_mm.z;
rec_FOVx = hdr.encoding.reconSpace.fieldOfView_mm.x;
rec_FOVy = hdr.encoding.reconSpace.fieldOfView_mm.y;
rec_FOVz = hdr.encoding.reconSpace.fieldOfView_mm.z;

%% Define Fixed Variables, Subject to Change with Sequence
%General information
NewImages = 1; %1 forces new images to be made

cd(GasDataLocation)
mkdir([GasDataLocation '\Gas Exchange Analysis']);
outputpath = [GasDataLocation '\Gas Exchange Analysis'];
GasExchange.outputpath = outputpath;
ReconVersion = 'TEST';%If not connected to git, can't determine hash so state test
NumPlotSlices = 7;%must be odd to work as expected
ScanVersion = 'XeCTC';

extraOvs = false;
OvsFactor = 1.0;
% dwell_ms = 0.64/68; %in ms - for 64 points with 0.64ms readout   
dwell_ms = D.head.sample_time_us(1) / 1000;
% if (exist('XeSin','var') && XeSin.non_cart_max_encoding_nrs.vals(1) ~= 63) %R59 oversampling
%     extraOvs = true;
%     OvsFactor = (XeSin.non_cart_max_encoding_nrs.vals(1)+1)/64;
% end
pw = 0.65; %pulse width in ms
DisFlipAngle = hdr.sequenceParameters.flipAngle_deg(2); %Dissolved flip angle
GasFlipAngle = hdr.sequenceParameters.flipAngle_deg(1); %Gas flip angle
try
    freq_jump = hdr.userParameters.userParameterLong(2);
    freq_jump = freq_jump.value;
catch
    try
        freq_jump = hdr.userParameters.userParameterDouble.value/0.02830;
    catch
        freq_jump = 7702; %change in freqeuncy for dissolved and off-res
    end
end
    
%% Import Acquisition Information
disp('Importing Acquisition Information...')
dwell_s = dwell_ms / 1000; %convert to s
Subject = 'Subject: ';   %folderName{1,2:end};
Notes = '';
try
    ScanDate = hdr.measurementInformation.seriesDate;
catch
    ScanDate = hdr.studyInformation.studyDate;
end
XeTR = hdr.sequenceParameters.TR(1)*2/1000;%account for 2 acqs, convert to s
if strcmp(MainInput.Scanner,'Siemens')
    XeTR = hdr.sequenceParameters.TR(1)*2/1000000;
end
TE90 = hdr.sequenceParameters.TE(1);
ActTE90 = TE90; %TE90 from center of pulse
Scanner = [hdr.acquisitionSystemInformation.systemVendor,'-',num2str(hdr.acquisitionSystemInformation.systemFieldStrength_T),'T'];
Xe_AcqMatrix = enc_Nx;
Xe_RecMatrix = 2*Xe_AcqMatrix;
disp('Importing Acquisition Information Completed.')

%% Import Data

cellSize = cell2mat(D.data(1,1));
XeData = zeros(size(cellSize,1),enc_Ny,2);
for i = 1:enc_Ny*2
    if mod(i, 2) == 0
        XeData(:,ceil(i/2),2) = cell2mat(D.data(1,i));
    else
        XeData(:,ceil(i/2),1) = cell2mat(D.data(1,i));
    end 
end

if strcmp(ScanVersion,'XeCTC') || strcmp(ScanVersion,'Duke')
    %Dissolved k-space
    DissolvedKSpace = XeData(:,:,2);
    %Gas k-space
    GasKSpace = XeData(:,:,1);

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
end
disp('Importing Data Completed.')


%% Calculate Trajectories
disp('Calculating Trajectories...')
if strcmp(ScanVersion,'XeCTC') || strcmp(ScanVersion,'Duke')
    del = 2.25;
    if (extraOvs)
        del = del * OvsFactor;
    end

    cellSize = cell2mat(D.traj(1,1));
    XeTraj = zeros(size(cellSize,1),size(cellSize,2),enc_Ny,2);
    for i = 1:enc_Ny*2
        if mod(i, 2) == 0
            XeTraj(:,:,ceil(i/2),2) = cell2mat(D.traj(1,i));
        else
            XeTraj(:,:,ceil(i/2),1) = cell2mat(D.traj(1,i));
        end 
    end
    XeTraj = squeeze(XeTraj(:,:,:,1)); % same Trajectories so use gas only 

    if (extraOvs)
        XeTraj = movmean(XeTraj,OvsFactor);
        XeTraj =  downsample(XeTraj,OvsFactor);
    end
end
disp('Calculating Trajectories Completed.')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Transform Data into Images and Spectra                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
MainInput.CalFileName = dir(fullfile(MainInput.XeDataLocation, '*Calibration.h5'));
MainInput.CalFileName = MainInput.CalFileName.name;
if strcmp(ScanVersion,'XeCTC') || strcmp(ScanVersion,'Duke')
    try
        [GasExResults, ~] = GasExchangeFunctions.XeCTC_Calibration_MRD(MainInput);
    catch
        error('Error. Could not analyze calibration Data file')
        success = 0;
    end
    time = GasExResults.DisFit.t;
    AppendedDissolvedNMRFit = GasExResults.DisFit;
    AppendedDissolvedFit = (time(2)-time(1))*fftshift(fft(AppendedDissolvedNMRFit.calcComponentTimeDomainSignal(time),[],1),1);
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
    if strcmp(ScanVersion,'XeCTC') || strcmp(ScanVersion,'Duke')
       %correction not possible
       CorrectedDissKSpace_SS = DissolvedKSpace_SS;
    else
       CorrectedDissKSpace_SS = GasExchangeFunctions.GasPhaseContaminationRemoval(DissolvedKSpace_SS,GasKSpace_SS,dwell_s,-freq_jump,AppendedDissolvedNMRFit.phase(3),AppendedDissolvedNMRFit.area(3),GasFlipAngle);
    end
    disp('Removing Gas Phase Contamination Completed.')
end
%        CorrectedDissKSpace_SS = GasExchangeFunctions.GasPhaseContaminationRemoval(DissolvedKSpace_SS,GasKSpace_SS,dwell_s,-freq_jump,AppendedDissolvedNMRFit.phase(3),AppendedDissolvedNMRFit.area(3),GasFlipAngle);


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
% sgtitle('Signal Dynamics','FontSize',22)
cd(outputpath)
savefig('SigDynamics.fig')
close(gcf)
% %H - Motion
% subplot(1,2,2);
% hold on
% set(gca,'FontSize',18)
% h1 = plot(HPhaseDynamics,'k','LineWidth',1);
% xl = xlim; yl = ylim;
% patch([motion_ind length(HPulses) length(HPulses) motion_ind], [max(ylim) max(ylim) min(ylim) min(ylim)], [0.5 0.5 0.5],'EdgeColor','none')%Remove due to motion
% uistack(h1, 'top');
% plot(HDynFit,'Color',[0 1 0 0.5],'LineWidth',3);
% plot(HDynFit+3*HDynRMSE,'Color',[0 1 0 0.5])
% plot(HDynFit-3*HDynRMSE,'Color',[0 1 0 0.5])
% xlim(xl);ylim(yl);
% title('H Signal Phase Dynamics')
% box on
% xlabel('Projection')
% ylabel('Signal Phase (deg)')
% legend({'Unused - Motion','Data','Fit'},'Location','best','FontSize',12)
% hold off
% sgtitle('Signal Dynamics','FontSize',22)


%% Reconstruct Xe Images
if(NewImages == 1)
    %Vent Image
    disp('Reconstructing Ventilation Image...')
    UncorrectedVentImage = GasExchangeFunctions.Dissolved_HighResImageRecon(Xe_RecMatrix,GasKSpace_SS,XeTraj_SS/2,PixelShift); %2x Resolution
    if strcmp(MainInput.Scanner,'Siemens')
        for i = 1:size(UncorrectedVentImage,1)
            img = UncorrectedVentImage(:,:,i);
            img = imrotate(img,90);
            UncorrectedVentImage(:,:,i) = flip(img,2);
        end
    end
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
    [VentImage, ~] = GasExchangeFunctions.Vent_BiasCorrection(UncorrectedVentImage,VentMask);
    disp('Ventilation Bias Corrected.')

    %Gas Image
    disp('Reconstructing Gas Image...')
    GasImage = GasExchangeFunctions.Dissolved_LowResImageRecon(Xe_RecMatrix,GasKSpace_SS,XeTraj_SS/2,PixelShift); %2x Resolution
    if strcmp(MainInput.Scanner,'Siemens')
        for i = 1:size(GasImage,1)
            img = GasImage(:,:,i);
            img = imrotate(img,90);
            GasImage(:,:,i) = flip(img,2);
        end
    end
    disp('Reconstructing Gas Image Completed.')

    %Dissolved Image
    disp('Reconstructing Dissolved Image...')
    DissolvedImage = GasExchangeFunctions.Dissolved_LowResImageRecon(Xe_RecMatrix,DissolvedKSpace_SS,XeTraj_SS/2,PixelShift); %2x Resolution
    if strcmp(MainInput.Scanner,'Siemens')
        for i = 1:size(DissolvedImage,1)
            img = DissolvedImage(:,:,i);
            img = imrotate(img,90);
            DissolvedImage(:,:,i) = flip(img,2);
        end
    end    
    disp('Reconstructing Dissolved Image Completed.')
    disp('Reconstructing Corrected Dissolved Image...')
    CorrDissolvedImage = GasExchangeFunctions.Dissolved_LowResImageRecon(Xe_RecMatrix,CorrectedDissKSpace_SS,XeTraj_SS/2,PixelShift); %2x Resolution
    disp('Reconstructing Corrected Dissolved Image Completed.')
elseif (VentMask == ones(size(VentMask)))%Vent mask already exists since not 
    %make success = 0 since indicates no signal
    success = 0;%mark as unsuccessful
end

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
clc
disp('Analyzing RBC Osclications...')
SpecRBCBarrierRatio = AppendedDissolvedNMRFit.area(1)/AppendedDissolvedNMRFit.area(2); %To determine expected ratio in image, only use appended to be closer to be at steady state
[RBCOsc_High_Key,RBCOsc_Low_Key,DissolvedKSpace_SS_Detrend_Normalized,RBC2Bar_struct,OscWorkFlow_Fig] = GasExchangeFunctions.dissolved_RBC_Osc_Keyhole(XeTraj_SS,GasKSpace_SS,CorrectedDissKSpace_SS,XeTR,SpecRBCBarrierRatio);
disp('Analyzing RBC Osclications Completed.')
saveas(OscWorkFlow_Fig,'OscWorkFlow_Fig.png')

%Reconstruct RBC Oscillation Images
disp('Reconstructing RBC keyhole Images...')
RBCOsc_High_Image = GasExchangeFunctions.Dissolved_RBCOscImageRecon(Xe_RecMatrix,RBCOsc_High_Key,XeTraj_SS/2,PixelShift);%2x Resolution
RBCOsc_Low_Image = GasExchangeFunctions.Dissolved_RBCOscImageRecon(Xe_RecMatrix,RBCOsc_Low_Key,XeTraj_SS/2,PixelShift);%2x Resolution
RBCOsc_Normalization = GasExchangeFunctions.Dissolved_RBCOscImageRecon(Xe_RecMatrix,DissolvedKSpace_SS_Detrend_Normalized,XeTraj_SS/2,PixelShift);%2x Resolution
disp('Reconstructing RBC keyhole Images Completed.')

%% store avriables
GasExchange.UncorrectedVentImage = UncorrectedVentImage;
GasExchange.VentImage = VentImage;
GasExchange.GasImage = GasImage;
GasExchange.DissolvedImage = DissolvedImage;
GasExchange.CorrDissolvedImage = CorrDissolvedImage;
GasExchange.AppendedDissolvedNMRFit = AppendedDissolvedNMRFit;
GasExchange.RBC2Bar_struct = RBC2Bar_struct;
GasExchange.RBCOsc_High_Image = RBCOsc_High_Image;
GasExchange.RBCOsc_Low_Image = RBCOsc_Low_Image;
GasExchange.RBCOsc_Normalization = RBCOsc_Normalization;
GasExchange.ActTE90 = ActTE90;
GasExchange.DisFlipAngle = DisFlipAngle;
GasExchange.PixelShift = PixelShift;
GasExchange.freq_jump = freq_jump;
GasExchange.DissolvedNMR = DissolvedNMR;
GasExchange.SigDynamics = SigDynamics;

close all;


end