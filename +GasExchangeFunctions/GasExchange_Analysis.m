
function [GasExchange] = GasExchange_Analysis(GasExchange,Proton,MainInput)
%   Inputs:
%      ProtonImage - 3D array:
%      UncorrectedVentImage:
%      VentImage:
%      H_RecMatrix:
%      ProtonSignal:
%      LungSignal:
%      LungStd:
%      NoiseMask:
%      ProtonMax:
%          
%   Outputs:
%      LungMask
%
%   Example: 
%   LoadData_Proton_GasExchange_Philips_Sin('C:\Users\mwillmering\Documents\Subject1')
%
% 
%   Package: https://github.com/cchmc-cpir/CCHMC-Gas-Exchange-Processing-Package
%
%   Author: Matthew Willmering
%   Work email: matthew.willmering@cchmc.org
%   Personal email: matt.willmering@gmail.com
%   Website: https://cpir.cchmc.org/

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f = waitbar(0,'Reading necessary variables...');
% pause(.5)

% SubjectAge = MainInput.SubjectAge;
%% 
VentImage = GasExchange.VentImage;
GasImage = GasExchange.GasImage;
DissolvedImage = GasExchange.DissolvedImage;
CorrDissolvedImage = GasExchange.CorrDissolvedImage;
try 
    ProtonImage = Proton.ProtonImageHR; % Proton.ProtonImageHR || Proton.Image;
    H_RecMatrix = Proton.H_RecMatrix;
    ProtonImageRegistered = Proton.ProtonHRRegistered; % ProtonHRRegistered || ProtonRegistered
    ProtonMaskRegistered = Proton.ProtonMaskRegistred > 0;
    ProtonMax = Proton.ProtonMax;
    Slices_Co = Proton.Slices_Co;
    Slices_Ax = Proton.Slices_Ax;
catch
    ProtonImage = Proton.Image; % Proton.ProtonImageHR || Proton.Image;
    H_RecMatrix = size(Proton.Image,1);
    ProtonImageRegistered = Proton.ProtonRegistered; % ProtonHRRegistered || ProtonRegistered
    ProtonMaskRegistered = Proton.ProtonMaskRegistred > 0;
    ProtonMax = max(Proton.Image(:));
    Slices_Ax = GasExchangeFunctions.pickSlices(ProtonMaskRegistered, 'ax', 7, 3); % Z-axis (XY planes)
    Slices_Co = GasExchangeFunctions.pickSlices(ProtonMaskRegistered, 'co', 7, 3); % Y-axis (XZ planes)
    Slices_Sa = GasExchangeFunctions.pickSlices(ProtonMaskRegistered, 'sa', 7, 3); % X-axis (YZ planes)

end
Slices_Ax = GasExchangeFunctions.pickSlices(ProtonMaskRegistered, 'ax', 7, 3); % Z-axis (XY planes)
% Slices_Co = GasExchangeFunctions.pickSlices(ProtonMaskRegistered, 'co', 7, 3); % Y-axis (XZ planes)
% Slices_Sa = GasExchangeFunctions.pickSlices(ProtonMaskRegistered, 'sa', 7, 3); % X-axis (YZ planes)

LungMask = GasExchange.LungMask > 0;
RBC2Bar_struct = GasExchange.RBC2Bar_struct;
RBCOsc_High_Image = GasExchange.RBCOsc_High_Image;
RBCOsc_Low_Image = GasExchange.RBCOsc_Low_Image;
RBCOsc_Normalization = GasExchange.RBCOsc_Normalization;

AppendedDissolvedNMRFit = GasExchange.AppendedDissolvedNMRFit;
ActTE90 = GasExchange.ActTE90;
freq_jump = GasExchange.freq_jump;

switch MainInput.HealthyReferenceType
    case 'Default' 
        MainInput.ImportHealthyCohort = 0;
        GasExchange.ImportHealthyCohort = 'yes';
        Scriptlocation = mfilename('fullpath');
        idcs = strfind(Scriptlocation,filesep);%determine location of file separators
        Scriptlocation = [Scriptlocation(1:idcs(end)-1),filesep];%remove file
        HealthyCohortFullPath = fullfile(Scriptlocation, 'HealthyCohort.mat');
    case 'Import' 
        MainInput.ImportHealthyCohort = 1;
        HealthyCohortFullPath = GasExchange.HealthyCohortFullPath;
        GasExchange.ImportHealthyCohort = 'yes';
    case 'Manual' 
        MainInput.ImportHealthyCohort = 2;
        GasExchange.ImportHealthyCohort = 'no';
        VentHealthyMean = GasExchange.VentHealthyMean;
        VentHealthyStd = GasExchange.VentHealthyStd;
        DissolvedHealthyMean = GasExchange.DissolvedHealthyMean;
        DissolvedHealthyStd = GasExchange.DissolvedHealthyStd;
        BarrierHealthyMean = GasExchange.BarrierHealthyMean;
        BarrierHealthyStd = GasExchange.MembraneHealthyMean;
        RBCHealthyMean = GasExchange.RBCHealthyMean;
        RBCtHealthyStd = GasExchange.RBCtHealthyStd;
        RBCBarrHealthyMean = GasExchange.RBCBarrHealthyMean;
        RBCBarrHealthyStd = GasExchange.RBCBarrHealthyStd;
        RBCOscHealthyMean = GasExchange.RBCOscHealthyMean;
        RBCOscHealthyStd = GasExchange.RBCOscHealthyStd;
end
if ~isfield(MainInput, 'CCHMC_DbGxAnalysis') 
    MainInput.CCHMC_DbGxAnalysis = 'no';
end
if strcmp(MainInput.CCHMC_DbGxAnalysis,'yes')
    outputpath = MainInput.gx_analysis_folder;
    GasExchange.outputpath = outputpath;
else
    if ~isfield(MainInput, 'OutputPath') || isempty(MainInput.OutputPath)
        DataLocation = MainInput.XeDataLocation;
        cd(DataLocation)
        mkdir([DataLocation '\GasExchange_Analysis']);
        outputpath = [DataLocation '\GasExchange_Analysis'];
    else
        outputpath = MainInput.OutputPath;
    end

end
cd(outputpath)
GasExchange.outputpath = outputpath;

ProtonMaskMontage = figure('Name','Lung Mask');set(ProtonMaskMontage,'WindowState','minimized');
montage(ProtonMaskRegistered,'DisplayRange',[0 1])%unregistered for these
disp('Lung Mask Completed.')
savefig('ProtonMaskMontage.fig')
close(gcf) 

%% General information
try
    MainInput.AnalysisCode_hash = Global.get_git_hashobject;%use git to find hash of commit used
    MainInput.AnalysisCode_hash = MainInput.AnalysisCode_hash(end-8:40);%remove enter at end
    ReconVersion = MainInput.AnalysisCode_hash;
catch
    ReconVersion = 'TEST';%If not connected to git, can't determine hash so state test
end

NumPlotSlices = 7;%must be odd to work as expected
DisFlipAngle = 20; %Dissolved flip angle
GasFlipAngle = 0.5; %Gas flip angle


%% Import Healthy Distribution, If Available
%Healthy Cohort 
% if(exist([DataLocation,'\HealthyCohort.mat'],'file') == 2) %if Healthy cohort data exists, use
if strcmp(GasExchange.ImportHealthyCohort, 'yes') == 1 && (MainInput.ImportHealthyCohort == 0 || MainInput.ImportHealthyCohort == 1)
%     HealthyFile = dir([DataLocation,'\HealthyCohort.mat']);%get file info
    HealthyData = load(HealthyCohortFullPath);%import all variable but figures
    VentThresh = HealthyData.thresholds.Vent;
    DissolvedThresh = HealthyData.thresholds.Dissolved;
    BarrierThresh = HealthyData.thresholds.Barrier;
    RBCThresh = HealthyData.thresholds.RBC;
    RBCBarrThresh = HealthyData.thresholds.RBCBarr;
    RBCOscThresh = HealthyData.thresholds.RBCOsc;
    HealthyCohortNum = HealthyData.CohortTag; %Healthy Cohort Data
    HealthyDistPresent = 1;
else %user input
    VentThresh = [VentHealthyMean-2*VentHealthyStd,VentHealthyMean-1*VentHealthyStd,VentHealthyMean,VentHealthyMean+1*VentHealthyStd,VentHealthyMean+2*VentHealthyStd];%6
    DissolvedThresh = [DissolvedHealthyMean-2*DissolvedHealthyStd,DissolvedHealthyMean-1*DissolvedHealthyStd,DissolvedHealthyMean,DissolvedHealthyMean+1*DissolvedHealthyStd,DissolvedHealthyMean+2*DissolvedHealthyStd];%6
    BarrierThresh = [BarrierHealthyMean-2*BarrierHealthyStd,BarrierHealthyMean-1*BarrierHealthyStd,BarrierHealthyMean,BarrierHealthyMean+1*BarrierHealthyStd,BarrierHealthyMean+2*BarrierHealthyStd,BarrierHealthyMean+3*BarrierHealthyStd,BarrierHealthyMean+4*BarrierHealthyStd];%8
    RBCThresh = [RBCHealthyMean-2*RBCtHealthyStd,RBCHealthyMean-1*RBCtHealthyStd,RBCHealthyMean,RBCHealthyMean+1*RBCtHealthyStd,RBCHealthyMean+2*RBCtHealthyStd];%6
    RBCBarrThresh = [RBCBarrHealthyMean-2*RBCBarrHealthyStd,RBCBarrHealthyMean-1*RBCBarrHealthyStd,RBCBarrHealthyMean,RBCBarrHealthyMean+1*RBCBarrHealthyStd,RBCBarrHealthyMean+2*RBCBarrHealthyStd];%6
    RBCOscThresh = [RBCOscHealthyMean-2*RBCOscHealthyStd,RBCOscHealthyMean-1*RBCOscHealthyStd,RBCOscHealthyMean,RBCOscHealthyMean+1*RBCOscHealthyStd,RBCOscHealthyMean+2*RBCOscHealthyStd,RBCOscHealthyMean+3*RBCOscHealthyStd,RBCOscHealthyMean+4*RBCOscHealthyStd];%8
    HealthyCohortNum = 'User Defind'; %Healthy Cohort Data
    HealthyDistPresent = 0;
end

%% Obtain Parameters of Interest from Spectrum
% waitbar(.5,f,'Obtain Parameters of Interest from Spectrum...');
% pause(1)
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
%% Create Dixon Images
% waitbar(.10,f,'Create Dixon Images...');
% pause(1)
%Dixon Separation of Dissolved Phases
disp('Dixon separation of Dissolved Phases...')
%Standard Dissolved Image
[BarrierImage,RBCImage,Delta_angle_deg,~] = GasExchangeFunctions.SinglePointDixon(CorrDissolvedImage,-SpecRBCBarrierRatio,GasImage,ProtonMaskRegistered);
[~,RBCOsc_High,~,~] = GasExchangeFunctions.SinglePointDixon(RBCOsc_High_Image,-RBC2Bar_struct.High,GasImage,ProtonMaskRegistered);
[~,RBCOsc_Low,~,~] = GasExchangeFunctions.SinglePointDixon(RBCOsc_Low_Image,-RBC2Bar_struct.Low,GasImage,ProtonMaskRegistered);
[~,RBCOsc_Norm,~,~] = GasExchangeFunctions.SinglePointDixon(RBCOsc_Normalization,-RBC2Bar_struct.Tot,GasImage,ProtonMaskRegistered);
disp('Dixon separation of Dissolved Phases Completed.')
% figure; imslice(abs(DissolvedImage))
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
% waitbar(.10,f,'Quantitative Corrections...');
% pause(1)
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
% waitbar(.15,f,'Calculate Ratios...');
% pause(1)
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
% waitbar(.20,f,'Calculate Binned Maps...');
% pause(1)
disp('Binning Ratios...')
%Create Binning Maps
SixBinMap = [1 0 0; 1 0.7143 0; 0.4 0.7 0.4; 0 1 0; 0 0.57 0.71; 0 0 1]; %Used for Vent and RBC
EightBinMap = [1 0 0; 1 0.7143 0; 0.4 0.7 0.4; 0 1 0; 184/255 226/255 145/255; 243/255 205/255 213/255; 225/255 129/255 162/255; 197/255 27/255 125/255]; %Used for barrier
SixBinRBCBarMap = [197/255 27/255 125/255; 225/255 129/255 162/255; 0.4 0.7 0.4; 0 1 0; 1 0.7143 0; 1 0 0]; %Used for RBC/Barrier ratio

%Vent Binning
VentBinMap = GasExchangeFunctions.BinImages(ScaledVentImage, VentThresh);
VentBinMap = VentBinMap.*ProtonMaskRegistered;%Mask to lung volume
VentBinMask = logical((VentBinMap-1).*ProtonMaskRegistered); %Mask for ventilated areas only

if ProtonMaskRegistered == ones(size(ProtonMaskRegistered))
    %if proton mask is entire FOV, dont mask
    success = 0;%mark as unsuccessful
    VentBinMask = ones(size(VentBinMask));
end

%Dissolved Binning
DissolvedBinMap = GasExchangeFunctions.BinImages(DissGasRatio, DissolvedThresh);
DissolvedBinMap = DissolvedBinMap.*VentBinMask;%Mask to ventilated volume

%Barrier Binning
BarrierBinMap = GasExchangeFunctions.BinImages(BarrGasRatio, BarrierThresh);
BarrierBinMap = BarrierBinMap.*VentBinMask;%Mask to ventilated volume

%RBC Binning
RBCBinMap = GasExchangeFunctions.BinImages(RBCGasRatio, RBCThresh);
RBCBinMap = RBCBinMap.*VentBinMask;%Mask to ventilated volume
RBCBinMask = logical((RBCBinMap-1).*VentBinMask); %Mask to ventilated volume

if VentBinMask == ones(size(VentBinMask))
    %if vent mask is entire FOV, dont mask
    success = 0;%mark as unsuccessful
    RBCBinMask = ones(size(RBCBinMask));
end

%RBC/Barrier Binning
RBCBarrierBinMap = GasExchangeFunctions.BinImages(RBCBarrRatio, RBCBarrThresh);
RBCBarrierBinMap = RBCBarrierBinMap.*VentBinMask;%Mask to ventilated volume

%RBC Oscillation Binning
RBCOscBinMap = GasExchangeFunctions.BinImages(RBCOsc, RBCOscThresh);
RBCOscBinMap = RBCOscBinMap.*RBCBinMask;%Mask to RBC-transfer volume
disp('Binning Ratios Completed.')


%% Calculate SNRs
% waitbar(.25,f,'Calculate SNRs...');
% pause(1)
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


%% Quantify Defects, Intensities, etc.
% waitbar(.30,f,'Quantify Defects, Intensities, etc....');
% pause(1)
disp('Quantifying Images...')
%Ventilation Quantification
VentMean = mean(ScaledVentImage(ProtonMaskRegistered(:)));
VentMedian = median(ScaledVentImage(ProtonMaskRegistered(:)));
VentStd = std(ScaledVentImage(ProtonMaskRegistered(:)));
VentBinPercents = zeros(6,1);
for bin = 1:6
    VentBinPercents(bin) = sum(VentBinMap(:)==bin)/sum(ProtonMaskRegistered(:)==1)*100;
end

%Dissolved/Gas Quantification
DissolvedMean = mean(DissGasRatio(VentBinMask(:)));
DissolvedMedian = median(DissGasRatio(VentBinMask(:)));
DissolvedStd = std(DissGasRatio(VentBinMask(:)));
DissolvedBinPercents = zeros(6,1);
for bin = 1:6
    DissolvedBinPercents(bin) = sum(DissolvedBinMap(:)==bin)/sum(VentBinMask(:)==1)*100;
end

%Barrier/Gas Quantification
BarrierUptakeMean = mean(BarrGasRatio(VentBinMask(:)));
BarrierUptakeMedian = median(BarrGasRatio(VentBinMask(:)));
BarrierUptakeStd = std(BarrGasRatio(VentBinMask(:)));
BarrierUptakeBinPercents = zeros(8,1);
for bin = 1:8
    BarrierUptakeBinPercents(bin) = sum(BarrierBinMap(:)==bin)/sum(VentBinMask(:)==1)*100;
end

%RBC/Gas Quantification
RBCTransferMean = mean(RBCGasRatio(VentBinMask(:)));
RBCTransferMedian = median(RBCGasRatio(VentBinMask(:)));
RBCTransferStd = std(RBCGasRatio(VentBinMask(:)));
RBCTransferBinPercents = zeros(6,1);
for bin = 1:6
    RBCTransferBinPercents(bin) = sum(RBCBinMap(:)==bin)/sum(VentBinMask(:)==1)*100;
end

%RBC/Barrier Quantification
validMask = VentBinMask(:) & isfinite(RBCBarrRatio(:));
validValues = RBCBarrRatio(validMask);
RBCBarrierMean   = nanmean(validValues);
RBCBarrierMedian = nanmedian(validValues);
RBCBarrierStd    = nanstd(validValues);
RBCBarrierBinPercents = zeros(6,1);
for bin = 1:6
    RBCBarrierBinPercents(bin) = sum(RBCBarrierBinMap(:)==bin)/sum(VentBinMask(:)==1)*100;
end

%RBC Oscillation Quantification
%RBCOscBin1Percent
RBCOscMean = mean(RBCOsc(RBCBinMask(:)));
RBCOscMedian = median(RBCOsc(RBCBinMask(:)));
RBCOscStd = std(RBCOsc(RBCBinMask(:)));
RBCOscBinPercents = zeros(8,1);
for bin = 1:8
    RBCOscBinPercents(bin) = sum(RBCOscBinMap(:)==bin)/sum(RBCBinMask(:)==1)*100;
end
disp('Quantifying Images Completed.')
%% Calculate Histograms
% waitbar(.40,f,'Calculate Histograms....');
% pause(1)
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
    VentHistFig = GasExchangeFunctions.CalculateDissolvedHistogram(ScaledVentImage(ProtonMaskRegistered(:)),VentEdges,VentThresh,SixBinMap,HealthyData.VentEdges,HealthyData.HealthyVentFit);
    DissHistFig = GasExchangeFunctions.CalculateDissolvedHistogram(DissGasRatio(VentBinMask(:)),DissEdges,DissolvedThresh,SixBinMap,HealthyData.DissolvedEdges,HealthyData.HealthyDissFit);
    BarHistFig = GasExchangeFunctions.CalculateDissolvedHistogram(BarrGasRatio(VentBinMask(:)),BarEdges,BarrierThresh,EightBinMap,HealthyData.BarsEdges,HealthyData.HealthyBarsFit);
    RBCHistFig = GasExchangeFunctions.CalculateDissolvedHistogram(RBCGasRatio(VentBinMask(:)),RBCEdges,RBCThresh,SixBinMap,HealthyData.RBCEdges,HealthyData.HealthyRBCsFit);
    RBCBarHistFig = GasExchangeFunctions.CalculateDissolvedHistogram(RBCBarrRatio(VentBinMask(:)),RBCBarEdges,RBCBarrThresh,SixBinRBCBarMap,HealthyData.RBCBarEdges,HealthyData.HealthyRBCBarFit);
    RBCOscHistFig = GasExchangeFunctions.CalculateDissolvedHistogram(RBCOsc(RBCBinMask(:)),RBCOscEdges,RBCOscThresh,EightBinMap,HealthyData.RBCOscEdges,HealthyData.HealthyRBCOscFit);
else
    VentHistFig = GasExchangeFunctions.CalculateDissolvedHistogram(ScaledVentImage(ProtonMaskRegistered(:)),VentEdges,VentThresh,SixBinMap,[],[]);
    DissHistFig = GasExchangeFunctions.CalculateDissolvedHistogram(DissGasRatio(VentBinMask(:)),DissEdges,DissolvedThresh,SixBinMap,[],[]);
    BarHistFig = GasExchangeFunctions.CalculateDissolvedHistogram(BarrGasRatio(VentBinMask(:)),BarEdges,BarrierThresh,EightBinMap,[],[]);
    RBCHistFig = GasExchangeFunctions.CalculateDissolvedHistogram(RBCGasRatio(VentBinMask(:)),RBCEdges,RBCThresh,SixBinMap,[],[]);
    RBCBarHistFig = GasExchangeFunctions.CalculateDissolvedHistogram(RBCBarrRatio(VentBinMask(:)),RBCBarEdges,RBCBarrThresh,SixBinRBCBarMap,[],[]);
    RBCOscHistFig = GasExchangeFunctions.CalculateDissolvedHistogram(RBCOsc(RBCBinMask(:)),RBCOscEdges,RBCOscThresh,EightBinMap,[],[]);
end
Ventilationtitle = ['99 Percentile Scaled Ventilation Histogram (Bins %: ',[num2str(VentBinPercents(1),'%1.2f')],', ',[num2str(VentBinPercents(2),'%1.2f')],', ',[num2str(VentBinPercents(3),'%1.2f')],', ',[num2str(VentBinPercents(4),'%1.2f')],', ',[num2str(VentBinPercents(5),'%1.2f')],', ',[num2str(VentBinPercents(6),'%1.2f')],')'];
Dissolvedtitle = ['Dissolved:Gas Histogram (Bins %: ',[num2str(DissolvedBinPercents(1),'%1.2f')],', ',[num2str(DissolvedBinPercents(2),'%1.2f')],', ',[num2str(DissolvedBinPercents(3),'%1.2f')],', ',[num2str(DissolvedBinPercents(4),'%1.2f')],', ',[num2str(DissolvedBinPercents(5),'%1.2f')],', ',[num2str(DissolvedBinPercents(6),'%1.2f')],')'];    
BarrierUptaketitle = ['Barrier:Gas Histogram (Bins %: ',[num2str(BarrierUptakeBinPercents(1),'%1.2f')],', ',[num2str(BarrierUptakeBinPercents(2),'%1.2f')],', ',[num2str(BarrierUptakeBinPercents(3),'%1.2f')],', ',[num2str(BarrierUptakeBinPercents(4),'%1.2f')],', ',[num2str(BarrierUptakeBinPercents(5),'%1.2f')],', ',[num2str(BarrierUptakeBinPercents(6),'%1.2f')],', ',[num2str(BarrierUptakeBinPercents(7),'%1.2f')],', ',[num2str(BarrierUptakeBinPercents(8),'%1.2f')],')'];
RBCTransfertitle = ['RBC:Gas Histogram (Bins %: ',[num2str(RBCTransferBinPercents(1),'%1.2f')],', ',[num2str(RBCTransferBinPercents(2),'%1.2f')],', ',[num2str(RBCTransferBinPercents(3),'%1.2f')],', ',[num2str(RBCTransferBinPercents(4),'%1.2f')],', ',[num2str(RBCTransferBinPercents(5),'%1.2f')],', ',[num2str(RBCTransferBinPercents(6),'%1.2f')],')'];
RBCBarriertitle = ['RBC:Barrier Histogram (Bins %: ',[num2str(RBCBarrierBinPercents(1),'%1.2f')],', ',[num2str(RBCBarrierBinPercents(2),'%1.2f')],', ',[num2str(RBCBarrierBinPercents(3),'%1.2f')],', ',[num2str(RBCBarrierBinPercents(4),'%1.2f')],', ',[num2str(RBCBarrierBinPercents(5),'%1.2f')],', ',[num2str(RBCBarrierBinPercents(6),'%1.2f')],')'];
RBCOsctitle = ['RBC Oscillation Histogram (Bins %: ',[num2str(RBCOscBinPercents(1),'%1.2f')],', ',[num2str(RBCOscBinPercents(2),'%1.2f')],', ',[num2str(RBCOscBinPercents(3),'%1.2f')],', ',[num2str(RBCOscBinPercents(4),'%1.2f')],', ',[num2str(RBCOscBinPercents(5),'%1.2f')],[num2str(RBCOscBinPercents(6),'%1.2f')],', ',[num2str(RBCOscBinPercents(7),'%1.2f')],', ',[num2str(RBCOscBinPercents(8),'%1.2f')],')'];

%Edit Names/Titles
set(VentHistFig,'Name','Ventilation Histogram');
    title(VentHistFig.CurrentAxes,Ventilationtitle);
    ylabel(VentHistFig.CurrentAxes,'Percentage of Voxels');
    xlabel(VentHistFig.CurrentAxes,'Scaled Ventilation (a.u.)');
saveas(VentHistFig,'Ventilation_Histogram.png')
set(DissHistFig,'Name','Dissolved Histogram');
    title(DissHistFig.CurrentAxes,Dissolvedtitle);
    ylabel(DissHistFig.CurrentAxes,'Percentage of Voxels');
    xlabel(DissHistFig.CurrentAxes,'Dissolved/Gas Ratio');
saveas(DissHistFig,'Dissolved_Histogram.png')
set(BarHistFig,'Name','Barrier Histogram');
    title(BarHistFig.CurrentAxes,BarrierUptaketitle);
    ylabel(BarHistFig.CurrentAxes,'Percentage of Voxels');
    xlabel(BarHistFig.CurrentAxes,'Barrier/Gas Ratio (Barrier-uptake)');
saveas(BarHistFig,'Barrier_Histogram.png')
set(RBCHistFig,'Name','RBC Histogram');
    title(RBCHistFig.CurrentAxes,RBCTransfertitle);
    ylabel(RBCHistFig.CurrentAxes,'Percentage of Voxels');
    xlabel(RBCHistFig.CurrentAxes,'RBC/Gas Ratio (RBC-transfer)');
saveas(RBCHistFig,'RBC_Histogram.png')
set(RBCBarHistFig,'Name','RBC:Barrier Histogram');
    title(RBCBarHistFig.CurrentAxes,RBCBarriertitle);
    ylabel(RBCBarHistFig.CurrentAxes,'Percentage of Voxels');
    xlabel(RBCBarHistFig.CurrentAxes,'RBC/Barrier Ratio');
saveas(RBCBarHistFig,'RBC_Barrier_Histogram.png')
set(RBCOscHistFig,'Name','RBC Oscillation Histogram');
    title(RBCOscHistFig.CurrentAxes,RBCOsctitle);
    ylabel(RBCOscHistFig.CurrentAxes,'Percentage of Voxels');
    xlabel(RBCOscHistFig.CurrentAxes,'RBC oscilations (% of RBC Signal)');
saveas(RBCOscHistFig,'RBC_Oscillation_Histogram.png')

disp('Calculating Histograms Completed.')

% save histograms
% print(VentHistFig,[outputpath , '\','Ventilation_Histogram'],'-dpng','-r300');
% print(DissHistFig,[outputpath , '\','Dissolved_Histogram'],'-dpng','-r300');
% print(BarHistFig,[outputpath , '\','Barrier_Histogram'],'-dpng','-r300');
% print(RBCHistFig,[outputpath , '\','RBC_Histogram'],'-dpng','-r300');
% print(RBCBarHistFig,[outputpath ,'\', 'RBC_Barrier_Histogram'],'-dpng','-r300');
% print(RBCOscHistFig,[outputpath , '\','RBC_Oscillation_Histogram'],'-dpng','-r300');
% 

%% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Export Data and Results                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Summary Figures
% waitbar(.50,f,'Summary Figures....');
% pause(1)
%Vent
SumVentFig = figure('Name','Ventitlation','units','normalized','outerposition',[0 0 1 4/NumPlotSlices]);set(SumVentFig,'WindowState','minimized');
set(SumVentFig,'color','white','Units','inches','Position',[0.25 0.25 2*NumPlotSlices 2*2])
tiledlayout(2,NumPlotSlices,'TileSpacing','none','Padding','compact');
for slice=1:NumPlotSlices %Plot Gas Phase Mask Bins slices
    nexttile
    [~,~] = Global.imoverlay(squeeze(abs(ProtonImageRegistered(:,:,Slices_Co(slice),1))),squeeze(VentBinMap(:,:,Slices_Co(slice),1)),[1,6],[0,0.99*ProtonMax],SixBinMap,1,gca);
    colormap(gca,SixBinMap)
    if slice == 1
        colorbar(gca,'north','Ticks',[])
    end
    if slice == round(NumPlotSlices/2)
        title('Ventilation','FontSize',24)
    end
end
for slice=1:NumPlotSlices %Plot Gas Phase Mask Bins slices
    nexttile
    [~,~] = Global.imoverlay(fliplr(rot90(squeeze(abs(ProtonImageRegistered(Slices_Ax(slice),:,:,1))),-1)),fliplr(rot90(squeeze(VentBinMap(Slices_Ax(slice),:,:,1)),-1)),[1,6],[0,0.99*ProtonMax],SixBinMap,1,gca);
    colormap(gca,SixBinMap)
    if slice == 1
        colorbar(gca,'north','Ticks',[])
    end
end
set(SumVentFig,'WindowState','minimized');

%Dissolved
SumDissFig = figure('Name','Dissolved','units','normalized','outerposition',[0 0 1 4/NumPlotSlices]);set(SumDissFig,'WindowState','minimized');
set(SumDissFig,'color','white','Units','inches','Position',[0.25 0.25 2*NumPlotSlices 2*2])
tiledlayout(2,NumPlotSlices,'TileSpacing','none','Padding','compact');
for slice=1:NumPlotSlices %Plot Dissolved Phase Mask Bins slices
    nexttile
    [~,~] = Global.imoverlay(squeeze(abs(ProtonImageRegistered(:,:,Slices_Co(slice),1))),squeeze(DissolvedBinMap(:,:,Slices_Co(slice),1)),[1,6],[0,0.99*ProtonMax],SixBinMap,1,gca);
    colormap(gca,SixBinMap)
    if slice == 1
        colorbar(gca,'north','Ticks',[])
    end
    if slice == round(NumPlotSlices/2)
        title('Dissolved','FontSize',24)
    end
end
for slice=1:NumPlotSlices %Plot Dissolved Phase Mask Bins slices
    nexttile
    [~,~] = Global.imoverlay(fliplr(rot90(squeeze(abs(ProtonImageRegistered(Slices_Ax(slice),:,:,1))),-1)),fliplr(rot90(squeeze(DissolvedBinMap(Slices_Ax(slice),:,:,1)),-1)),[1,6],[0,0.99*ProtonMax],SixBinMap,1,gca);
    colormap(gca,SixBinMap)
    if slice == 1
        colorbar(gca,'north','Ticks',[])
    end
end
set(SumDissFig,'WindowState','minimized');

%Barrier
SumBarrFig = figure('Name','Membrane','units','normalized','outerposition',[0 0 1 4/NumPlotSlices]);set(SumBarrFig,'WindowState','minimized');
set(SumBarrFig,'color','white','Units','inches','Position',[0.25 0.25 2*NumPlotSlices 2*2])
tiledlayout(2,NumPlotSlices,'TileSpacing','none','Padding','compact');
for slice=1:NumPlotSlices %Plot Barrier Phase Mask Bins slices
    nexttile
    [~,~] = Global.imoverlay(squeeze(abs(ProtonImageRegistered(:,:,Slices_Co(slice),1))),squeeze(BarrierBinMap(:,:,Slices_Co(slice),1)),[1,8],[0,0.99*ProtonMax],EightBinMap,1,gca);
    colormap(gca,EightBinMap)
    if slice == 1
        colorbar(gca,'north','Ticks',[])
    end
    if slice == round(NumPlotSlices/2)
        title('Membrane','FontSize',24)
    end
end
for slice=1:NumPlotSlices %Plot Barrier Phase Mask Bins slices
    nexttile
    [~,~] = Global.imoverlay(fliplr(rot90(squeeze(abs(ProtonImageRegistered(Slices_Ax(slice),:,:,1))),-1)),fliplr(rot90(squeeze(BarrierBinMap(Slices_Ax(slice),:,:,1)),-1)),[1,8],[0,0.99*ProtonMax],EightBinMap,1,gca);
    colormap(gca,EightBinMap)
    if slice == 1
        colorbar(gca,'north','Ticks',[])
    end
end
set(SumBarrFig,'WindowState','minimized');

%RBC
SumRBCFig = figure('Name','RBC','units','normalized','outerposition',[0 0 1 4/NumPlotSlices]);set(SumRBCFig,'WindowState','minimized');
set(SumRBCFig,'color','white','Units','inches','Position',[0.25 0.25 2*NumPlotSlices 2*2])
tiledlayout(2,NumPlotSlices,'TileSpacing','none','Padding','compact');
for slice=1:NumPlotSlices %Plot RBC Phase Mask Bins slices
    nexttile
    [~,~] = Global.imoverlay(squeeze(abs(ProtonImageRegistered(:,:,Slices_Co(slice),1))),squeeze(RBCBinMap(:,:,Slices_Co(slice),1)),[1,6],[0,0.99*ProtonMax],SixBinMap,1,gca);
    colormap(gca,SixBinMap)
    if slice == 1
        colorbar(gca,'north','Ticks',[])
    end
    if slice == round(NumPlotSlices/2)
        title('RBC','FontSize',24)
    end
end
for slice=1:NumPlotSlices %Plot RBC Phase Mask Bins slices
    nexttile
    [~,~] = Global.imoverlay(fliplr(rot90(squeeze(abs(ProtonImageRegistered(Slices_Ax(slice),:,:,1))),-1)),fliplr(rot90(squeeze(RBCBinMap(Slices_Ax(slice),:,:,1)),-1)),[1,6],[0,0.99*ProtonMax],SixBinMap,1,gca);
    colormap(gca,SixBinMap)
    if slice == 1
        colorbar(gca,'north','Ticks',[])
    end
end
set(SumRBCFig,'WindowState','minimized');

%RBC:Barrier
SumRBCBarFig = figure('Name','RBC:Barrier Ratio','units','normalized','outerposition',[0 0 1 4/NumPlotSlices]);set(SumRBCBarFig,'WindowState','minimized');
set(SumRBCBarFig,'color','white','Units','inches','Position',[0.25 0.25 2*NumPlotSlices 2*2])
tiledlayout(2,NumPlotSlices,'TileSpacing','none','Padding','compact');
for slice=1:NumPlotSlices %Plot RBC Barrier Mask Bins slices
    nexttile
    [~,~] = Global.imoverlay(squeeze(abs(ProtonImageRegistered(:,:,Slices_Co(slice),1))),squeeze(RBCBarrierBinMap(:,:,Slices_Co(slice),1)),[1,6],[0,0.99*ProtonMax],SixBinRBCBarMap,1,gca);
    colormap(gca,SixBinRBCBarMap)
    if slice == 1
        colorbar(gca,'north','Ticks',[])
    end
    if slice == round(NumPlotSlices/2)
        title('RBC:Membrane','FontSize',24)
    end
end
for slice=1:NumPlotSlices %Plot RBC Barrier Mask Bins slices
    nexttile
    [~,~] = Global.imoverlay(fliplr(rot90(squeeze(abs(ProtonImageRegistered(Slices_Ax(slice),:,:,1))),-1)),fliplr(rot90(squeeze(RBCBarrierBinMap(Slices_Ax(slice),:,:,1)),-1)),[1,6],[0,0.99*ProtonMax],SixBinRBCBarMap,1,gca);
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
    [~,~] = Global.imoverlay(squeeze(abs(ProtonImageRegistered(:,:,Slices_Co(slice),1))),squeeze(RBCOscBinMap(:,:,Slices_Co(slice),1)),[1,8],[0,0.99*ProtonMax],EightBinMap,1,gca);
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
    [~,~] = Global.imoverlay(fliplr(rot90(squeeze(abs(ProtonImageRegistered(Slices_Ax(slice),:,:,1))),-1)),fliplr(rot90(squeeze(RBCOscBinMap(Slices_Ax(slice),:,:,1)),-1)),[1,8],[0,0.99*ProtonMax],EightBinMap,1,gca);
    colormap(gca,EightBinMap)
    if slice == 1
        colorbar(gca,'north','Ticks',[])
    end
end
set(SumRBCOscFig,'WindowState','minimized');
pause(0.001)


%% Save Images
% waitbar(.70,f,'Saving Images....');
% pause(1)
disp('Saving Niftis...')
%Niftis (Intensity Images)
niftiwrite(abs(fliplr(rot90(ProtonImage,-1))),[outputpath,'\ProtonImage.nii'],'Compressed',true);
info = niftiinfo([outputpath,'\ProtonImage.nii.gz']);
info.Description = strcat('Package Version: ', ReconVersion);
niftiwrite(abs(fliplr(rot90(ProtonImage,-1))),[outputpath,'\ProtonImage.nii'],info,'Compressed',true);

niftiwrite(abs(fliplr(rot90(Proton.ProtonHRRegistered,-1))),[outputpath,'\ProtonImageHRReg.nii'],'Compressed',true);
info = niftiinfo([outputpath,'\ProtonImageHRReg.nii.gz']);
info.Description = strcat('Package Version: ', ReconVersion);
niftiwrite(abs(fliplr(rot90(Proton.ProtonHRRegistered,-1))),[outputpath,'\ProtonImageHRReg.nii'],info,'Compressed',true);

niftiwrite(abs(fliplr(rot90(Proton.ProtonRegistered,-1))),[outputpath,'\ProtonRegistered.nii'],'Compressed',true);
info = niftiinfo([outputpath,'\ProtonRegistered.nii.gz']);
info.Description = strcat('Package Version: ', ReconVersion);
niftiwrite(abs(fliplr(rot90(Proton.ProtonRegistered,-1))),[outputpath,'\ProtonRegistered.nii'],info,'Compressed',true);

niftiwrite(abs(fliplr(rot90(LungMask,-1))),[outputpath,'\Mask.nii'],'Compressed',true);
info = niftiinfo([outputpath,'\Mask.nii.gz']);
info.Description = strcat('Package Version: ', ReconVersion);
niftiwrite(abs(fliplr(rot90(LungMask,-1))),[outputpath,'\Mask.nii'],info,'Compressed',true);

niftiwrite(abs(fliplr(rot90(ProtonMaskRegistered,-1))),[outputpath,'\ProtonMaskRegistered.nii'],'Compressed',true);%for obtaining other variables later
info = niftiinfo([outputpath,'\ProtonMaskRegistered.nii.gz']);
info.Description = strcat('Package Version: ', ReconVersion);
niftiwrite(abs(fliplr(rot90(ProtonMaskRegistered,-1))),[outputpath,'\ProtonMaskRegistered.nii'],info,'Compressed',true);

niftiwrite(abs(fliplr(rot90(VentImage,-1))),[outputpath,'\VentImage.nii'],'Compressed',true);
info = niftiinfo([outputpath,'\VentImage.nii.gz']);
info.Description = strcat('Package Version: ', ReconVersion);
niftiwrite(abs(fliplr(rot90(VentImage,-1))),[outputpath,'\VentImage.nii'],info,'Compressed',true);

niftiwrite(abs(fliplr(rot90(GasImage,-1))),[outputpath,'\GasImage.nii'],'Compressed',true);
info = niftiinfo([outputpath,'\GasImage.nii.gz']);
info.Description = strcat('Package Version: ', ReconVersion);
niftiwrite(abs(fliplr(rot90(GasImage,-1))),[outputpath,'\GasImage.nii'],info,'Compressed',true);

niftiwrite(abs(fliplr(rot90(CorrDissolvedImage,-1))),[outputpath,'\DissolvedImage.nii'],'Compressed',true);
info = niftiinfo([outputpath,'\DissolvedImage.nii.gz']);
info.Description = strcat('Package Version: ', ReconVersion);
niftiwrite(abs(fliplr(rot90(CorrDissolvedImage,-1))),[outputpath,'\DissolvedImage.nii'],info,'Compressed',true);

niftiwrite(abs(fliplr(rot90(BarrierImage,-1))),[outputpath,'\BarrierImage.nii'],'Compressed',true);
info = niftiinfo([outputpath,'\BarrierImage.nii.gz']);
info.Description = strcat('Package Version: ', ReconVersion);
niftiwrite(abs(fliplr(rot90(BarrierImage,-1))),[outputpath,'\BarrierImage.nii'],info,'Compressed',true);

niftiwrite(abs(fliplr(rot90(RBCImage,-1))),[outputpath,'\RBCImage.nii'],'Compressed',true);
info = niftiinfo([outputpath,'\RBCImage.nii.gz']);
info.Description = strcat('Package Version: ', ReconVersion);
niftiwrite(abs(fliplr(rot90(RBCImage,-1))),[outputpath,'\RBCImage.nii'],info,'Compressed',true);

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
%     X = print('-RGBImage',['-r',num2str(H_RecMatrix/2)]);%2 inches
     Xdata = getframe(gcf);
     X = Xdata.cdata;    
    if (slice == 1)
        imwrite(X,[outputpath,'\BinnedVent.tif'],'Description',strcat('Package Version: ', ReconVersion,'; Cohort: ', HealthyCohortNum));%write new/ overwrite tiff
    else
        imwrite(X,[outputpath,'\BinnedVent.tif'],'WriteMode','append','Description',strcat('Package Version: ', ReconVersion,'; Cohort: ', HealthyCohortNum));%append tiff
    end
end
% read tiff
cd(outputpath)
tiff_info = imfinfo('BinnedVent.tif'); % return tiff structure, one element per image
BinnedVentmap = uint8(zeros(tiff_info(1).Height ,tiff_info(1).Width ,3,length(tiff_info)));
%concatenate each successive tiff to tiff_stack
for ii = 2 : size(tiff_info, 1)
    temp_tiff = imread('BinnedVent.tif', ii);
    BinnedVentmap(:,:,:,ii) = temp_tiff;
end
BinnedVentmap = permute(BinnedVentmap,[1 2 4 3]);
GasExchange.BinnedVentmap = BinnedVentmap;
disp('Saving Vent Tiff Completed. Saving Dissolved Tiff...')
% S = orthosliceViewer(BinnedVentmap);

%Dissolved Binned
for slice=1:H_RecMatrix %repeat for rest of slices
    Xe = imshow(DissolvedBinMap(:,:,slice),[1,6],'Parent',ax2);%plot Xe image
    set(Xe, 'AlphaData', VentBinMask(:,:,slice));%mask Xe image
    colormap(ax2,SixBinMap);%change colors
    imshow(abs(ProtonImageRegistered(:,:,slice)),[0,0.99*ProtonMax],'Parent',ax1);%plot H
%     X = print('-RGBImage',['-r',num2str(H_RecMatrix/2)]);%2 inches
     Xdata = getframe(gcf);
     X = Xdata.cdata;    
    if (slice == 1)
        imwrite(X,[outputpath,'\BinnedDissolved.tif'],'Description',strcat('Package Version: ', ReconVersion,'; Cohort: ', HealthyCohortNum));%write new/ overwrite tiff
    else
        imwrite(X,[outputpath,'\BinnedDissolved.tif'],'WriteMode','append','Description',strcat('Package Version: ', ReconVersion,'; Cohort: ', HealthyCohortNum));%append tiff
    end
end
% read tiff
cd(outputpath)
tiff_info = imfinfo('BinnedDissolved.tif'); % return tiff structure, one element per image
BinnedDissolvedmap = uint8(zeros(tiff_info(1).Height ,tiff_info(1).Width ,3,length(tiff_info)));
%concatenate each successive tiff to tiff_stack
for ii = 2 : size(tiff_info, 1)
    temp_tiff = imread('BinnedDissolved.tif', ii);
    BinnedDissolvedmap(:,:,:,ii) = temp_tiff;
end
BinnedDissolvedmap = permute(BinnedDissolvedmap,[1 2 4 3]);
GasExchange.BinnedDissolvedmap = BinnedDissolvedmap;
disp('Saving Dissolved Tiff Completed. Saving Barrier Tiff...')
% S = orthosliceViewer(BinnedDissolvedmap);

%Barrier Binned
for slice=1:H_RecMatrix %repeat for rest of slices
    Xe = imshow(BarrierBinMap(:,:,slice),[1,8],'Parent',ax2);%plot Xe image
    set(Xe, 'AlphaData', VentBinMask(:,:,slice));%mask Xe image
    colormap(ax2,EightBinMap);%change colors
    imshow(abs(ProtonImageRegistered(:,:,slice)),[0,0.99*ProtonMax],'Parent',ax1);%plot H
%     X = print('-RGBImage',['-r',num2str(H_RecMatrix/2)]);%2 inches
     Xdata = getframe(gcf);
     X = Xdata.cdata;     
    if (slice == 1)
        imwrite(X,[outputpath,'\BinnedBarrierUptake.tif'],'Description',strcat('Package Version: ', ReconVersion,'; Cohort: ', HealthyCohortNum));%write new/ overwrite tiff
    else
        imwrite(X,[outputpath,'\BinnedBarrierUptake.tif'],'WriteMode','append','Description',strcat('Package Version: ', ReconVersion,'; Cohort: ', HealthyCohortNum));%append tiff
    end
end
% read tiff
cd(outputpath)
tiff_info = imfinfo('BinnedBarrierUptake.tif'); % return tiff structure, one element per image
BinnedBarrierUptakemap = uint8(zeros(tiff_info(1).Height ,tiff_info(1).Width ,3,length(tiff_info)));
%concatenate each successive tiff to tiff_stack
for ii = 2 : size(tiff_info, 1)
    temp_tiff = imread('BinnedBarrierUptake.tif', ii);
    BinnedBarrierUptakemap(:,:,:,ii) = temp_tiff;
end
BinnedBarrierUptakemap = permute(BinnedBarrierUptakemap,[1 2 4 3]);
GasExchange.BinnedBarrierUptakemap = BinnedBarrierUptakemap;
disp('Saving Barrier Tiff Completed. Saving RBC Tiff...')

%RBC Binned
for slice=1:H_RecMatrix %repeat for rest of slices
    Xe = imshow(RBCBinMap(:,:,slice),[1,6],'Parent',ax2);%plot Xe image
    set(Xe, 'AlphaData', VentBinMask(:,:,slice));%mask Xe image
    colormap(ax2,SixBinMap);%change colors
    imshow(abs(ProtonImageRegistered(:,:,slice)),[0,0.99*ProtonMax],'Parent',ax1);%plot H
%     X = print('-RGBImage',['-r',num2str(H_RecMatrix/2)]);%2 inches
     Xdata = getframe(gcf);
     X = Xdata.cdata;    
    if (slice == 1)
        imwrite(X,[outputpath,'\BinnedRBCTransfer.tif'],'Description',strcat('Package Version: ', ReconVersion,'; Cohort: ', HealthyCohortNum));%write new/ overwrite tiff
    else
        imwrite(X,[outputpath,'\BinnedRBCTransfer.tif'],'WriteMode','append','Description',strcat('Package Version: ', ReconVersion,'; Cohort: ', HealthyCohortNum));%append tiff
    end
end
% read tiff
cd(outputpath)
tiff_info = imfinfo('BinnedRBCTransfer.tif'); % return tiff structure, one element per image
BinnedRBCTransfermap = uint8(zeros(tiff_info(1).Height ,tiff_info(1).Width ,3,length(tiff_info)));
%concatenate each successive tiff to tiff_stack
for ii = 2 : size(tiff_info, 1)
    temp_tiff = imread('BinnedRBCTransfer.tif', ii);
    BinnedRBCTransfermap(:,:,:,ii) = temp_tiff;
end
BinnedRBCTransfermap = permute(BinnedRBCTransfermap,[1 2 4 3]);
GasExchange.BinnedRBCTransfermap = BinnedRBCTransfermap;
disp('Saving RBC Tiff Completed. Saving RBC:Barrier Tiff...')

%RBC:Barrier Binned
for slice=1:H_RecMatrix %repeat for rest of slices
    Xe = imshow(RBCBarrierBinMap(:,:,slice),[1,6],'Parent',ax2);%plot Xe image
    set(Xe, 'AlphaData', VentBinMask(:,:,slice));%mask Xe image
    colormap(ax2,SixBinRBCBarMap);%change colors
    imshow(abs(ProtonImageRegistered(:,:,slice)),[0,0.99*ProtonMax],'Parent',ax1);%plot H
%     X = print('-RGBImage',['-r',num2str(H_RecMatrix/2)]);%2 inches
     Xdata = getframe(gcf);
     X = Xdata.cdata;    
    if (slice == 1)
        imwrite(X,[outputpath,'\BinnedRBCBarrier.tif'],'Description',strcat('Package Version: ', ReconVersion,'; Cohort: ', HealthyCohortNum));%write new/ overwrite tiff
    else
        imwrite(X,[outputpath,'\BinnedRBCBarrier.tif'],'WriteMode','append','Description',strcat('Package Version: ', ReconVersion,'; Cohort: ', HealthyCohortNum));%append tiff
    end
end
% read tiff
cd(outputpath)
tiff_info = imfinfo('BinnedRBCBarrier.tif'); % return tiff structure, one element per image
BinnedRBCBarriermap = uint8(zeros(tiff_info(1).Height ,tiff_info(1).Width ,3,length(tiff_info)));
%concatenate each successive tiff to tiff_stack
for ii = 2 : size(tiff_info, 1)
    temp_tiff = imread('BinnedRBCBarrier.tif', ii);
    BinnedRBCBarriermap(:,:,:,ii) = temp_tiff;
end
BinnedRBCBarriermap = permute(BinnedRBCBarriermap,[1 2 4 3]);
GasExchange.BinnedRBCBarriermap = BinnedRBCBarriermap;
disp('Saving RBC:Barrier Tiff Completed. Saving RBC Oscillation Tiff...')

%RBC Oscillation Binned
for slice=1:H_RecMatrix %repeat for rest of slices
    Xe = imshow(RBCOscBinMap(:,:,slice),[1,8],'Parent',ax2);%plot Xe image
    set(Xe, 'AlphaData', RBCBinMask(:,:,slice));%mask Xe image
    colormap(ax2,EightBinMap);%change colors
    imshow(abs(ProtonImageRegistered(:,:,slice)),[0,0.99*ProtonMax],'Parent',ax1);%plot H
%     X = print('-RGBImage',['-r',num2str(H_RecMatrix/2)]);%2 inches
     Xdata = getframe(gcf);
     X = Xdata.cdata;     
    if (slice == 1)
        imwrite(X,[outputpath,'\BinnedRBCOscillation.tif'],'Description',strcat('Package Version: ', ReconVersion,'; Cohort: ', HealthyCohortNum));%write new/ overwrite tiff
    else
        imwrite(X,[outputpath,'\BinnedRBCOscillation.tif'],'WriteMode','append','Description',strcat('Package Version: ', ReconVersion,'; Cohort: ', HealthyCohortNum));%append tiff
    end
end
% read tiff
cd(outputpath)
tiff_info = imfinfo('BinnedRBCOscillation.tif'); % return tiff structure, one element per image
BinnedRBCOscillationmap = uint8(zeros(tiff_info(1).Height ,tiff_info(1).Width ,3,length(tiff_info)));
%concatenate each successive tiff to tiff_stack
for ii = 2 : size(tiff_info, 1)
    temp_tiff = imread('BinnedRBCOscillation.tif', ii);
    BinnedRBCOscillationmap(:,:,:,ii) = temp_tiff;
end
BinnedRBCOscillationmap = permute(BinnedRBCOscillationmap,[1 2 4 3]);
GasExchange.BinnedRBCOscillationmap = BinnedRBCOscillationmap;
disp('Saving RBC Oscillation Tiff Completed.')
close(tiff)


%% Save Report
% waitbar(.90,f,'Save Report....');
% pause(1)
disp('Exporting Results to Powerpoint Summary...')
%Start new presentation
isOpen  = Global.exportToPPTX();
if ~isempty(isOpen) %If PowerPoint already started, then close first and then open a new one
    Global.exportToPPTX('close');
end
ReportTitle = 'GasExchange_Analysis';
Global.exportToPPTX('new','Dimensions',[16 9], ...
    'Title',ReportTitle, ...
    'Author','CPIR @ CCHMC');
Subject='Subject';
%Add slides
DissolvedNMR = openfig('DissolvedNMR.fig');
Global.exportToPPTX('addslide'); %Spectrum
Global.exportToPPTX('addpicture',DissolvedNMR,'Scale','maxfixed');
Global.exportToPPTX('addtext',sprintf('129Xe Dissolved-Phase Spectrum'),'Position',[0 0 5 3]);

SigDynamics = openfig('SigDynamics.fig');
Global.exportToPPTX('addslide'); %Signal Dynamics
Global.exportToPPTX('addpicture',SigDynamics);
Global.exportToPPTX('addtext',sprintf('Signal Dynamics'),'Position',[0 0 5 3]);
try
    cd([Proton.folder, '\GasExchange_Analysis'])
    ProtonMontage = openfig('ProtonMontage.fig');
catch
    if MainInput.NoProtonImage == 0
        folderName = 'GasExchange_Analysis';
        newFolder = fullfile(MainInput.HDataLocation, folderName);
        mkdir(newFolder);
        cd(newFolder);
        ProtonMontage = figure('Name','Proton Image');set(ProtonMontage,'WindowState','minimized');
        montage(Proton.ProtonRegistered,'DisplayRange',[0 max(Proton.ProtonRegistered(:))*0.8])%unregistered for these
        savefig('ProtonMontage.fig')
        close(gcf)
    else
        ProtonMontage = figure('Name','Proton Image');set(ProtonMontage,'WindowState','minimized');
        montage(Proton.ProtonRegistered,'DisplayRange',[])%unregistered for these
        savefig('ProtonMontage.fig')
        close(gcf)
    end
end
if MainInput.NoProtonImage == 0
    Global.exportToPPTX('addslide'); %Proton Image
    Global.exportToPPTX('addpicture',ProtonMontage);
    Global.exportToPPTX('addtext',sprintf('Proton Image'),'Position',[0 0 5 3]);
end
cd(outputpath)
ProtonMaskMontage = openfig('ProtonMaskMontage.fig');
Global.exportToPPTX('addslide'); %Proton Mask
Global.exportToPPTX('addpicture',ProtonMaskMontage);
Global.exportToPPTX('addtext',sprintf('Proton Mask'),'Position',[0 0 5 3]);

VentMontage = openfig('VentMontage.fig');
Global.exportToPPTX('addslide'); %Uncorrected Vent Image
Global.exportToPPTX('addpicture',VentMontage);
Global.exportToPPTX('addtext',sprintf(['Uncorrected Ventilation Image-\nSNR: ', num2str(VentSNR)]),'Position',[0 0 5 3]);

CorrVentMontage = openfig('CorrVentMontage.fig');
Global.exportToPPTX('addslide'); %Vent Image
Global.exportToPPTX('addpicture',CorrVentMontage);
Global.exportToPPTX('addtext',sprintf(['Corrected Ventilation Image-\nUncorrected SNR: ', num2str(VentSNR)]),'Position',[0 0 5 3]);

try
    Registrationfig = openfig('Registerationfig.fig');
catch
    Registrationfig = figure('Name','Registration','units','normalized','outerposition',[0 0 1 4/NumPlotSlices]);set(Registrationfig,'WindowState','minimized');
    set(Registrationfig,'color','white','Units','inches','Position',[0.25 0.25 2*NumPlotSlices 4*2+1])
    tiledlayout(4,NumPlotSlices,'TileSpacing','none','Padding','compact');
    for slice=1:NumPlotSlices
        nexttile
        imshowpair(abs(ProtonImageRegistered(:,:,Slices_Co(slice)))/ProtonMax,abs(VentImage(:,:,Slices_Co(slice)))/max(abs(VentImage(ProtonMaskRegistered(:)))),'Scaling','none','ColorChannels','red-cyan')
        if slice == round(NumPlotSlices/2)
            title('Before Registration','FontSize',24)
        end
    end
    for slice=1:NumPlotSlices
        nexttile
        imshowpair(fliplr(rot90(squeeze(abs(ProtonImageRegistered(Slices_Ax(slice),:,:))),-1))/ProtonMax,fliplr(rot90(squeeze(abs(VentImage(Slices_Ax(slice),:,:))),-1))/max(abs(VentImage(ProtonMaskRegistered(:)))),'Scaling','none','ColorChannels','red-cyan')
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
    savefig('Registerationfig.fig')
end
Registration = openfig('Registerationfig.fig');
Global.exportToPPTX('addslide'); %Registration
Global.exportToPPTX('addpicture',Registration);
Global.exportToPPTX('addtext',sprintf('Registration'),'Position',[0 0 5 3]);

GasMontage = openfig('GasMontage.fig');
Global.exportToPPTX('addslide'); %Gas Image
Global.exportToPPTX('addpicture',GasMontage);
Global.exportToPPTX('addtext',sprintf(['Gas Image-\nSNR: ', num2str(GasSNR)]),'Position',[0 0 5 3]);

DissolvedMontage = openfig('DissolvedMontage.fig');
Global.exportToPPTX('addslide'); %Dissolved Image
Global.exportToPPTX('addpicture',DissolvedMontage);
Global.exportToPPTX('addtext',sprintf(['Dissolved Image-\nSNR: ', num2str(DissolvedSNR)]),'Position',[0 0 5 3]);

CorrDissolvedMontage = openfig('CorrDissolvedMontage.fig');
Global.exportToPPTX('addslide'); %Corrected Dissolved Image
Global.exportToPPTX('addpicture',CorrDissolvedMontage);
Global.exportToPPTX('addtext',sprintf(['Corrected Dissolved Image-\nSNR: ', num2str(CorrDissolvedSNR)]),'Position',[0 0 5 3]);

Global.exportToPPTX('addslide'); %Barrier Image
Global.exportToPPTX('addpicture',BarrierMontage);
Global.exportToPPTX('addtext',sprintf(['Barrier Image-\nSNR: ', num2str(BarrierSNR)]),'Position',[0 0 5 3]);

Global.exportToPPTX('addslide'); %RBC Image
Global.exportToPPTX('addpicture',RBCMontage);
Global.exportToPPTX('addtext',sprintf(['RBC Image-\nSNR: ', num2str(RBCSNR)]),'Position',[0 0 5 3]);

Global.exportToPPTX('addslide'); %Vent Summary
Global.exportToPPTX('addpicture',SumVentFig,'Position',[(16-2*NumPlotSlices)/2 0 2*NumPlotSlices 2*2]);
Global.exportToPPTX('addpicture',VentHistFig,'Position',[0 4.5 8 4.5]);
if (HealthyDistPresent)
    [~,pm,~,zm] = ztest(VentMean,HealthyData.HealthyMeans.MeanVent,HealthyData.HealthyMeans.StdVent);
    Global.exportToPPTX('addtext',sprintf(['Mean: ',num2str(VentMean,'%1.3f'),' (Healthy = ',num2str(HealthyData.HealthyMeans.MeanVent,'%1.3f'),'±',num2str(HealthyData.HealthyMeans.StdVent,'%1.3f'),'); Z = ',num2str(zm,'%1.3f'),', P = ',num2str(pm,'%1.3f')]),'Position',[8 4.5 8 4.5],'HorizontalAlignment','center');
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
    Global.exportToPPTX('addtext',sprintf(['Mean: ',num2str(VentMean,'%1.3f'),'±',num2str(VentStd,'%1.3f')]),'Position',[8 4.5 8 4.5],'HorizontalAlignment','center');
    BinTable   = { ...
    '',{'Defect','BackgroundColor',SixBinMap(1,:),'FontWeight','bold'},{'Low','BackgroundColor',SixBinMap(2,:),'FontWeight','bold'},{'Healthy','BackgroundColor',SixBinMap(3,:),'FontWeight','bold'},{'Healthy','BackgroundColor',SixBinMap(4,:),'FontWeight','bold'},{'High','BackgroundColor',SixBinMap(5,:),'FontWeight','bold','Color','w'},{'High','BackgroundColor',SixBinMap(6,:),'FontWeight','bold','Color','w'};...
    'Subject','Bin 1','Bin 2','Bin 3','Bin 4','Bin 5','Bin 6';...
    Subject,num2str(VentBinPercents(1),'%1.2f'),num2str(VentBinPercents(2),'%1.2f'),num2str(VentBinPercents(3),'%1.2f'),num2str(VentBinPercents(4),'%1.2f'),num2str(VentBinPercents(5),'%1.2f'),num2str(VentBinPercents(6),'%1.2f');...
    };
end
Global.exportToPPTX('addtable',BinTable,'Position',[8 5 7.75 2],'Vert','middle','Horiz','center','FontSize',12);
Global.exportToPPTX('addtext',sprintf('Ventilation Summary'),'Position',[0 0 5 3]);


Global.exportToPPTX('addslide'); %Dissolved Summary
Global.exportToPPTX('addpicture',SumDissFig,'Position',[(16-2*NumPlotSlices)/2 0 2*NumPlotSlices 2*2]);
Global.exportToPPTX('addpicture',DissHistFig,'Position',[0 4.5 8 4.5]);
if (HealthyDistPresent)
    [~,pm,~,zm] = ztest(DissolvedMean,HealthyData.HealthyMeans.MeanDissolved,HealthyData.HealthyMeans.StdDissolved);
    Global.exportToPPTX('addtext',sprintf(['Mean: ',num2str(DissolvedMean,'%1.5f'),' (Healthy = ',num2str(HealthyData.HealthyMeans.MeanDissolved,'%1.5f'),'±',num2str(HealthyData.HealthyMeans.StdDissolved,'%1.5f'),'); Z = ',num2str(zm,'%1.3f'),', P = ',num2str(pm,'%1.3f')]),'Position',[8 4.5 8 4.5],'HorizontalAlignment','center');
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
    Global.exportToPPTX('addtext',sprintf(['Mean: ',num2str(DissolvedMean,'%1.5f'),'±',num2str(DissolvedStd,'%1.5f')]),'Position',[8 4.5 8 4.5],'HorizontalAlignment','center');
    BinTable   = { ...
    '',{'Defect','BackgroundColor',SixBinMap(1,:),'FontWeight','bold'},{'Low','BackgroundColor',SixBinMap(2,:),'FontWeight','bold'},{'Healthy','BackgroundColor',SixBinMap(3,:),'FontWeight','bold'},{'Healthy','BackgroundColor',SixBinMap(4,:),'FontWeight','bold'},{'High','BackgroundColor',SixBinMap(5,:),'FontWeight','bold','Color','w'},{'High','BackgroundColor',SixBinMap(6,:),'FontWeight','bold','Color','w'};...
    'Subject','Bin 1','Bin 2','Bin 3','Bin 4','Bin 5','Bin 6';...
    Subject,num2str(DissolvedBinPercents(1),'%1.2f'),num2str(DissolvedBinPercents(2),'%1.2f'),num2str(DissolvedBinPercents(3),'%1.2f'),num2str(DissolvedBinPercents(4),'%1.2f'),num2str(DissolvedBinPercents(5),'%1.2f'),num2str(DissolvedBinPercents(6),'%1.2f');...
    };
end
Global.exportToPPTX('addtable',BinTable,'Position',[8 5 7.75 2],'Vert','middle','Horiz','center','FontSize',12);
Global.exportToPPTX('addtext',sprintf('Dissolved Summary'),'Position',[0 0 5 3]);

Global.exportToPPTX('addslide'); %Barrier Summary
Global.exportToPPTX('addpicture',SumBarrFig,'Position',[(16-2*NumPlotSlices)/2 0 2*NumPlotSlices 2*2]);
Global.exportToPPTX('addpicture',BarHistFig,'Position',[0 4.5 8 4.5]);
if (HealthyDistPresent)
    [~,pm,~,zm] = ztest(BarrierUptakeMean,HealthyData.HealthyMeans.MeanBarrier,HealthyData.HealthyMeans.StdBarrier);
    Global.exportToPPTX('addtext',sprintf(['Mean: ',num2str(BarrierUptakeMean,'%1.5f'),' (Healthy = ',num2str(HealthyData.HealthyMeans.MeanBarrier,'%1.5f'),'±',num2str(HealthyData.HealthyMeans.StdBarrier,'%1.5f'),'); Z = ',num2str(zm,'%1.3f'),', P = ',num2str(pm,'%1.3f')]),'Position',[8 4.5 8 4.5],'HorizontalAlignment','center');
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
    Global.exportToPPTX('addtext',sprintf(['Mean: ',num2str(BarrierUptakeMean,'%1.5f'),'±',num2str(BarrierUptakeStd,'%1.5f')]),'Position',[8 4.5 8 4.5],'HorizontalAlignment','center');
        BinTable   = { ...
    '',{'Defect','BackgroundColor',EightBinMap(1,:),'FontWeight','bold'},{'Low','BackgroundColor',EightBinMap(2,:),'FontWeight','bold'},{'Healthy','BackgroundColor',EightBinMap(3,:),'FontWeight','bold'},{'Healthy','BackgroundColor',EightBinMap(4,:),'FontWeight','bold'},{'Elevated','BackgroundColor',EightBinMap(5,:),'FontWeight','bold'},{'Elevated','BackgroundColor',EightBinMap(6,:),'FontWeight','bold'},{'High','BackgroundColor',EightBinMap(7,:),'FontWeight','bold','Color','w'},{'High','BackgroundColor',EightBinMap(8,:),'FontWeight','bold','Color','w'};...
    'Subject','Bin 1','Bin 2','Bin 3','Bin 4','Bin 5','Bin 6','Bin 7','Bin 8';...
    Subject,num2str(BarrierUptakeBinPercents(1),'%1.2f'),num2str(BarrierUptakeBinPercents(2),'%1.2f'),num2str(BarrierUptakeBinPercents(3),'%1.2f'),num2str(BarrierUptakeBinPercents(4),'%1.2f'),num2str(BarrierUptakeBinPercents(5),'%1.2f'),num2str(BarrierUptakeBinPercents(6),'%1.2f'),num2str(BarrierUptakeBinPercents(7),'%1.2f'),num2str(BarrierUptakeBinPercents(8),'%1.2f');...
    };
end
Global.exportToPPTX('addtable',BinTable,'Position',[8 5 7.75 2],'Vert','middle','Horiz','center','FontSize',12);
Global.exportToPPTX('addtext',sprintf('Barrier-uptake Summary'),'Position',[0 0 5 3]);

Global.exportToPPTX('addslide'); %RBC Summary
Global.exportToPPTX('addpicture',SumRBCFig,'Position',[(16-2*NumPlotSlices)/2 0 2*NumPlotSlices 2*2]);
Global.exportToPPTX('addpicture',RBCHistFig,'Position',[0 4.5 8 4.5]);
if (HealthyDistPresent)
    [~,pm,~,zm] = ztest(RBCTransferMean,HealthyData.HealthyMeans.MeanRBC,HealthyData.HealthyMeans.StdRBC);
    Global.exportToPPTX('addtext',sprintf(['Mean: ',num2str(RBCTransferMean,'%1.3e'),' (Healthy = ',num2str(HealthyData.HealthyMeans.MeanRBC,'%1.3e'),'±',num2str(HealthyData.HealthyMeans.StdRBC,'%1.3e'),'); Z = ',num2str(zm,'%1.3f'),', P = ',num2str(pm,'%1.3f')]),'Position',[8 4.5 8 4.5],'HorizontalAlignment','center');
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
    Global.exportToPPTX('addtext',sprintf(['Mean: ',num2str(RBCTransferMean,'%1.3e'),'±',num2str(RBCTransferStd,'%1.3f')]),'Position',[8 4.5 8 4.5],'HorizontalAlignment','center');
    BinTable   = { ...
    '',{'Defect','BackgroundColor',SixBinMap(1,:),'FontWeight','bold'},{'Low','BackgroundColor',SixBinMap(2,:),'FontWeight','bold'},{'Healthy','BackgroundColor',SixBinMap(3,:),'FontWeight','bold'},{'Healthy','BackgroundColor',SixBinMap(4,:),'FontWeight','bold'},{'High','BackgroundColor',SixBinMap(5,:),'FontWeight','bold','Color','w'},{'High','BackgroundColor',SixBinMap(6,:),'FontWeight','bold','Color','w'};...
    'Subject','Bin 1','Bin 2','Bin 3','Bin 4','Bin 5','Bin 6';...
    Subject,num2str(RBCTransferBinPercents(1),'%1.2f'),num2str(RBCTransferBinPercents(2),'%1.2f'),num2str(RBCTransferBinPercents(3),'%1.2f'),num2str(RBCTransferBinPercents(4),'%1.2f'),num2str(RBCTransferBinPercents(5),'%1.2f'),num2str(RBCTransferBinPercents(6),'%1.2f');...
    };
end
Global.exportToPPTX('addtable',BinTable,'Position',[8 5 7.75 2],'Vert','middle','Horiz','center','FontSize',12);
Global.exportToPPTX('addtext',sprintf('RBC-Transfer Summary'),'Position',[0 0 5 3]);

Global.exportToPPTX('addslide'); %RBC:Barrier Summary
Global.exportToPPTX('addpicture',SumRBCBarFig,'Position',[(16-2*NumPlotSlices)/2 0 2*NumPlotSlices 2*2]);
Global.exportToPPTX('addpicture',RBCBarHistFig,'Position',[0 4.5 8 4.5]);
if (HealthyDistPresent)
    [~,pm,~,zm] = ztest(RBCBarrierMean,HealthyData.HealthyMeans.MeanRBCBar,HealthyData.HealthyMeans.StdRBCBar);
    Global.exportToPPTX('addtext',sprintf(['Mean: ',num2str(RBCBarrierMean,'%1.3f'),' (Healthy = ',num2str(HealthyData.HealthyMeans.MeanRBCBar,'%1.3f'),'±',num2str(HealthyData.HealthyMeans.StdRBCBar,'%1.3f'),'); Z = ',num2str(zm,'%1.3f'),', P = ',num2str(pm,'%1.3f')]),'Position',[8 4.5 8 4.5],'HorizontalAlignment','center');
    p=zeros(6,1);
    z=zeros(6,1);
    for bin = 1:6
        [~,p(bin),~,z(bin)] = ztest(RBCBarrierBinPercents(bin),HealthyData.BinPercentMeans.RBCBar(bin),HealthyData.BinPercentStds.RBCBar(bin));
    end
    BinTable   = { ...
    '',{'High Barrier','BackgroundColor',SixBinRBCBarMap(1,:),'FontWeight','bold','Color','w'},{'Elevated Barrier','BackgroundColor',SixBinRBCBarMap(2,:),'FontWeight','bold','Color','w'},{'Healthy','BackgroundColor',SixBinRBCBarMap(3,:),'FontWeight','bold'},{'Healthy','BackgroundColor',SixBinRBCBarMap(4,:),'FontWeight','bold'},{'Elevated RBC','BackgroundColor',SixBinRBCBarMap(5,:),'FontWeight','bold'},{'High RBC','BackgroundColor',SixBinRBCBarMap(6,:),'FontWeight','bold'};...
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
    Global.exportToPPTX('addtext',sprintf(['Mean: ',num2str(RBCBarrierMean,'%1.3f'),'±',num2str(RBCBarrierStd,'%1.3f')]),'Position',[8 4.5 8 4.5],'HorizontalAlignment','center');
    BinTable   = { ...
    '',{'Defect','BackgroundColor',SixBinRBCBarMap(1,:),'FontWeight','bold','Color','w'},{'Low','BackgroundColor',SixBinRBCBarMap(2,:),'FontWeight','bold','Color','w'},{'Healthy','BackgroundColor',SixBinRBCBarMap(3,:),'FontWeight','bold'},{'Healthy','BackgroundColor',SixBinRBCBarMap(4,:),'FontWeight','bold'},{'High','BackgroundColor',SixBinRBCBarMap(5,:),'FontWeight','bold'},{'High','BackgroundColor',SixBinRBCBarMap(6,:),'FontWeight','bold'};...
    'Subject','Bin 1','Bin 2','Bin 3','Bin 4','Bin 5','Bin 6';...
    Subject,num2str(RBCBarrierBinPercents(1),'%1.2f'),num2str(RBCBarrierBinPercents(2),'%1.2f'),num2str(RBCBarrierBinPercents(3),'%1.2f'),num2str(RBCBarrierBinPercents(4),'%1.2f'),num2str(RBCBarrierBinPercents(5),'%1.2f'),num2str(RBCBarrierBinPercents(6),'%1.2f');...
    };
end
Global.exportToPPTX('addtable',BinTable,'Position',[8 5 7.75 2],'Vert','middle','Horiz','center','FontSize',12);
Global.exportToPPTX('addtext',sprintf('RBC:Barrier Summary'),'Position',[0 0 5 3]);

Global.exportToPPTX('addslide'); %RBC Osc Summary
Global.exportToPPTX('addpicture',SumRBCOscFig,'Position',[(16-2*NumPlotSlices)/2 0 2*NumPlotSlices 2*2]);
Global.exportToPPTX('addpicture',RBCOscHistFig,'Position',[0 4.5 8 4.5]);
if (HealthyDistPresent)
    [~,pm,~,zm] = ztest(RBCOscMean,HealthyData.HealthyMeans.MeanRBCOsc,HealthyData.HealthyMeans.StdRBCOsc);
    Global.exportToPPTX('addtext',sprintf(['Mean: ',num2str(RBCOscMean,'%2.1f'),' (Healthy = ',num2str(HealthyData.HealthyMeans.MeanRBCOsc,'%2.1f'),'±',num2str(HealthyData.HealthyMeans.StdRBCOsc,'%2.1f'),'); Z = ',num2str(zm,'%1.3f'),', P = ',num2str(pm,'%1.3f')]),'Position',[8 4.5 8 4.5],'HorizontalAlignment','center');
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
    Global.exportToPPTX('addtext',sprintf(['Mean: ',num2str(RBCOscMean,'%2.1f'),'±',num2str(RBCOscStd,'%2.1f')]),'Position',[8 4.5 8 4.5],'HorizontalAlignment','center');
    BinTable   = { ...
    '',{'Defect','BackgroundColor',EightBinMap(1,:),'FontWeight','bold'},{'Low','BackgroundColor',EightBinMap(2,:),'FontWeight','bold'},{'Healthy','BackgroundColor',EightBinMap(3,:),'FontWeight','bold'},{'Healthy','BackgroundColor',EightBinMap(4,:),'FontWeight','bold'},{'Elevated','BackgroundColor',EightBinMap(5,:),'FontWeight','bold'},{'Elevated','BackgroundColor',EightBinMap(6,:),'FontWeight','bold'},{'High','BackgroundColor',EightBinMap(7,:),'FontWeight','bold','Color','w'},{'High','BackgroundColor',EightBinMap(8,:),'FontWeight','bold','Color','w'};...
    'Subject','Bin 1','Bin 2','Bin 3','Bin 4','Bin 5','Bin 6','Bin 7','Bin 8';...
    Subject,num2str(RBCOscBinPercents(1),'%1.2f'),num2str(RBCOscBinPercents(2),'%1.2f'),num2str(RBCOscBinPercents(3),'%1.2f'),num2str(RBCOscBinPercents(4),'%1.2f'),num2str(RBCOscBinPercents(5),'%1.2f'),num2str(RBCOscBinPercents(6),'%1.2f'),num2str(RBCOscBinPercents(7),'%1.2f'),num2str(RBCOscBinPercents(8),'%1.2f');...
    };
end
Global.exportToPPTX('addtable',BinTable,'Position',[8 5 7.75 2],'Vert','middle','Horiz','center','FontSize',12);
Global.exportToPPTX('addtext',sprintf('RBC Oscillation Summary'),'Position',[0 0 5 3]);

Global.exportToPPTX('addslide'); %RBC Osc Workflow
OscWorkFlow_Fig = imread('OscWorkFlow_Fig.png');
Global.exportToPPTX('addpicture',OscWorkFlow_Fig);
Global.exportToPPTX('addtext',sprintf('RBC Oscillation Workflow'),'Position',[0 0 5 3]);

%Save presentation and close presentation -- overwrite file if it already exists
Global.exportToPPTX('save',fullfile(outputpath, ReportTitle));
Global.exportToPPTX('close');
disp('Exporting Results to Powerpoint Summary Completed.')

%% Image quality
fid = GasExchange.CorrectedDissKSpace_SS;
noise_kspace = std(abs(fid(end,end-24:end)));
signal_kspace_end = max(abs(fid(:,end)));
spectralSNR = signal_kspace_end/noise_kspace;

k0 = GasExchange.GasKSpace_SS(1,:);
k0 = smoothdata(k0,"movmedian",20);
ft = fittype('exp1');
opts = fitoptions(ft);
opts.Lower = [0, -Inf];
opts.Upper = [2*max(abs(k0)), 0];
opts.StartPoint = [max(abs(k0)), -0.001];
opts.Weights = normpdf(linspace(0,length(k0),length(k0)),0,length(k0)/4); %weight fit towards start of breath hold to detect change
[~, k0_gof] = fit((1:length(k0))',abs(k0)',ft, opts);
k0_fit = k0_gof.adjrsquare;
% determine quality
if (spectralSNR > 25.0 && k0_fit > 0.99)
    quality = '5-Excellent'; qualityscore = 5;
elseif (spectralSNR > 15.0 && k0_fit > 0.985)
    quality = '4-Good'; qualityscore = 4;
elseif (spectralSNR > 15.0 && k0_fit > 0.98)
    quality = '3-Fair'; qualityscore = 3;
elseif (spectralSNR > 10.0)
    quality = '2-Poor'; qualityscore = 2;
else
    quality = '1-Failed'; qualityscore = 1;
end

%% Save varaibles to GasExchange class

GasExchange.VentThresh = VentThresh;
GasExchange.DissolvedThresh = DissolvedThresh;
GasExchange.BarrierThresh = BarrierThresh;
GasExchange.RBCThresh = RBCThresh;
GasExchange.RBCBarrThresh = RBCBarrThresh;
GasExchange.RBCOscThresh = RBCOscThresh;
GasExchange.HealthyCohortNum = HealthyCohortNum;
GasExchange.VentMean = VentMean;
GasExchange.VentStd = VentStd;
GasExchange.VentBinPercents = VentBinPercents;
GasExchange.DissolvedMean = DissolvedMean;
GasExchange.DissolvedStd = DissolvedStd;
GasExchange.DissolvedBinPercents = DissolvedBinPercents;
GasExchange.BarrierUptakeMean = BarrierUptakeMean;
GasExchange.BarrierUptakeStd = BarrierUptakeStd;
GasExchange.BarrierUptakeBinPercents = BarrierUptakeBinPercents;
GasExchange.RBCTransferMean = RBCTransferMean;
GasExchange.RBCTransferStd = RBCTransferStd;
GasExchange.RBCTransferBinPercents = RBCTransferBinPercents;
GasExchange.RBCBarrierMean = RBCBarrierMean;
GasExchange.RBCBarrierStd = RBCBarrierStd;
GasExchange.RBCBarrierBinPercents = RBCBarrierBinPercents;
GasExchange.RBCOscMean = RBCOscMean;
GasExchange.RBCOscStd = RBCOscStd;
GasExchange.RBCOscBinPercents = RBCOscBinPercents;
GasExchange.imagequality = quality;
GasExchange.image_qualityscore = qualityscore;

%% write report
% ------------------------------------------------------------------------
%  Gas Exchange Analysis Report 
% -------------------------------------------------------------------------
GasExchange.writereport = 'yes'; 
% GasExchange.outputpath
MainInput.ImageQuality = quality;
cd(GasExchange.outputpath)
if strcmp(GasExchange.writereport,'yes')
    clc;
    if ~isfield(MainInput, 'N4Bias') 
        MainInput.N4Bias = 'yes';
    end
    if ~isfield(MainInput, 'AnalysisMethod') 
        MainInput.AnalysisMethod = '1-point Dixon';
    end
    if ~isfield(MainInput, 'AgeCor') 
        MainInput.AgeCor = 'no';
    end
    if ~isfield(MainInput, 'Note') 
        MainInput.Note = '';
    end    
    % DEFINE INPUTS
    pptDir         = GasExchange.outputpath;
    todayStr       = datestr(now,'yyyymmdd');
    if strcmp(MainInput.CCHMC_DbGxAnalysis,'yes')
        pptxFileName = 'GasExchange_Report';
    else
        pptxFileName = ['GasExchange_Report_' todayStr];
    end
    pptxName       = fullfile(pptDir, [pptxFileName,'.pptx']);

    pptoutputPattern = fullfile(pptDir, 'GasExchange_Report_*.pptx');
    pptFiles = dir(pptoutputPattern);   
    for k = 1:length(pptFiles)
        delete(fullfile(pptDir, pptFiles(k).name));  % delete each matching file
    end

    pptoutputPattern = fullfile(pptDir, 'GasExchange_Report_*.pdf');
    pptFiles = dir(pptoutputPattern);   
    for k = 1:length(pptFiles)
        delete(fullfile(pptDir, pptFiles(k).name));  % delete each matching file
    end
    % colored rows
    rowColors = [
        .7 0.7 0.7;     % header
        1 1 1;          % Mean
        0.9 0.4 0.4;    % Defect
        0.9 0.7 0;    % Low
        0.3 0.6 0.9;    % High
    ];
    rowColors2 = [
        .7 0.7 0.7;     % header
        1 1 1;          % Mean
        0.9 0.4 0.4;    % defect M             
        0.9 0.7 0;    % low M
        225/255 129/255 162/255;    % elevalted M            
        197/255 27/255 125/255;    % High M
    ];
    rowColors3 = [
        .7 0.7 0.7;     % header
        1 1 1;          % Mean
        197/255 27/255 125/255;    % high M
        225/255 129/255 162/255;    % elevated M 
        0.9 0.7 0;    % elevalted RBC
        0.9 0.4 0.4;    % High RBC
    ];        
    % Open/Create PPTX
    isOpen = Global.exportToPPTX();
    if ~isempty(isOpen), Global.exportToPPTX('close'); end
    if isfile(pptxName)
        Global.exportToPPTX('open', pptxName);
        Global.exportToPPTX('switchslide', 1);
    else
        Global.exportToPPTX('new','Dimensions',[16 9],'Title','Gas Exchange Analysis','Author','CPIR @ CCHMC');
    end
    
    % Add Slide
    Global.exportToPPTX('addslide');
    
    % Add Title Bar
    Global.exportToPPTX('addtext','Gas Exchange Analysis Report', 'Position',[0 -0.1 5 0.5], 'FontSize',20,'FontWeight','bold','Color',[1 0 0],'BackgroundColor',[1 1 1],'HorizontalAlignment','center');
    
    %Signal Dynamics
    SigDynamics = openfig('SigDynamics.fig');
    Global.exportToPPTX('addpicture',SigDynamics,'Position',[-0.3 6.35 6 2.9]);
    
    % Settings Table
    settinglabels = {'Subject ID','Age(y)','Sex','Disease','Scan Date','Scanner','Scan Software','Sequence','Recon','XIPline Cmt','Denoise','N4Bias','Method','Age Cor.','Image Quality','Note','ProcessDate', 'Analyst Initials'};
    settings = {MainInput.SubjectID,MainInput.Age,MainInput.Sex,MainInput.Disease, MainInput.ScanDate, MainInput.Scanner,MainInput.ScannerSoftware,MainInput.SequenceType,MainInput.Recon,MainInput.AnalysisCode_hash,MainInput.denoiseXe,MainInput.N4Bias,MainInput.AnalysisMethod,MainInput.AgeCor,MainInput.ImageQuality,MainInput.Note,todayStr,MainInput.Analyst};
    settingsssummary = cell(18, 2);
    settingsssummary(1,:) = {
        {'Setting','BackgroundColor',rowColors(1,:),'FontWeight','bold'},
        {'Value','BackgroundColor',rowColors(1,:),'FontWeight','bold'}
    };
    for i = 1:numel(settinglabels)
        settingsssummary{i+1,1} = {settinglabels{i}, 'BackgroundColor', rowColors(2,:), 'FontWeight', 'bold'};
        settingsssummary{i+1,2} = settings{i};
    end
    Global.exportToPPTX('addtable', settingsssummary, 'Position', [0.05 0.35 2.5 3], ...
        'Vert', 'middle', 'Horiz', 'center', 'FontSize', 10);
    
    % images
    Global.exportToPPTX('addpicture',SumVentFig,'Position',[4.8 0 8.5 2.5]);
    Global.exportToPPTX('addpicture',SumBarrFig,'Position',[4.8 2.1 8.5 2.5]);
    Global.exportToPPTX('addpicture',SumRBCFig,'Position',[4.8 4.2 8.5 2.5]);
    Global.exportToPPTX('addpicture',SumRBCBarFig,'Position',[4.8 6.5 8.5 2.5]);
    
    % ---------- helper: build a value cell with conditional background ----------
    baseWhite = [1 1 1];
    lightRed  = [1 0.8 0.8];
    makeValCell = @(val, mu, sd) { val, 'BackgroundColor', ...
        (val < (mu - 1.96*sd) || val > (mu + 1.96*sd)) * lightRed + ...
        ~(val < (mu - 1.96*sd) || val > (mu + 1.96*sd)) * baseWhite };

    makeValCell2 = @(val, mu, sd) { val, 'BackgroundColor', ...
        ( val > (mu + 1.645*sd)) * lightRed + ...
        ~( val > (mu + 1.645*sd)) * baseWhite };
    

    safeHypot = @(a,b) sqrt(a.^2 + b.^2);
    safeHypot2 = @(a, b, c, d) sqrt(a.^2 + b.^2 + c.^2 + d.^2);  

    
    % vent
    % Histogram
    VentHist2Fig = GasExchangeFunctions.CalculateDissolvedHistogram2(...
        ScaledVentImage(ProtonMaskRegistered(:)),...
        VentEdges,VentThresh,SixBinMap,HealthyData.VentEdges,...
        HealthyData.HealthyVentFit,round(VentSNR,2),round(VentMean,2),round(VentMedian,2),round(VentStd,3),VentBinPercents,'Ventilation');
    Global.exportToPPTX('addpicture',VentHist2Fig,'Position',[12.8 0.2 3.2 2]);
            
    % --- Patient values --- VentMean
    valVec = [ round(VentMean,2), ...
               round(VentBinPercents(1),2), ...
               round(VentBinPercents(2),2), ...
               round(VentBinPercents(5)+VentBinPercents(6),2) ];
    
    % --- Reference means ---
    muVec  = [ round(HealthyData.HealthyMeans.MeanVent,2), ...
               round(HealthyData.BinPercentMeans.Vent(1),2), ...
               round(HealthyData.BinPercentMeans.Vent(2),2), ...
               round(HealthyData.BinPercentMeans.Vent(5)+HealthyData.BinPercentMeans.Vent(6),2) ];
    
    % --- Reference SDs (adjust field names if needed) ---
    sdVec  = [ HealthyData.HealthyMeans.StdVent, ...
               HealthyData.BinPercentStds.Vent(1), ...
               HealthyData.BinPercentStds.Vent(2), ...
               safeHypot(HealthyData.BinPercentStds.Vent(5), HealthyData.BinPercentStds.Vent(6)) ];
    
    labels = {'Mean','Defect%','Low%','High%'};
    
    resultssummary = cell(4,3);
    resultssummary(1,:) = {
       {'Metric','BackgroundColor',rowColors(1,:),'FontWeight','bold'}, ...
       {'Value','BackgroundColor',rowColors(1,:),'FontWeight','bold'}, ...
       {'Ref.','BackgroundColor',rowColors(1,:),'FontWeight','bold'}};
    
    for i = 1:numel(labels)
        % Metric label cell (keep your row shading here if you want)
        resultssummary{i+1,1} = {labels{i}, 'BackgroundColor', rowColors(i+1,:), 'FontWeight', 'bold'};
    
        % Value cell with conditional background: white if within LLN–ULN, light red otherwise
        resultssummary{i+1,2} = makeValCell2(valVec(i), muVec(i), sdVec(i));
    
        % Reference column (mean); optionally show LLN–ULN instead:
        % lln = muVec(i)-1.96*sdVec(i); uln = muVec(i)+1.96*sdVec(i);
        % resultssummary{i+1,3} = sprintf('%.2f (%.2f–%.2f)', muVec(i), lln, uln);
        resultssummary{i+1,3} = muVec(i);
    end
    Global.exportToPPTX('addtable', resultssummary, 'Position', [2.6 0.35 2.7 1.5], ...
        'Vert','middle','Horiz','center','FontSize',14);
    
    % Membrane
    % Histogram
    MemHist2Fig = GasExchangeFunctions.CalculateDissolvedHistogram2(...
        BarrGasRatio(VentBinMask(:)),BarEdges,BarrierThresh,...
        EightBinMap,HealthyData.BarsEdges,...
        HealthyData.HealthyBarsFit,round(BarrierSNR,2),round(BarrierUptakeMean,4),...
        round(BarrierUptakeMedian,4),round(BarrierUptakeStd,4),BarrierUptakeBinPercents,'Membrane');
    Global.exportToPPTX('addpicture',MemHist2Fig,'Position',[12.8 2.3 3.2 2.1]);
        
    valVec = [ round(BarrierUptakeMean,4), ...
               round(BarrierUptakeBinPercents(1),2), ...
               round(BarrierUptakeBinPercents(2),2), ...
               round(BarrierUptakeBinPercents(5)+BarrierUptakeBinPercents(6),2),...
               round(BarrierUptakeBinPercents(7)+BarrierUptakeBinPercents(8),2) ];
    
    muVec  = [ round(HealthyData.HealthyMeans.MeanBarrier,4), ...
               round(HealthyData.BinPercentMeans.Barrier(1),2), ...
               round(HealthyData.BinPercentMeans.Barrier(2),2), ...
               round(HealthyData.BinPercentMeans.Barrier(5)+HealthyData.BinPercentMeans.Barrier(6),2),...
               round(HealthyData.BinPercentMeans.Barrier(7)+HealthyData.BinPercentMeans.Barrier(8),2) ];
    
    sdVec  = [ HealthyData.HealthyMeans.StdBarrier, ...
               HealthyData.BinPercentStds.Barrier(1), ...
               HealthyData.BinPercentStds.Barrier(2), ...
               safeHypot(HealthyData.BinPercentStds.Barrier(5), HealthyData.BinPercentStds.Barrier(6)) , ...
               safeHypot(HealthyData.BinPercentStds.Barrier(7), HealthyData.BinPercentStds.Barrier(8)) ];
    
    labels = {'Mean','Defect%','Low%','Elvted%','High%'};
    
    resultssummary = cell(4,3);
    resultssummary(1,:) = {
       {'Metric','BackgroundColor',rowColors2(1,:),'FontWeight','bold'}, ...
       {'Value','BackgroundColor',rowColors2(1,:),'FontWeight','bold'}, ...
       {'Ref.','BackgroundColor',rowColors2(1,:),'FontWeight','bold'}};
    
    for i = 1:numel(labels)
        resultssummary{i+1,1} = {labels{i}, 'BackgroundColor', rowColors2(i+1,:), 'FontWeight', 'bold'};
        resultssummary{i+1,2} = makeValCell(valVec(i), muVec(i), sdVec(i));
        resultssummary{i+1,3} = muVec(i);
    end
    
    Global.exportToPPTX('addtable', resultssummary, 'Position', [2.6 2.45 2.7 1.5], ...
        'Vert','middle','Horiz','center','FontSize',14);
    
    
    % RBC
    % Histogram
    RBCHist2Fig = GasExchangeFunctions.CalculateDissolvedHistogram2(...
        RBCGasRatio(VentBinMask(:)),RBCEdges,RBCThresh,...
        SixBinMap,HealthyData.RBCEdges,HealthyData.HealthyRBCsFit,...
        round(RBCSNR,2),round(RBCTransferMean,4), round(RBCTransferMedian,4),round(RBCTransferStd,4),...
        RBCTransferBinPercents,'RBC');
    Global.exportToPPTX('addpicture',RBCHist2Fig,'Position',[12.8 4.5 3.2 2.2]);
        
    valVec = [ round(RBCTransferMean,4), ...
               round(RBCTransferBinPercents(1),2), ...
               round(RBCTransferBinPercents(2),2), ...
               round(RBCTransferBinPercents(5)+RBCTransferBinPercents(6),2) ];
    
    muVec  = [ round(HealthyData.HealthyMeans.MeanRBC,4), ...
               round(HealthyData.BinPercentMeans.RBC(1),2), ...
               round(HealthyData.BinPercentMeans.RBC(2),2), ...
               round(HealthyData.BinPercentMeans.RBC(5)+HealthyData.BinPercentMeans.RBC(6),2) ];
    
    sdVec  = [ HealthyData.HealthyMeans.StdRBC, ...
               HealthyData.BinPercentStds.RBC(1), ...
               HealthyData.BinPercentStds.RBC(2), ...
               safeHypot(HealthyData.BinPercentStds.RBC(5), HealthyData.BinPercentStds.RBC(6)) ];
    
    labels = {'Mean','Defect%','Low%','High%'};
    
    resultssummary = cell(4,3);
    resultssummary(1,:) = {
       {'Metric','BackgroundColor',rowColors(1,:),'FontWeight','bold'}, ...
       {'Value','BackgroundColor',rowColors(1,:),'FontWeight','bold'}, ...
       {'Ref.','BackgroundColor',rowColors(1,:),'FontWeight','bold'}};
    
    for i = 1:numel(labels)
        resultssummary{i+1,1} = {labels{i}, 'BackgroundColor', rowColors(i+1,:), 'FontWeight', 'bold'};
        resultssummary{i+1,2} = makeValCell(valVec(i), muVec(i), sdVec(i));
        resultssummary{i+1,3} = muVec(i);
    end
    
    Global.exportToPPTX('addtable', resultssummary, 'Position', [2.6 4.55 2.7 1.5], ...
        'Vert','middle','Horiz','center','FontSize',14);
    
    
    % RBC:M
    % Histogram
    RBCBarHist2Fig = GasExchangeFunctions.CalculateDissolvedHistogram2(...
        RBCBarrRatio(VentBinMask(:)),RBCBarEdges,RBCBarrThresh,SixBinRBCBarMap,...
        HealthyData.RBCBarEdges,HealthyData.HealthyRBCBarFit,...
        '',round(RBCBarrierMean,2), round(RBCBarrierMedian,2),round(RBCBarrierStd,2),...
        RBCBarrierBinPercents,'RBC:M');
    Global.exportToPPTX('addpicture',RBCBarHist2Fig,'Position',[12.8 6.7 3.2 2.2]);
        
    valVec = [ round(RBCBarrierMean,2), ...
               round(RBCBarrierBinPercents(1),2), ...
               round(RBCBarrierBinPercents(2),2), ...
               round(RBCBarrierBinPercents(5),2),...
               round(RBCBarrierBinPercents(6),2) ];
    
    muVec  = [ round(HealthyData.HealthyMeans.MeanRBCBar,2), ...
               round(HealthyData.BinPercentMeans.RBCBar(1),2), ...
               round(HealthyData.BinPercentMeans.RBCBar(2),2), ...
               round(HealthyData.BinPercentMeans.RBCBar(5),2),...
               round(HealthyData.BinPercentMeans.RBCBar(6),2) ];
    
    sdVec  = [ HealthyData.HealthyMeans.StdRBCBar, ...
               HealthyData.BinPercentStds.RBCBar(1), ...
               HealthyData.BinPercentStds.RBCBar(2), ...
               HealthyData.BinPercentStds.RBCBar(5), ...
               HealthyData.BinPercentStds.RBCBar(6) ];
    
    labels = {'Mean','Hi-M%','El-M%','El-RBC%','Hi-RBC%'};        
    resultssummary = cell(4,3);
    resultssummary(1,:) = {
       {'Metric','BackgroundColor',rowColors3(1,:),'FontWeight','bold'}, ...
       {'Value','BackgroundColor',rowColors3(1,:),'FontWeight','bold'}, ...
       {'Ref.','BackgroundColor',rowColors3(1,:),'FontWeight','bold'}};
    
    for i = 1:numel(labels)
        resultssummary{i+1,1} = {labels{i}, 'BackgroundColor', rowColors3(i+1,:), 'FontWeight', 'bold'};
        resultssummary{i+1,2} = makeValCell(valVec(i), muVec(i), sdVec(i));
        resultssummary{i+1,3} = muVec(i);
    end
    
    Global.exportToPPTX('addtable', resultssummary, 'Position', [2.6 6.85 2.7 1.5], ...
        'Vert','middle','Horiz','center','FontSize',14);


    % Red box (as an empty text box with red fill/border)
    Global.exportToPPTX('addtext',' ', ...
        'Position',[3.52 2.08 0.3 0.2], ...     % [x y w h] – adjust as needed
        'BackgroundColor',lightRed, ...
        'LineColor',[1 1 1], 'LineWidth',1);
    
    % Label next to it
    Global.exportToPPTX('addtext','Abnormal', ...
        'Position',[3.77 2 1.5 0.3], ...
        'FontSize',14, 'FontWeight','bold', ...
        'Color',[0 0 0], 'HorizontalAlignment','left');


     Global.exportToPPTX('addtext',['SpectSNR = ',num2str(round(spectralSNR,1))], ...
        'Position',[1.17 6.75 1.5 0.3], ...
        'FontSize',10, 'FontWeight','bold', ...
        'Color',[0 0 1], 'HorizontalAlignment','left');   

     Global.exportToPPTX('addtext',['K0-fit = ',num2str(round(k0_fit,3))], ...
        'Position',[1.17 6.95 1.5 0.3], ...
        'FontSize',10, 'FontWeight','bold', ...
        'Color',[0 0 1], 'HorizontalAlignment','left'); 

    % Save & close
    Global.exportToPPTX('save',pptxName);
    Global.exportToPPTX('close');
    fprintf('Report saved to %s\n',pptxName);

     close all;
    
    % save report as a PDF
    ppt = actxserver('PowerPoint.Application');
    
    presentation = ppt.Presentations.Open(pptxName);
    PDFoutputPath = fullfile(pptDir,[pptxFileName,'.pdf']);
    if exist(PDFoutputPath, 'file')
        delete(PDFoutputPath);  % remove existing PDF to avoid overwrite conflict
    end    
    presentation.SaveAs(PDFoutputPath, 32);
    pause(2);  % <-- allow time for file to be written
    presentation.Close();
    ppt.Quit();
    delete(ppt);   
 
end
close all;
%% write JSON file 

if strcmp(MainInput.CCHMC_DbGxAnalysis,'yes')
clc
   scan_num = 1;
   AnalysisCode_path = 'https://github.com/aboodbdaiwi/XIPline';
   AnalysisDate = str2double(datestr(datetime('today'), 'yyyymmdd'));
    % JSON file
    clc;
    headers = {'SUBJECT_ID','Age','Sex','Disease','ScanDate','ScanTime','scan_num','CalibrationDataPath','DissolvedDataPath','AnatomicalDataPath','gx_analysis_folder','Note','ImageQuality','Scanner', ...
               'ScannerSoftware','SequenceType','ReconType','ScanVersion','AnalysisCode_path','AnalysisDate','AnalysisVersion','AnalysisCode_hash',...
               'Denoise','N4Bias','AnalysisMethod','AgeCorrection','Spectral_SNR','k0_fit','TE90','Actual_TE90','Flip_Angle_Factor', ...
               'Package_Version','Process_Date','Healthy_Cohort', ...
               'Gas_Frequency','Barrier_Frequency','RBC_Frequency', ...
               'Gas_T2star','Barrier_T2star','Barrier_LtoG_ratio','RBC_T2star', ...
               'Gas_Phase','Barrier_Phase','RBC_Phase', ...
               'Dixon_Phase_Separation','RBC_Barrier_Ratio','Gas_Contamination','Dixon_Phase_Shift', ...
               'Ventilation_SNR','Gas_SNR','Dissolved_SNR','CorrDissolved_SNR','Barrier_SNR','RBC_SNR', ...
               'Ventilation_Mean','Ventilation_Median','Ventilation_Std_Dev', ...
               'Dissolved_Mean','Dissolved_Median','Dissolved_Std_Dev', ...
               'Barrier_Uptake_Mean','Barrier_Uptake_Median','Barrier_Uptake_Std_Dev', ...
               'RBC_Transfer_Mean','RBC_Transfer_Median','RBC_Transfer_Std_Dev', ...
               'RBC2Barrier_Mean','RBCBarrier_Median','RBC2Barrier_Std_Dev', ...
               'RBC_Osc_K0','RBC_Osc_Mean','RBC_Osc_Std_Dev', ...
               'Ventilation_Bin1_Percent','Ventilation_Bin2_Percent','Ventilation_Bin3_Percent','Ventilation_Bin4_Percent','Ventilation_Bin5_Percent','Ventilation_Bin6_Percent', ...
               'Dissolved_Bin1_Percent','Dissolved_Bin2_Percent','Dissolved_Bin3_Percent','Dissolved_Bin4_Percent','Dissolved_Bin5_Percent','Dissolved_Bin6_Percent', ...
               'Barrier_Uptake_Bin1_Percent','Barrier_Uptake_Bin2_Percent','Barrier_Uptake_Bin3_Percent','Barrier_Uptake_Bin4_Percent','Barrier_Uptake_Bin5_Percent','Barrier_Uptake_Bin6_Percent','Barrier_Uptake_Bin7_Percent','Barrier_Uptake_Bin8_Percent', ...
               'RBC_Transfer_Bin1_Percent','RBC_Transfer_Bin2_Percent','RBC_Transfer_Bin3_Percent','RBC_Transfer_Bin4_Percent','RBC_Transfer_Bin5_Percent','RBC_Transfer_Bin6_Percent', ...
               'RBC2Barrier_Bin1_Percent','RBC2Barrier_Bin2_Percent','RBC2Barrier_Bin3_Percent','RBC2Barrier_Bin4_Percent','RBC2Barrier_Bin5_Percent','RBC2Barrier_Bin6_Percent', ...
               'RBC_Osc_Bin1_Percent','RBC_Osc_Bin2_Percent','RBC_Osc_Bin3_Percent','RBC_Osc_Bin4_Percent','RBC_Osc_Bin5_Percent','RBC_Osc_Bin6_Percent','RBC_Osc_Bin7_Percent','RBC_Osc_Bin8_Percent',...
               'woods_analysis_folder'};
    
    % --- Collect values in same order ---
    NewData = {MainInput.SubjectID, MainInput.Age,MainInput.Sex, MainInput.Disease,MainInput.ScanDate,MainInput.timeStr,scan_num,MainInput.cal_file,MainInput.gx_file, MainInput.anat_file,MainInput.gx_analysis_folder,MainInput.Note,MainInput.ImageQuality, MainInput.Scanner, ...
               MainInput.ScannerSoftware,MainInput.SequenceType,MainInput.Recon,MainInput.Institute,AnalysisCode_path,string(AnalysisDate),MainInput.analysisversion,MainInput.AnalysisCode_hash,...
               MainInput.denoiseXe,MainInput.N4Bias,MainInput.AnalysisMethod,MainInput.AgeCor,spectralSNR,k0_fit, GasExchange.TE90, GasExchange.ActTE90, GasExchange.XeSin.flip_angles.vals(1,1)/DisFlipAngle, ...
               ReconVersion, datestr(date,29), HealthyCohortNum, ...
               GasFreq, BarrierFreq, RBCFreq, ...
               GasT2Star, BarrierT2Star, BarrierLtoGRatio, RBCT2Star, ...
               AppendedDissolvedNMRFit.phase(3), AppendedDissolvedNMRFit.phase(2), AppendedDissolvedNMRFit.phase(1),...%Spectro Fitting Results - Phase
               SpecDissolvedPhaseDifference, SpecRBCBarrierRatio, ...
               SpecGasPhaseContaminationPercentage, Delta_angle_deg, ...
               VentSNR, GasSNR, DissolvedSNR, CorrDissolvedSNR, BarrierSNR, RBCSNR, ...
               VentMean,VentMedian, VentStd, ...
               DissolvedMean, DissolvedMedian, DissolvedStd, ...
               BarrierUptakeMean, BarrierUptakeMedian, BarrierUptakeStd, ...
               RBCTransferMean, RBCTransferMedian, RBCTransferStd, ...
               RBCBarrierMean,RBCBarrierMedian, RBCBarrierStd, ...
               RBC2Bar_struct.RBCk0Osc, RBCOscMean, RBCOscStd, ...
               VentBinPercents(1), VentBinPercents(2), VentBinPercents(3), VentBinPercents(4), VentBinPercents(5), VentBinPercents(6), ...
               DissolvedBinPercents(1), DissolvedBinPercents(2), DissolvedBinPercents(3), DissolvedBinPercents(4), DissolvedBinPercents(5), DissolvedBinPercents(6), ...
               BarrierUptakeBinPercents(1), BarrierUptakeBinPercents(2), BarrierUptakeBinPercents(3), BarrierUptakeBinPercents(4), BarrierUptakeBinPercents(5), BarrierUptakeBinPercents(6), BarrierUptakeBinPercents(7), BarrierUptakeBinPercents(8), ...
               RBCTransferBinPercents(1), RBCTransferBinPercents(2), RBCTransferBinPercents(3), RBCTransferBinPercents(4), RBCTransferBinPercents(5), RBCTransferBinPercents(6), ...
               RBCBarrierBinPercents(1), RBCBarrierBinPercents(2), RBCBarrierBinPercents(3), RBCBarrierBinPercents(4), RBCBarrierBinPercents(5), RBCBarrierBinPercents(6), ...
               RBCOscBinPercents(1), RBCOscBinPercents(2), RBCOscBinPercents(3), RBCOscBinPercents(4), RBCOscBinPercents(5), RBCOscBinPercents(6), RBCOscBinPercents(7), RBCOscBinPercents(8), ''};
    
    % --- Convert into struct ---
    Outputs = cell2struct(NewData, headers, 2);
    OutputJSONFile = fullfile(MainInput.analysisfolder, ['GxAnalysis_','ser-',MainInput.timeStr,'.json']);
    Global.exportStructToJSON(Outputs, OutputJSONFile);

end

%% write out summary as excel files
GasExchangeExcelFile = fullfile(outputpath, 'GasExchange_workspace.xlsx');
Global.exportStructToExcel(GasExchange, GasExchangeExcelFile);
MainInputExcelFile = fullfile(outputpath, 'MainInput_workspace.xlsx');
Global.exportStructToExcel(MainInput, MainInputExcelFile);
ProtonExcelFile = fullfile(outputpath, 'Proton_workspace.xlsx');
Global.exportStructToExcel(Proton, ProtonExcelFile);

%% Export Mat File for Cohort/Additional Analysis
disp('Exporting Workspace...')
close all; % no need to save open figures as they will reopen automatically. can replot if needed
save(fullfile(outputpath, 'GasExchangeAnalysis.mat'),'GasExchange');
save(fullfile(outputpath, 'workspace.mat'));
disp('Exporting Workspace Completed.')

%% Completed
disp('Gas Exchange Imaging Process Completed')
% close(f)
end


% Please add updates here
% 
% 
% 
% 


