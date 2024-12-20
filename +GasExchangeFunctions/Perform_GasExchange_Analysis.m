function Perform_GasExchange_Analysis(...
          VentImage,...
          GasImage,...
          DissolvedImage,...
          CorrDissolvedImage,...
          ProtonImage,...
          H_RecMatrix,...
          LungMask,...
          RBC2Bar_struct,...
          RBCOsc_High_Image,...
          RBCOsc_Low_Image,...
          RBCOsc_Normalization,...
          ProtonImageRegistered,...
          ProtonMaskRegistered,...
          ProtonMax,...
          AppendedDissolvedNMRFit,...
          ActTE90,...
          freq_jump,...
          Slices_Co,...
          Slices_Ax,...
          VentHealthyMean,...
          VentHealthyStd,...
          DissolvedHealthyMean,...
          DissolvedHealthyStd ,...
          BarrierHealthyMean,...
          BarrierHealthyStd ,...
          RBCHealthyMean,...
          RBCtHealthyStd ,...
          RBCBarrHealthyMean ,...
          RBCBarrHealthyStd ,...
          RBCOscHealthyMean ,...
          RBCOscHealthyStd,...
          DataLocation)
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
%% General information
try
    ReconVersion = GasExchange.get_git_hashobject;%use git to find hash of commit used
    ReconVersion = ReconVersion(1:40);%remove enter at end
catch
    ReconVersion = 'TEST';%If not connected to git, can't determine hash so state test
end

NumPlotSlices = 7;%must be odd to work as expected
DisFlipAngle = 20; %Dissolved flip angle
GasFlipAngle = 0.5; %Gas flip angle


%% Import Healthy Distribution, If Available
%Healthy Cohort
if(exist([DataLocation,'\HealthyCohort.mat'],'file') == 2) %if Healthy cohort data exists, use
    HealthyFile = dir([DataLocation,'\HealthyCohort.mat']);%get file info
    HealthyData = load([HealthyFile.folder,'\',HealthyFile.name], '-regexp','^(?!Comb|Ind)...');%import all variable but figures
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
set(BarHistFig,'Name','Barrier Histogram');title(BarHistFig.CurrentAxes,'Barrier:Gas Histogram');ylabel(BarHistFig.CurrentAxes,'Percentage of Voxels');xlabel(BarHistFig.CurrentAxes,'Barrier/Gas Ratio (Barrier-uptake)');
set(RBCHistFig,'Name','RBC Histogram');title(RBCHistFig.CurrentAxes,'RBC:Gas Histogram');ylabel(RBCHistFig.CurrentAxes,'Percentage of Voxels');xlabel(RBCHistFig.CurrentAxes,'RBC/Gas Ratio (RBC-transfer)');
set(RBCBarHistFig,'Name','RBC:Barrier Histogram');title(RBCBarHistFig.CurrentAxes,'RBC:Barrier Histogram');ylabel(RBCBarHistFig.CurrentAxes,'Percentage of Voxels');xlabel(RBCBarHistFig.CurrentAxes,'RBC/Barrier Ratio');
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


%Barrier/Gas Quantification
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
SumBarrFig = figure('Name','Barrier Binned','units','normalized','outerposition',[0 0 1 4/NumPlotSlices]);set(SumBarrFig,'WindowState','minimized');
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
        title('Barrier Binned','FontSize',24)
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
SumRBCFig = figure('Name','RBC Binned','units','normalized','outerposition',[0 0 1 4/NumPlotSlices]);set(SumRBCFig,'WindowState','minimized');
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
        title('RBC Binned','FontSize',24)
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
SumRBCBarFig = figure('Name','RBC:Barrier Ratio Binned','units','normalized','outerposition',[0 0 1 4/NumPlotSlices]);set(SumRBCBarFig,'WindowState','minimized');
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
        title('RBC:Barrier Binned','FontSize',24)
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
ReportTitle = 'GasExchangeAnalysis';
GasExchange_exportToPPTX.exportToPPTX('new','Dimensions',[16 9], ...
    'Title',ReportTitle, ...
    'Author','CPIR @ CCHMC');
Subject='Subject';
%Add slides

GasExchange_exportToPPTX.exportToPPTX('addslide'); %Barrier Image
GasExchange_exportToPPTX.exportToPPTX('addpicture',BarrierMontage);
GasExchange_exportToPPTX.exportToPPTX('addtext',sprintf(['Barrier Image-\nSNR: ', num2str(BarrierSNR)]),'Position',[0 0 5 3]);

GasExchange_exportToPPTX.exportToPPTX('addslide'); %RBC Image
GasExchange_exportToPPTX.exportToPPTX('addpicture',RBCMontage);
GasExchange_exportToPPTX.exportToPPTX('addtext',sprintf(['RBC Image-\nSNR: ', num2str(RBCSNR)]),'Position',[0 0 5 3]);

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
GasExchange_exportToPPTX.exportToPPTX('addtext',sprintf('Barrier-uptake Summary'),'Position',[0 0 5 3]);

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
    GasExchange_exportToPPTX.exportToPPTX('addtext',sprintf(['Mean: ',num2str(RBCBarrierMean,'%1.3f')]),'Position',[8 4.5 8 4.5],'HorizontalAlignment','center');
    BinTable   = { ...
    '',{'Defect','BackgroundColor',SixBinRBCBarMap(1,:),'FontWeight','bold','Color','w'},{'Low','BackgroundColor',SixBinRBCBarMap(2,:),'FontWeight','bold','Color','w'},{'Healthy','BackgroundColor',SixBinRBCBarMap(3,:),'FontWeight','bold'},{'Healthy','BackgroundColor',SixBinRBCBarMap(4,:),'FontWeight','bold'},{'High','BackgroundColor',SixBinRBCBarMap(5,:),'FontWeight','bold'},{'High','BackgroundColor',SixBinRBCBarMap(6,:),'FontWeight','bold'};...
    'Subject','Bin 1','Bin 2','Bin 3','Bin 4','Bin 5','Bin 6';...
    Subject,num2str(RBCBarrierBinPercents(1),'%1.2f'),num2str(RBCBarrierBinPercents(2),'%1.2f'),num2str(RBCBarrierBinPercents(3),'%1.2f'),num2str(RBCBarrierBinPercents(4),'%1.2f'),num2str(RBCBarrierBinPercents(5),'%1.2f'),num2str(RBCBarrierBinPercents(6),'%1.2f');...
    };
end
GasExchange_exportToPPTX.exportToPPTX('addtable',BinTable,'Position',[8 5 7.75 2],'Vert','middle','Horiz','center','FontSize',12);
GasExchange_exportToPPTX.exportToPPTX('addtext',sprintf('RBC:Barrier Summary'),'Position',[0 0 5 3]);

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

%Save presentation and close presentation -- overwrite file if it already exists
GasExchange_exportToPPTX.exportToPPTX('save',fullfile(DataLocation, ReportTitle));
GasExchange_exportToPPTX.exportToPPTX('close');
disp('Exporting Results to Powerpoint Summary Completed.')
%% 






end



