function thresholds = GasExchange_CreateHealthyCohort(DataFolderLocation)
%GasExchange_CreateHealthyCohort - Obatins healthy distributions and thresholds
%   Requires folder to search through and obtain the calculated
%   GasExchangeValues.mat created when anaylzed via
%   GasExchangeImagingProcess for the list in section below
%
%   Syntax:  GasExchange_CreateHealthyCohort(DataLocation)
%
%   Inputs:
%      DataLocation - (optional)string containing the path of the folder
%      containing all subjects with Dissolved Phase Data analysis completed
%
%   Outputs:
%      thresholds - thresholds for the vent, dissolved:gas, barrier-uptake, rbc-transfer,
%      rbc:barrier ratio, and rbc oscillation
%
%   Example: 
%      GasExchange_CreateHealthyCohort()
%      GasExchange_CreateHealthyCohort('C:\Users\mwillmering\Documents\DissolvedData')
%
%   See also: GasExchange_ProcessSingleSubject; GasExchangeV3_ProcessBatch
% 
%   Author: Matthew Willmering
%   Work email: matthew.willmering@cchmc.org
%   Personal email: matt.willmering@gmail.com
%   Website: https://cpir.cchmc.org/

clearvars -except DataFolderLocation; close all;
disp('Gas Exchange Healthy Cohort Analysis Starting...')


%% Healthy Subjects List
[HealthyCohort, CohortTag] = GasExchange.GetHealthyCohort;
HealthyCohortSize = size(HealthyCohort,1);


%% Get Data Folder Location
%Determine folder to use
disp('Looking for necessary subjects...')
if exist('DataFolderLocation','var') == 0
    DataFolderLocation = uigetdir('C:\', 'Select Folder containing Gas Exchange Data Folders');
end


%% Import Data
folders = dir(DataFolderLocation); %obtain directory
numFolders = size(folders,1);
HealthyErrors = 0;

for ii = 1:numFolders
    folder = strsplit(folders(ii).name); %split file name at spaces
    
    if(size(folder,2) ~= 2) %only use folder if date and subject, no notes
        continue;
    else
        subject = folder(1,2);
    end
    
    Healthy = sum(contains(HealthyCohort,subject)); %determine if matched any healthys
    fileName = [folders(ii).folder,'\',folders(ii).name,'\GasExchangeWorkspace.mat'];
    if (Healthy) %import data from folder if healthy
        if(~isfile(fileName))%see if GasExchangeValues file exists
            HealthyErrors = HealthyErrors + 1;
            disp(['Error: GasExchangeWorkspace does not exist for healthy subject ', folders(ii).name])
            continue;
        else
            subject = strrep(subject{1,1},'-','_');
            HealthyGasExchangeStruct.(subject) = load(fileName);
        end
    end
end
NumHealthyData = numel(fieldnames(HealthyGasExchangeStruct));

%% Run Vent First to Remake Vent Mask
%% Organize Data for Easier Processing
HealthyGasExchangeCell = struct2cell(HealthyGasExchangeStruct);
VentHealthy{NumHealthyData,1} = [];
for ii = 1 : NumHealthyData
    VentHealthy{ii,1} = HealthyGasExchangeCell{ii,1}.ScaledVentImage(HealthyGasExchangeCell{ii,1}.ProtonMaskRegistered(:));
end

%% Combine Data
CombVents = cat(1,VentHealthy{:});
%Remove nan's &inf's
removes = isnan(CombVents) | isinf(CombVents);
CombVents(removes) = [];

%% Box-Cox transform all data
%Will use two-parameter model to deal with zeros and negatives
%https://en.wikipedia.org/wiki/Power_transform#Box%E2%80%93Cox_transformation
%requires shifting of data, lambda2 to be min+0.01
CombVentsShifted = CombVents + 0.01 + abs(min(CombVents));
%calculate transforms of shifted data
disp('Calculating Vent Transform started...')
[CombVents_Trans, LambdaVent] = boxcox(abs(CombVentsShifted));
disp('Calculating Vent Transform Completed.')

%% Calculate Mean and Std Dev of Transformed Data
removes = isnan(CombVents_Trans) | isinf(CombVents_Trans);
CombVentMean_Trans = mean(CombVents_Trans(~removes));
CombVentStd_Trans = std(CombVents_Trans(~removes));

%% Transform Thresholds Back
%Vent Thresholds is 6 Bins, uses 5 values mean, mean+-std, mean+-2*std
transThresh = [CombVentMean_Trans-2*CombVentStd_Trans, CombVentMean_Trans-1*CombVentStd_Trans, CombVentMean_Trans, CombVentMean_Trans+1*CombVentStd_Trans, CombVentMean_Trans+2*CombVentStd_Trans];
if LambdaVent ~= 0
    realThresh = (transThresh*LambdaVent+1).^(1/LambdaVent) - (0.01 + abs(min(CombVents)));
else
    realThresh = exp(transThresh) - (0.01 + abs(min(CombVents)));
end
thresholds.Vent = real(realThresh);

%% Calculate Fit to Transformed Data
VentEdges_Trans = linspace(prctile(real(CombVents_Trans(:)),0.01),prctile(real(CombVents_Trans(:)),99.99),75);
HealthyVentFitTrans = normpdf(VentEdges_Trans,CombVentMean_Trans,CombVentStd_Trans); %assume normal distribution now to get fit
TransNormFactor = sum(HealthyVentFitTrans(:));%calculate area under curve
HealthyVentFitTrans = HealthyVentFitTrans./TransNormFactor;%normalize area to one

%% Transform Fit to Real Data
%Can't simply transform fit back since box-cox will be around geometric mean
%Thus, will use a kernel fit as guides to healthy distribution
smoothFactor = 0.2;
VentEdges = linspace(0,1,100);
HealthyVentFit = fitdist(CombVents,'Kernel','Bandwidth',smoothFactor*(thresholds.Vent(4)-thresholds.Vent(2)));
HealthyVentFit = pdf(HealthyVentFit,VentEdges);
NormFactor = sum(HealthyVentFit(:));
HealthyVentFit = HealthyVentFit./NormFactor;

%% Recalculate Ventilation Masks
VentDefectLevel = thresholds.Vent(1);
for ii = 1 : NumHealthyData
    HealthyGasExchangeCell{ii,1}.VentBinMask = logical((HealthyGasExchangeCell{ii,1}.ScaledVentImage > VentDefectLevel).*HealthyGasExchangeCell{ii,1}.ProtonMaskRegistered);
end


%% Run Dissolved, Barrier, and RBC to Remake RBC Mask
%% Organize Data for Easier Processing
DissHealthy{NumHealthyData,1} = [];
BarrierHealthy{NumHealthyData,1} = [];
RBCHealthy{NumHealthyData,1} = [];
RBCBarrHealthy{NumHealthyData,1} = [];
for ii = 1 : NumHealthyData
    DissHealthy{ii,1} = HealthyGasExchangeCell{ii,1}.DissGasRatio(HealthyGasExchangeCell{ii,1}.VentBinMask(:));
    BarrierHealthy{ii,1} = HealthyGasExchangeCell{ii,1}.BarrGasRatio(HealthyGasExchangeCell{ii,1}.VentBinMask(:));
    RBCHealthy{ii,1} = HealthyGasExchangeCell{ii,1}.RBCGasRatio(HealthyGasExchangeCell{ii,1}.VentBinMask(:));
    RBCBarrHealthy{ii,1} = HealthyGasExchangeCell{ii,1}.RBCBarrRatio(HealthyGasExchangeCell{ii,1}.VentBinMask(:));
end

%% Combine Data
CombDiss = cat(1,DissHealthy{:});
CombBars = cat(1,BarrierHealthy{:});
CombRBCs = cat(1,RBCHealthy{:});
CombRBCBars = cat(1,RBCBarrHealthy{:});
%Remove nan's &inf's
removes = isnan(CombDiss) | isinf(CombDiss);
CombDiss(removes) = [];
removes = isnan(CombBars) | isinf(CombBars);
CombBars(removes) = [];
removes = isnan(CombRBCs) | isinf(CombRBCs);
CombRBCs(removes) = [];
removes = isnan(CombRBCBars) | isinf(CombRBCBars);
CombRBCBars(removes) = [];

%% Box-Cox transform all data
%Will use two-parameter model to deal with zeros and negatives
%https://en.wikipedia.org/wiki/Power_transform#Box%E2%80%93Cox_transformation
%requires shifting of data, lambda2 to be min+0.001
CombDissShifted = CombDiss + 0.01 + abs(min(CombDiss));
CombBarsShifted = CombBars + 0.01 + abs(min(CombBars));
CombRBCsShifted = CombRBCs + 0.01 + abs(min(CombRBCs));
CombRBCBarsShifted = CombRBCBars + 0.01 + abs(min(CombRBCBars));
%calculate transforms of shifted data
disp('Calculating Transforms started...')
[CombDiss_Trans, LambdaDiss] = boxcox(abs(CombDissShifted));
[CombBars_Trans, LambdaBar] = boxcox(abs(CombBarsShifted));
[CombRBCs_Trans, LambdaRBC] = boxcox(abs(CombRBCsShifted));
[CombRBCBars_Trans, LambdaRBCBar] = boxcox(abs(CombRBCBarsShifted));
disp('Calculating Transforms Completed.')

%% Calculate Mean and Std Dev of Transformed Data
%Dissolved
removes = isnan(CombDiss_Trans) | isinf(CombDiss_Trans);
CombDissMean_Trans = mean(CombDiss_Trans(~removes));
CombDissStd_Trans = std(CombDiss_Trans(~removes));
%Barrier
removes = isnan(CombBars_Trans) | isinf(CombBars_Trans);
CombBarsMean_Trans = mean(CombBars_Trans(~removes));
CombBarsStd_Trans = std(CombBars_Trans(~removes));
%RBC
removes = isnan(CombRBCs_Trans) | isinf(CombRBCs_Trans);
CombRBCsMean_Trans = mean(CombRBCs_Trans(~removes));
CombRBCsStd_Trans = std(CombRBCs_Trans(~removes));
%RBC:Barrier
removes = isnan(CombRBCBars_Trans) | isinf(CombRBCBars_Trans);
CombRBCBarsMean_Trans = mean(CombRBCBars_Trans(~removes));
CombRBCBarsStd_Trans = std(CombRBCBars_Trans(~removes));

%% Transform Thresholds Back
%Dissolved Thresholds is 6 Bins, uses 5 values mean, mean+-std, mean+-2*std
transThresh = [CombDissMean_Trans-2*CombDissStd_Trans, CombDissMean_Trans-1*CombDissStd_Trans, CombDissMean_Trans, CombDissMean_Trans+1*CombDissStd_Trans, CombDissMean_Trans+2*CombDissStd_Trans];
if LambdaDiss ~= 0
    realThresh = (transThresh*LambdaDiss+1).^(1/LambdaDiss) - (0.01 + abs(min(CombDiss)));
else
    realThresh = exp(transThresh - (0.01 + abs(min(CombDiss))));
end
thresholds.Dissolved = real(realThresh);
%Barrier Thresholds is 8 Bins, uses 7 values mean, mean+-std, mean+-2*std, mean+3std, mean+4std
transThresh = [CombBarsMean_Trans-2*CombBarsStd_Trans, CombBarsMean_Trans-1*CombBarsStd_Trans, CombBarsMean_Trans, CombBarsMean_Trans+1*CombBarsStd_Trans, CombBarsMean_Trans+2*CombBarsStd_Trans, CombBarsMean_Trans+3*CombBarsStd_Trans, CombBarsMean_Trans+4*CombBarsStd_Trans];
if LambdaBar ~= 0
    realThresh = (transThresh*LambdaBar+1).^(1/LambdaBar) - (0.01 + abs(min(CombBars)));
else
    realThresh = exp(transThresh) - (0.01 + abs(min(CombBars)));
end
thresholds.Barrier = real(realThresh);
%RBC Thresholds is 6 Bins, uses 5 values mean, mean+-std, mean+-2*std
transThresh = [CombRBCsMean_Trans-2*CombRBCsStd_Trans, CombRBCsMean_Trans-1*CombRBCsStd_Trans, CombRBCsMean_Trans, CombRBCsMean_Trans+1*CombRBCsStd_Trans, CombRBCsMean_Trans+2*CombRBCsStd_Trans];
if LambdaRBC ~= 0
    realThresh = (transThresh*LambdaRBC+1).^(1/LambdaRBC) - (0.01 + abs(min(CombRBCs)));
else
    realThresh = exp(transThresh) - (0.01 + abs(min(CombRBCs)));
end
thresholds.RBC = real(realThresh);
%RBC:Barr Thresholds is 6 Bins, uses 5 values mean, mean+-std, mean+-2*std
transThresh = [CombRBCBarsMean_Trans-2*CombRBCBarsStd_Trans, CombRBCBarsMean_Trans-1*CombRBCBarsStd_Trans, CombRBCBarsMean_Trans, CombRBCBarsMean_Trans+1*CombRBCBarsStd_Trans, CombRBCBarsMean_Trans+2*CombRBCBarsStd_Trans];
if LambdaRBCBar ~= 0
    realThresh = (transThresh*LambdaRBCBar+1).^(1/LambdaRBCBar) - (0.01 + abs(min(CombRBCBars)));
else
    realThresh = exp(transThresh) - (0.01 + abs(min(CombRBCBars)));
end
thresholds.RBCBarr = real(realThresh);

%% Calculate Fit to Transformed Data
%Dissolved
DissEdges_Trans = linspace(prctile(real(CombDiss_Trans(:)),0.01),prctile(real(CombDiss_Trans(:)),99.99),100);
HealthyDissFitTrans = normpdf(DissEdges_Trans,CombDissMean_Trans,CombDissStd_Trans); %assume normal distribution now to get fit
TransNormFactor = sum(HealthyDissFitTrans(:));%calculate area under curve
HealthyDissFitTrans = HealthyDissFitTrans./TransNormFactor;%normalize area to one
%Barrier
BarsEdges_Trans = linspace(prctile(real(CombBars_Trans(:)),0.1),prctile(real(CombBars_Trans(:)),99.9),100);
HealthyBarsFitTrans = normpdf(BarsEdges_Trans,CombBarsMean_Trans,CombBarsStd_Trans);
TransNormFactor = sum(HealthyBarsFitTrans(:));
HealthyBarsFitTrans = HealthyBarsFitTrans./TransNormFactor;
%RBC
RBCEdges_Trans = linspace(prctile(real(CombRBCs_Trans(:)),0.1),prctile(real(CombRBCs_Trans(:)),99.9),100);
HealthyRBCsFitTrans = normpdf(RBCEdges_Trans,CombRBCsMean_Trans,CombRBCsStd_Trans);
TransNormFactor = sum(HealthyRBCsFitTrans(:));
HealthyRBCsFitTrans = HealthyRBCsFitTrans./TransNormFactor;
%RBC:Barrier
RBCBarsEdges_Trans = linspace(prctile(real(CombRBCBars_Trans(:)),0.1),prctile(real(CombRBCBars_Trans(:)),99.9),100);
HealthyRBCBarFitTrans = normpdf(RBCBarsEdges_Trans,CombRBCBarsMean_Trans,CombRBCBarsStd_Trans);
TransNormFactor = sum(HealthyRBCBarFitTrans(:));
HealthyRBCBarFitTrans = HealthyRBCBarFitTrans./TransNormFactor;

%% Transform Fit to Real Data
%Can't simply transform fit back since box-cox will be around geometric mean
%Thus, will use a kernel fit as guides to healthy distribution
%Dissolved
DissolvedEdges = linspace(0,3*thresholds.Dissolved(3),100);
HealthyDissFit = fitdist(CombDiss,'Kernel','Bandwidth',smoothFactor*(thresholds.Dissolved(4)-thresholds.Dissolved(2)));
HealthyDissFit = pdf(HealthyDissFit,DissolvedEdges);
NormFactor = sum(HealthyDissFit(:));
HealthyDissFit = HealthyDissFit./NormFactor;
%Barrier
BarsEdges = linspace(0,3*thresholds.Barrier(3),100);
HealthyBarsFit = fitdist(CombBars,'Kernel','Bandwidth',smoothFactor*(thresholds.Barrier(5)-thresholds.Barrier(3)));
HealthyBarsFit = pdf(HealthyBarsFit,BarsEdges);
NormFactor = sum(HealthyBarsFit(:));
HealthyBarsFit = HealthyBarsFit./NormFactor;
%RBC
RBCEdges = linspace(0,3*thresholds.RBC(3),100);
HealthyRBCsFit = fitdist(CombRBCs,'Kernel','Bandwidth',smoothFactor*(thresholds.RBC(4)-thresholds.RBC(2)));
HealthyRBCsFit = pdf(HealthyRBCsFit,RBCEdges);
NormFactor = sum(HealthyRBCsFit(:));
HealthyRBCsFit = HealthyRBCsFit./NormFactor;
%RBC:Barrier
RBCBarEdges = linspace(0,3*thresholds.RBCBarr(3),100);
HealthyRBCBarFit = fitdist(CombRBCBars,'Kernel','Bandwidth',smoothFactor*(thresholds.RBCBarr(4)-thresholds.RBCBarr(2)));
HealthyRBCBarFit = pdf(HealthyRBCBarFit,RBCBarEdges);
NormFactor = sum(HealthyRBCBarFit(:));
HealthyRBCBarFit = HealthyRBCBarFit./NormFactor;

%% Recalculate RBC Masks
RBCDefectLevel = thresholds.RBC(1);
for ii = 1 : NumHealthyData
    HealthyGasExchangeCell{ii,1}.RBCBinMask = logical((HealthyGasExchangeCell{ii,1}.RBCGasRatio > RBCDefectLevel).*HealthyGasExchangeCell{ii,1}.VentBinMask);
end


%% Run RBC Oscillation with New RBC Mask
% Organize Data for Easier Processing
RBCOscHealthy{NumHealthyData,1} = [];
for ii = 1 : NumHealthyData
    RBCOscHealthy{ii,1} = HealthyGasExchangeCell{ii,1}.RBCOsc(HealthyGasExchangeCell{ii,1}.RBCBinMask(:));
end

%% Combine Data
CombRBCOsc = cat(1,RBCOscHealthy{:});
%Remove nan's &inf's
removes = isnan(CombRBCOsc) | isinf(CombRBCOsc);
CombVents(removes) = [];

%% Box-Cox transform all data
%Will use two-parameter model to deal with zeros and negatives
%https://en.wikipedia.org/wiki/Power_transform#Box%E2%80%93Cox_transformation
%requires shifting of data, lambda2 to be min+0.001
CombRBCOscShifted = CombRBCOsc + 0.01 + abs(min(CombRBCOsc));
%calculate transforms of shifted data
disp('Calculating Transforms started...')
[CombRBCOsc_Trans, LambdaRBCOsc] = boxcox(abs(CombRBCOscShifted));
disp('Calculating Transforms Completed.')

%% Calculate Mean and Std Dev of Transformed Data
%RBC Osc
removes = isnan(CombRBCOsc_Trans) | isinf(CombRBCOsc_Trans);
CombRBCOscMean_Trans = mean(CombRBCOsc_Trans(~removes));
CombRBCOscStd_Trans = std(CombRBCOsc_Trans(~removes));

%% Transform Thresholds Back
%RBCOsc Thresholds is 8 Bins, uses 7 values mean, mean+-std, mean+-2*std, mean+3*std, mean+4*std
transThresh = [CombRBCOscMean_Trans-2*CombRBCOscStd_Trans, CombRBCOscMean_Trans-1*CombRBCOscStd_Trans, CombRBCOscMean_Trans, CombRBCOscMean_Trans+1*CombRBCOscStd_Trans, CombRBCOscMean_Trans+2*CombRBCOscStd_Trans, CombRBCOscMean_Trans+3*CombRBCOscStd_Trans, CombRBCOscMean_Trans+4*CombRBCOscStd_Trans];
if LambdaRBCOsc ~= 0
    realThresh = (transThresh*LambdaRBCOsc+1).^(1/LambdaRBCOsc) - (0.01 + abs(min(CombRBCOsc)));
else
    realThresh = exp(transThresh) - (0.01 + abs(min(CombRBCOsc)));
end
thresholds.RBCOsc = real(realThresh);

%% Calculate Fit to Transformed Data
%RBC Osc
RBCOscEdges_Trans = linspace(prctile(abs(CombRBCOsc_Trans(:)),0.1),prctile(abs(CombRBCOsc_Trans(:)),99.9),100);
HealthyRBCOscFitTrans = normpdf(RBCOscEdges_Trans,CombRBCOscMean_Trans,CombRBCOscStd_Trans);
TransNormFactor = sum(HealthyRBCOscFitTrans(:));
HealthyRBCOscFitTrans = HealthyRBCOscFitTrans./TransNormFactor;

%% Transform Fit to Real Data
%Can't simply transform fit back since box-cox will be around geometric mean
%Thus, will use a kernel fit as guides to healthy distribution
%RBC Osc
RBCOscEdges = linspace(prctile(CombRBCOsc(:),0.01),100,100);
HealthyRBCOscFit = fitdist(CombRBCOsc,'Kernel','Bandwidth',smoothFactor*(thresholds.RBCOsc(4)-thresholds.RBCOsc(2)));
HealthyRBCOscFit = pdf(HealthyRBCOscFit,RBCOscEdges);
NormFactor = sum(HealthyRBCOscFit(:));
HealthyRBCOscFit = HealthyRBCOscFit./NormFactor;

%% Calculate Cohort Mean and Std
for ii = 1 : NumHealthyData
    %Vent
    tempData = VentHealthy{ii};%remove from struct
    tempData(isinf(tempData)|isnan(tempData)) = [];%delete infs and nans
    HealthyMeans.Vent(ii) = mean(tempData);%get and save mean
    %Dissolved
    tempData = DissHealthy{ii};%remove from struct
    tempData(isinf(tempData)|isnan(tempData)) = [];%delete infs and nans
    HealthyMeans.Dissolved(ii) = mean(tempData);%get and save mean    
    %Barrier
    tempData = BarrierHealthy{ii};%remove from struct
    tempData(isinf(tempData)|isnan(tempData)) = [];%delete infs and nans
    HealthyMeans.Barrier(ii) = mean(tempData);%get and save mean    
    %RBC
    tempData = RBCHealthy{ii};%remove from struct
    tempData(isinf(tempData)|isnan(tempData)) = [];%delete infs and nans
    HealthyMeans.RBC(ii) = mean(tempData);%get and save mean    
    %RBC:bar
    tempData = RBCBarrHealthy{ii};%remove from struct
    tempData(isinf(tempData)|isnan(tempData)) = [];%delete infs and nans
    HealthyMeans.RBCBar(ii) = mean(tempData);%get and save mean    
    %RBC Osc
    tempData = RBCOscHealthy{ii};%remove from struct
    tempData(isinf(tempData)|isnan(tempData)) = [];%delete infs and nans
    HealthyMeans.RBCOsc(ii) = mean(tempData);%get and save mean    
end

HealthyMeans.MeanVent = mean(HealthyMeans.Vent,'all');
HealthyMeans.MeanDissolved = mean(HealthyMeans.Dissolved,'all');
HealthyMeans.MeanBarrier = mean(HealthyMeans.Barrier,'all');
HealthyMeans.MeanRBC = mean(HealthyMeans.RBC,'all');
HealthyMeans.MeanRBCBar = mean(HealthyMeans.RBCBar,'all');
HealthyMeans.MeanRBCOsc = mean(HealthyMeans.RBCOsc,'all');

HealthyMeans.StdVent = std(HealthyMeans.Vent,[],'all');
HealthyMeans.StdDissolved = std(HealthyMeans.Dissolved,[],'all');
HealthyMeans.StdBarrier = std(HealthyMeans.Barrier,[],'all');
HealthyMeans.StdRBC = std(HealthyMeans.RBC,[],'all');
HealthyMeans.StdRBCBar = std(HealthyMeans.RBCBar,[],'all');
HealthyMeans.StdRBCOsc = std(HealthyMeans.RBCOsc,[],'all');


%% Calcuate Bin Means and Std. Devs. for Comparison
for ii = 1 : NumHealthyData
    BinPercents.Vent(ii,1:6) = diff([0, sum(VentHealthy{ii,1}<thresholds.Vent)/length(VentHealthy{ii,1})*100, 100]);
    BinPercents.Dissolved(ii,1:6) = diff([0, sum(DissHealthy{ii,1}<thresholds.Dissolved)/length(DissHealthy{ii,1})*100, 100]);
    BinPercents.Barrier(ii,1:8) = diff([0, sum(BarrierHealthy{ii,1}<thresholds.Barrier)/length(BarrierHealthy{ii,1})*100, 100]);
    BinPercents.RBC(ii,1:6) = diff([0, sum(RBCHealthy{ii,1}<thresholds.RBC)/length(RBCHealthy{ii,1})*100, 100]);
    BinPercents.RBCBar(ii,1:6) = diff([0, sum(RBCBarrHealthy{ii,1}<thresholds.RBCBarr)/length(RBCBarrHealthy{ii,1})*100, 100]);
    BinPercents.RBCOsc(ii,1:8) = diff([0, sum(RBCOscHealthy{ii,1}<thresholds.RBCOsc)/length(RBCOscHealthy{ii,1})*100, 100]);
end

if (NumHealthyData > 1)
    %mean
    BinPercentMeans.Vent = mean(BinPercents.Vent);
    BinPercentMeans.Dissolved = mean(BinPercents.Dissolved);
    BinPercentMeans.Barrier = mean(BinPercents.Barrier);
    BinPercentMeans.RBC = mean(BinPercents.RBC);
    BinPercentMeans.RBCBar = mean(BinPercents.RBCBar);
    BinPercentMeans.RBCOsc = mean(BinPercents.RBCOsc);
    %std
    BinPercentStds.Vent = std(BinPercents.Vent);
    BinPercentStds.Dissolved = std(BinPercents.Dissolved);
    BinPercentStds.Barrier = std(BinPercents.Barrier);
    BinPercentStds.RBC = std(BinPercents.RBC);
    BinPercentStds.RBCBar = std(BinPercents.RBCBar);
    BinPercentStds.RBCOsc = std(BinPercents.RBCOsc);
else
    %mean
    BinPercentMeans.Vent = BinPercents.Vent;
    BinPercentMeans.Dissolved = BinPercents.Dissolved;
    BinPercentMeans.Barrier = BinPercents.Barrier;
    BinPercentMeans.RBC = BinPercents.RBC;
    BinPercentMeans.RBCBar = BinPercents.RBCBar;
    BinPercentMeans.RBCOsc = BinPercents.RBCOsc;
    %std
    BinPercentStds.Vent = zeros(size(BinPercents.Vent));
    BinPercentStds.Dissolved = zeros(size(BinPercents.Dissolved));
    BinPercentStds.Barrier = zeros(size(BinPercents.Barrier));
    BinPercentStds.RBC = zeros(size(BinPercents.RBC));
    BinPercentStds.RBCBar = zeros(size(BinPercents.RBCBar));
    BinPercentStds.RBCOsc = zeros(size(BinPercents.RBCOsc));
end


%% Plot Combined Data w/ Fits
SixBinColor = [1 0 0; 1 0.7143 0; 0.4 0.7 0.4; 0 1 0; 0 0.57 0.71; 0 0 1];
EightBinColor = [1 0 0; 1 0.7143 0; 0.4 0.7 0.4; 0 1 0; 184/255 226/255 145/255; 243/255 205/255 213/255; 225/255 129/255 162/255; 197/255 27/255 125/255];
SixBinRBCBarColor = [197/255 27/255 125/255; 225/255 129/255 162/255; 0.4 0.7 0.4; 0 1 0; 1 0.7143 0; 1 0 0];

%Vent
CombVentHistFig = figure('Name','Ventilation Histogram');
set(CombVentHistFig,'color','white','Units','inches','Position',[0.25 0.25 8 9])
tiledlayout(2,1,'TileSpacing','none')
%[ha, ~] = tight_subplot(2, 1, [0.1 0.05], [0.05 0.025], [0.05 0.01]);
    %Plot Standard
%axes(ha(1));
nexttile
hold on
histogram(CombVents,VentEdges,'Normalization','probability','EdgeColor',[0.5 0.5 0.5],'FaceColor',[1 1 1],'FaceAlpha',0.75);%histogram
plot(VentEdges,abs(HealthyVentFit),'k--');%healthy fit
legend('Subjects Distribution','Healthy Control Distribution','Location','best');
    %bins
ylims = ylim;
xlims = xlim;
pat = patch([0 0 thresholds.Vent(1,1) thresholds.Vent(1,1)],[0 ylims(2) ylims(2) 0],SixBinColor(1,:),'HandleVisibility','off');
pat.LineStyle = 'none';
pat.FaceVertexAlphaData = 0.1;
pat.FaceAlpha = 'flat';
uistack(pat,'bottom');
for bin = 2:size(thresholds.Vent,2)
    pat = patch([thresholds.Vent(1,bin-1) thresholds.Vent(1,bin-1) thresholds.Vent(1,bin) thresholds.Vent(1,bin)],[0 ylims(2) ylims(2) 0],SixBinColor(bin,:),'HandleVisibility','off');
    pat.LineStyle = 'none';
    pat.FaceVertexAlphaData = 0.1;
    pat.FaceAlpha = 'flat';
    uistack(pat,'bottom');
end
pat = patch([thresholds.Vent(1,5) thresholds.Vent(1,5) xlims(2) xlims(2)],[0 ylims(2) ylims(2) 0],SixBinColor(6,:),'HandleVisibility','off');
pat.LineStyle = 'none';
pat.FaceVertexAlphaData = 0.1;
pat.FaceAlpha = 'flat';
uistack(pat,'bottom');
axis([0 1 -inf inf])
title('Healthy Ventilation Cohort')
hold off
    %Plot Transformed
%axes(ha(2));
nexttile
hold on
histogram(CombVents_Trans,VentEdges_Trans,'Normalization','probability','EdgeColor',[0.5 0.5 0.5],'FaceColor',[1 1 1],'FaceAlpha',0.75);
plot(VentEdges_Trans,HealthyVentFitTrans,'k--');
legend('Subjects Distribution','Healthy Control Distribution','Location','best');
title('Healthy Ventilation Cohort - Transformed')
hold off

%Dissolved
CombDissHistFig = figure('Name','Dissolved Histogram');
set(CombDissHistFig,'color','white','Units','inches','Position',[0.25 0.25 8 9])
tiledlayout(2,1,'TileSpacing','none')
    %Plot Standard
nexttile
hold on
histogram(CombDiss,DissolvedEdges,'Normalization','probability','EdgeColor',[0.5 0.5 0.5],'FaceColor',[1 1 1],'FaceAlpha',0.75);%histogram
plot(DissolvedEdges,abs(HealthyDissFit),'k--');%healthy fit
legend('Subjects Distribution','Healthy Control Distribution','Location','best');
    %bins
ylims = ylim;
xlims = xlim;
pat = patch([0 0 thresholds.Dissolved(1,1) thresholds.Dissolved(1,1)],[0 ylims(2) ylims(2) 0],SixBinColor(1,:),'HandleVisibility','off');
pat.LineStyle = 'none';
pat.FaceVertexAlphaData = 0.1;
pat.FaceAlpha = 'flat';
uistack(pat,'bottom');
for bin = 2:size(thresholds.Dissolved,2)
    pat = patch([thresholds.Dissolved(1,bin-1) thresholds.Dissolved(1,bin-1) thresholds.Dissolved(1,bin) thresholds.Dissolved(1,bin)],[0 ylims(2) ylims(2) 0],SixBinColor(bin,:),'HandleVisibility','off');
    pat.LineStyle = 'none';
    pat.FaceVertexAlphaData = 0.1;
    pat.FaceAlpha = 'flat';
    uistack(pat,'bottom');
end
pat = patch([thresholds.Dissolved(1,5) thresholds.Dissolved(1,5) xlims(2) xlims(2)],[0 ylims(2) ylims(2) 0],SixBinColor(6,:),'HandleVisibility','off');
pat.LineStyle = 'none';
pat.FaceVertexAlphaData = 0.1;
pat.FaceAlpha = 'flat';
uistack(pat,'bottom');
axis([0 3*thresholds.Dissolved(3) -inf inf])
title('Healthy Dissolved Cohort')
hold off
    %Plot Transformed
nexttile
hold on
histogram(CombDiss_Trans,DissEdges_Trans,'Normalization','probability','EdgeColor',[0.5 0.5 0.5],'FaceColor',[1 1 1],'FaceAlpha',0.75);
plot(DissEdges_Trans,HealthyDissFitTrans,'k--');
legend('Subjects Distribution','Healthy Control Distribution','Location','best');
title('Healthy Dissolved Cohort - Transformed')
hold off

%Barrier
CombBarsHistFig = figure('Name','Barrier Histogram');
set(CombBarsHistFig,'color','white','Units','inches','Position',[0.25 0.25 8 9])
tiledlayout(2,1,'TileSpacing','none')
    %Plot Standard
nexttile
hold on
histogram(CombBars,BarsEdges,'Normalization','probability','EdgeColor',[0.5 0.5 0.5],'FaceColor',[1 1 1],'FaceAlpha',0.75);%histogram
plot(BarsEdges,abs(HealthyBarsFit),'k--');%healthy fit
legend('Subjects Distribution','Healthy Control Distribution','Location','best');
    %bins
ylims = ylim;
xlims = [0 3*thresholds.Barrier(3)];
pat = patch([xlims(1) xlims(1) thresholds.Barrier(1,1) thresholds.Barrier(1,1)],[0 ylims(2) ylims(2) 0],EightBinColor(1,:),'HandleVisibility','off');
pat.LineStyle = 'none';
pat.FaceVertexAlphaData = 0.1;
pat.FaceAlpha = 'flat';
uistack(pat,'bottom');
for bin = 2:size(thresholds.Barrier,2)
    pat = patch([thresholds.Barrier(1,bin-1) thresholds.Barrier(1,bin-1) thresholds.Barrier(1,bin) thresholds.Barrier(1,bin)],[0 ylims(2) ylims(2) 0],EightBinColor(bin,:),'HandleVisibility','off');
    pat.LineStyle = 'none';
    pat.FaceVertexAlphaData = 0.1;
    pat.FaceAlpha = 'flat';
    uistack(pat,'bottom');
end
pat = patch([thresholds.Barrier(1,7) thresholds.Barrier(1,7) xlims(2) xlims(2)],[0 ylims(2) ylims(2) 0],EightBinColor(8,:),'HandleVisibility','off');
pat.LineStyle = 'none';
pat.FaceVertexAlphaData = 0.1;
pat.FaceAlpha = 'flat';
uistack(pat,'bottom');
axis([0 3*thresholds.Barrier(3) -inf inf])
title('Healthy Barrier Cohort')
hold off
    %Plot Transformed
nexttile
hold on
histogram(CombBars_Trans,BarsEdges_Trans,'Normalization','probability','EdgeColor',[0.5 0.5 0.5],'FaceColor',[1 1 1],'FaceAlpha',0.75);
plot(BarsEdges_Trans,HealthyBarsFitTrans,'k--');
legend('Subjects Distribution','Healthy Control Distribution','Location','best');
title('Healthy Barrier Cohort - Transformed')
hold off

%RBC
CombRBCsHistFig = figure('Name','RBC Histogram');
set(CombRBCsHistFig,'color','white','Units','inches','Position',[0.25 0.25 8 9])
tiledlayout(2,1,'TileSpacing','none')
    %Plot Standard
nexttile
hold on
histogram(CombRBCs,RBCEdges,'Normalization','probability','EdgeColor',[0.5 0.5 0.5],'FaceColor',[1 1 1],'FaceAlpha',0.75);%histogram
plot(RBCEdges,abs(HealthyRBCsFit),'k--');%healthy fit
legend('Subjects Distribution','Healthy Control Distribution','Location','best');
    %bins
ylims = ylim;
xlims = [0 3*thresholds.RBC(3)];
pat = patch([xlims(1) xlims(1) thresholds.RBC(1,1) thresholds.RBC(1,1)],[0 ylims(2) ylims(2) 0],SixBinColor(1,:),'HandleVisibility','off');
pat.LineStyle = 'none';
pat.FaceVertexAlphaData = 0.1;
pat.FaceAlpha = 'flat';
uistack(pat,'bottom');
for bin = 2:size(thresholds.RBC,2)
    pat = patch([thresholds.RBC(1,bin-1) thresholds.RBC(1,bin-1) thresholds.RBC(1,bin) thresholds.RBC(1,bin)],[0 ylims(2) ylims(2) 0],SixBinColor(bin,:),'HandleVisibility','off');
    pat.LineStyle = 'none';
    pat.FaceVertexAlphaData = 0.1;
    pat.FaceAlpha = 'flat';
    uistack(pat,'bottom');
end
pat = patch([thresholds.RBC(1,5) thresholds.RBC(1,5) xlims(2) xlims(2)],[0 ylims(2) ylims(2) 0],SixBinColor(6,:),'HandleVisibility','off');
pat.LineStyle = 'none';
pat.FaceVertexAlphaData = 0.1;
pat.FaceAlpha = 'flat';
uistack(pat,'bottom');
axis([0 3*thresholds.RBC(3) -inf inf])
title('Healthy RBC Cohort')
hold off
    %Plot Transformed
nexttile
hold on
histogram(CombRBCs_Trans,RBCEdges_Trans,'Normalization','probability','EdgeColor',[0.5 0.5 0.5],'FaceColor',[1 1 1],'FaceAlpha',0.75);
plot(RBCEdges_Trans,HealthyRBCsFitTrans,'k--');
legend('Subjects Distribution','Healthy Control Distribution','Location','best');
title('Healthy RBC Cohort - Transformed')
hold off

%RBC:Barrier
CombRBCBarHistFig = figure('Name','RBC:Barrier Histogram');
set(CombRBCBarHistFig,'color','white','Units','inches','Position',[0.25 0.25 8 9])
tiledlayout(2,1,'TileSpacing','none')
    %Plot Standard
nexttile
hold on
histogram(CombRBCBars,RBCBarEdges,'Normalization','probability','EdgeColor',[0.5 0.5 0.5],'FaceColor',[1 1 1],'FaceAlpha',0.75);%histogram
plot(RBCBarEdges,abs(HealthyRBCBarFit),'k--');%healthy fit
legend('Subjects Distribution','Healthy Control Distribution','Location','best');
    %bins
ylims = ylim;
xlims = [0 3*thresholds.RBCBarr(3)];
pat = patch([xlims(1) xlims(1) thresholds.RBCBarr(1,1) thresholds.RBCBarr(1,1)],[0 ylims(2) ylims(2) 0],SixBinRBCBarColor(1,:),'HandleVisibility','off');
pat.LineStyle = 'none';
pat.FaceVertexAlphaData = 0.1;
pat.FaceAlpha = 'flat';
uistack(pat,'bottom');
for bin = 2:size(thresholds.RBCBarr,2)
    pat = patch([thresholds.RBCBarr(1,bin-1) thresholds.RBCBarr(1,bin-1) thresholds.RBCBarr(1,bin) thresholds.RBCBarr(1,bin)],[0 ylims(2) ylims(2) 0],SixBinRBCBarColor(bin,:),'HandleVisibility','off');
    pat.LineStyle = 'none';
    pat.FaceVertexAlphaData = 0.1;
    pat.FaceAlpha = 'flat';
    uistack(pat,'bottom');
end
pat = patch([thresholds.RBCBarr(1,5) thresholds.RBCBarr(1,5) xlims(2) xlims(2)],[0 ylims(2) ylims(2) 0],SixBinRBCBarColor(6,:),'HandleVisibility','off');
pat.LineStyle = 'none';
pat.FaceVertexAlphaData = 0.1;
pat.FaceAlpha = 'flat';
uistack(pat,'bottom');
axis([0 3*thresholds.RBCBarr(3) -inf inf])
title('Healthy RBC:Barrier Cohort')
hold off
    %Plot Transformed
nexttile
hold on
histogram(real(CombRBCBars_Trans),real(RBCBarsEdges_Trans),'Normalization','probability','EdgeColor',[0.5 0.5 0.5],'FaceColor',[1 1 1],'FaceAlpha',0.75);
plot(real(RBCBarsEdges_Trans),real(HealthyRBCBarFitTrans),'k--');
legend('Subjects Distribution','Healthy Control Distribution','Location','best');
title('Healthy RBC:Barrier Cohort - Transformed')
hold off

%RBC Osc
CombRBCOscHistFig = figure('Name','RBC Oscillation Histogram');
set(CombRBCOscHistFig,'color','white','Units','inches','Position',[0.25 0.25 8 9])
tiledlayout(2,1,'TileSpacing','none')
    %Plot Standard
nexttile
hold on
histogram(CombRBCOsc,RBCOscEdges,'Normalization','probability','EdgeColor',[0.5 0.5 0.5],'FaceColor',[1 1 1],'FaceAlpha',0.75);%histogram
plot(RBCOscEdges,abs(HealthyRBCOscFit),'k--');%healthy fit
legend('Subjects Distribution','Healthy Control Distribution','Location','best');
    %bins
ylims = ylim;
xlims = [-100 100];
pat = patch([xlims(1) xlims(1) thresholds.RBCOsc(1,1) thresholds.RBCOsc(1,1)],[0 ylims(2) ylims(2) 0],EightBinColor(1,:),'HandleVisibility','off');
pat.LineStyle = 'none';
pat.FaceVertexAlphaData = 0.1;
pat.FaceAlpha = 'flat';
uistack(pat,'bottom');
for bin = 2:size(thresholds.RBCOsc,2)
    pat = patch([thresholds.RBCOsc(1,bin-1) thresholds.RBCOsc(1,bin-1) thresholds.RBCOsc(1,bin) thresholds.RBCOsc(1,bin)],[0 ylims(2) ylims(2) 0],EightBinColor(bin,:),'HandleVisibility','off');
    pat.LineStyle = 'none';
    pat.FaceVertexAlphaData = 0.1;
    pat.FaceAlpha = 'flat';
    uistack(pat,'bottom');
end
pat = patch([thresholds.RBCOsc(1,7) thresholds.RBCOsc(1,7) xlims(2) xlims(2)],[0 ylims(2) ylims(2) 0],EightBinColor(8,:),'HandleVisibility','off');
pat.LineStyle = 'none';
pat.FaceVertexAlphaData = 0.1;
pat.FaceAlpha = 'flat';
uistack(pat,'bottom');
axis([-100 100 -inf inf])
title('Healthy RBC Oscillation Cohort')
hold off
    %Plot Transformed
nexttile
hold on
histogram(CombRBCOsc_Trans,RBCOscEdges_Trans,'Normalization','probability','EdgeColor',[0.5 0.5 0.5],'FaceColor',[1 1 1],'FaceAlpha',0.75);
plot(RBCOscEdges_Trans,HealthyRBCOscFitTrans,'k--');
legend('Subjects Distribution','Healthy Control Distribution','Location','best');
title('Healthy RBC Oscillation Cohort - Transformed')
hold off


%% Plot Individual Data w/ Fits
%Vent
IndVentHistFig = figure('Name','Ventilation Histogram');
set(IndVentHistFig,'color','white','Units','inches','Position',[0.25 0.25 8 4.5])
    %Plot Standard
hold on
for ii = 1:NumHealthyData
    histogram(VentHealthy{ii,:},VentEdges,'Normalization','probability','EdgeColor','auto','Linewidth',1.5,'DisplayStyle','stairs');%histogram
end
plot(VentEdges,abs(HealthyVentFit),'k--','Linewidth',1.5);%healthy fit
legend(fieldnames(HealthyGasExchangeStruct),'Location','best','Interpreter','none')
    %bins
ylims = ylim;
xlims = xlim;
pat = patch([0 0 thresholds.Vent(1,1) thresholds.Vent(1,1)],[0 ylims(2) ylims(2) 0],SixBinColor(1,:),'HandleVisibility','off','FaceAlpha',0.5);
pat.LineStyle = 'none';
pat.FaceVertexAlphaData = 0.1;
pat.FaceAlpha = 'flat';
uistack(pat,'bottom');
for bin = 2:size(thresholds.Vent,2)
    pat = patch([thresholds.Vent(1,bin-1) thresholds.Vent(1,bin-1) thresholds.Vent(1,bin) thresholds.Vent(1,bin)],[0 ylims(2) ylims(2) 0],SixBinColor(bin,:),'HandleVisibility','off','FaceAlpha',0.75);
    pat.LineStyle = 'none';
    pat.FaceVertexAlphaData = 0.1;
    pat.FaceAlpha = 'flat';
    uistack(pat,'bottom');
end
pat = patch([thresholds.Vent(1,5) thresholds.Vent(1,5) xlims(2) xlims(2)],[0 ylims(2) ylims(2) 0],SixBinColor(6,:),'HandleVisibility','off','FaceAlpha',0.5);
pat.LineStyle = 'none';
pat.FaceVertexAlphaData = 0.1;
pat.FaceAlpha = 'flat';
uistack(pat,'bottom');
axis([0 1 -inf inf])
alim([0 1])
title('Healthy Ventilation Cohort')
hold off

%Dissolved
IndDissHistFig = figure('Name','Dissolved Histogram');
set(IndDissHistFig,'color','white','Units','inches','Position',[0.25 0.25 8 4.5])
    %Plot Standard
hold on
for ii = 1:NumHealthyData
    histogram(DissHealthy{ii,:},DissolvedEdges,'Normalization','probability','EdgeColor','auto','Linewidth',1.5,'DisplayStyle','stairs');%histogram
end
plot(DissolvedEdges,abs(HealthyDissFit),'k--','Linewidth',1.5);%healthy fit
legend(fieldnames(HealthyGasExchangeStruct),'Location','best','Interpreter','none')
    %bins
ylims = ylim;
xlims = xlim;
pat = patch([0 0 thresholds.Dissolved(1,1) thresholds.Dissolved(1,1)],[0 ylims(2) ylims(2) 0],SixBinColor(1,:),'HandleVisibility','off','FaceAlpha',0.5);
pat.LineStyle = 'none';
pat.FaceVertexAlphaData = 0.1;
pat.FaceAlpha = 'flat';
uistack(pat,'bottom');
for bin = 2:size(thresholds.Dissolved,2)
    pat = patch([thresholds.Dissolved(1,bin-1) thresholds.Dissolved(1,bin-1) thresholds.Dissolved(1,bin) thresholds.Dissolved(1,bin)],[0 ylims(2) ylims(2) 0],SixBinColor(bin,:),'HandleVisibility','off','FaceAlpha',0.75);
    pat.LineStyle = 'none';
    pat.FaceVertexAlphaData = 0.1;
    pat.FaceAlpha = 'flat';
    uistack(pat,'bottom');
end
pat = patch([thresholds.Dissolved(1,5) thresholds.Dissolved(1,5) xlims(2) xlims(2)],[0 ylims(2) ylims(2) 0],SixBinColor(6,:),'HandleVisibility','off','FaceAlpha',0.5);
pat.LineStyle = 'none';
pat.FaceVertexAlphaData = 0.1;
pat.FaceAlpha = 'flat';
uistack(pat,'bottom');
axis([0 3*thresholds.Dissolved(3) -inf inf])
alim([0 1])
title('Healthy Dissolved Cohort')
hold off

%Barrier
IndBarsHistFig = figure('Name','Barrier Histogram');
set(IndBarsHistFig,'color','white','Units','inches','Position',[0.25 0.25 8 4.5])
    %Plot Standard
hold on
for ii = 1:NumHealthyData
    histogram(BarrierHealthy{ii,:},BarsEdges,'Normalization','probability','EdgeColor','auto','Linewidth',1.5,'DisplayStyle','stairs');%histogram
end
plot(BarsEdges,abs(HealthyBarsFit),'k--','Linewidth',1.5);%healthy fit
legend(fieldnames(HealthyGasExchangeStruct),'Location','best','Interpreter','none')
    %bins
ylims = ylim;
xlims = [0 3*thresholds.Barrier(4)];
pat = patch([xlims(1) xlims(1) thresholds.Barrier(1,1) thresholds.Barrier(1,1)],[0 ylims(2) ylims(2) 0],EightBinColor(1,:),'HandleVisibility','off','FaceAlpha',1,'AlphaDataMapping','none');
pat.LineStyle = 'none';
pat.FaceVertexAlphaData = 0.1;
pat.FaceAlpha = 'flat';
uistack(pat,'bottom');
for bin = 2:size(thresholds.Barrier,2)
    pat = patch([thresholds.Barrier(1,bin-1) thresholds.Barrier(1,bin-1) thresholds.Barrier(1,bin) thresholds.Barrier(1,bin)],[0 ylims(2) ylims(2) 0],EightBinColor(bin,:),'HandleVisibility','off','FaceAlpha',1,'AlphaDataMapping','none');
    pat.LineStyle = 'none';
    pat.FaceVertexAlphaData = 0.1;
    pat.FaceAlpha = 'flat';
    uistack(pat,'bottom');
end
pat = patch([thresholds.Barrier(1,7) thresholds.Barrier(1,7) xlims(2) xlims(2)],[0 ylims(2) ylims(2) 0],EightBinColor(8,:),'HandleVisibility','off','FaceAlpha',1,'AlphaDataMapping','none');
pat.LineStyle = 'none';
pat.FaceVertexAlphaData = 0.1;
pat.FaceAlpha = 'flat';
uistack(pat,'bottom');
axis([0 3*thresholds.Barrier(3) -inf inf])
title('Healthy Barrier Cohort')
hold off

%RBC
IndRBCsHistFig = figure('Name','RBC Histogram');
set(IndRBCsHistFig,'color','white','Units','inches','Position',[0.25 0.25 8 4.5])
%Plot Standard
hold on
for ii = 1:NumHealthyData
    histogram(RBCHealthy{ii,:},RBCEdges,'Normalization','probability','EdgeColor','auto','Linewidth',1.5,'DisplayStyle','stairs');%histogram
end
plot(RBCEdges,abs(HealthyRBCsFit),'k--','Linewidth',1.5);%healthy fit
legend(fieldnames(HealthyGasExchangeStruct),'Location','best','Interpreter','none')
%bins
ylims = ylim;
xlims = [0 3*thresholds.RBC(3)];
pat = patch([xlims(1) xlims(1) thresholds.RBC(1,1) thresholds.RBC(1,1)],[0 ylims(2) ylims(2) 0],SixBinColor(1,:),'HandleVisibility','off','FaceAlpha',1,'AlphaDataMapping','none');
pat.LineStyle = 'none';
pat.FaceVertexAlphaData = 0.1;
pat.FaceAlpha = 'flat';
uistack(pat,'bottom');
for bin = 2:size(thresholds.RBC,2)
    pat = patch([thresholds.RBC(1,bin-1) thresholds.RBC(1,bin-1) thresholds.RBC(1,bin) thresholds.RBC(1,bin)],[0 ylims(2) ylims(2) 0],SixBinColor(bin,:),'HandleVisibility','off','FaceAlpha',1,'AlphaDataMapping','none');
    pat.LineStyle = 'none';
    pat.FaceVertexAlphaData = 0.1;
    pat.FaceAlpha = 'flat';
    uistack(pat,'bottom');
end
pat = patch([thresholds.RBC(1,5) thresholds.RBC(1,5) xlims(2) xlims(2)],[0 ylims(2) ylims(2) 0],SixBinColor(6,:),'HandleVisibility','off','FaceAlpha',1,'AlphaDataMapping','none');
pat.LineStyle = 'none';
pat.FaceVertexAlphaData = 0.1;
pat.FaceAlpha = 'flat';
uistack(pat,'bottom');
axis([0 3*thresholds.RBC(3) -inf inf])
title('Healthy RBC Cohort')
hold off

%RBC:Barrier
IndRBC_BarsHistFig = figure('Name','RBC:Barrier Histogram');
set(IndRBC_BarsHistFig,'color','white','Units','inches','Position',[0.25 0.25 8 4.5])
%Plot Standard
hold on
for ii = 1:NumHealthyData
    histogram(RBCBarrHealthy{ii,:},RBCBarEdges,'Normalization','probability','EdgeColor','auto','Linewidth',1.5,'DisplayStyle','stairs');%histogram
end
plot(RBCBarEdges,abs(HealthyRBCBarFit),'k--','Linewidth',1.5);%healthy fit
legend(fieldnames(HealthyGasExchangeStruct),'Location','best','Interpreter','none')
%bins
ylims = ylim;
xlims = [0 3*thresholds.RBCBarr(3)];
pat = patch([xlims(1) xlims(1) thresholds.RBCBarr(1,1) thresholds.RBCBarr(1,1)],[0 ylims(2) ylims(2) 0],SixBinRBCBarColor(1,:),'HandleVisibility','off','FaceAlpha',1,'AlphaDataMapping','none');
pat.LineStyle = 'none';
pat.FaceVertexAlphaData = 0.1;
pat.FaceAlpha = 'flat';
uistack(pat,'bottom');
for bin = 2:size(thresholds.RBCBarr,2)
    pat = patch([thresholds.RBCBarr(1,bin-1) thresholds.RBCBarr(1,bin-1) thresholds.RBCBarr(1,bin) thresholds.RBCBarr(1,bin)],[0 ylims(2) ylims(2) 0],SixBinRBCBarColor(bin,:),'HandleVisibility','off','FaceAlpha',1,'AlphaDataMapping','none');
    pat.LineStyle = 'none';
    pat.FaceVertexAlphaData = 0.1;
    pat.FaceAlpha = 'flat';
    uistack(pat,'bottom');
end
pat = patch([thresholds.RBCBarr(1,5) thresholds.RBCBarr(1,5) xlims(2) xlims(2)],[0 ylims(2) ylims(2) 0],SixBinRBCBarColor(6,:),'HandleVisibility','off','FaceAlpha',1,'AlphaDataMapping','none');
pat.LineStyle = 'none';
pat.FaceVertexAlphaData = 0.1;
pat.FaceAlpha = 'flat';
uistack(pat,'bottom');
axis([0 3*thresholds.RBCBarr(3) -inf inf])
title('Healthy RBC:Barrier Cohort')
hold off

%RBC Osc
IndRBCOscHistFig = figure('Name','RBC Oscillation Histogram');
set(IndRBCOscHistFig,'color','white','Units','inches','Position',[0.25 0.25 8 4.5])
    %Plot Standard
hold on
for ii = 1:NumHealthyData
    histogram(RBCOscHealthy{ii,:},RBCOscEdges,'Normalization','probability','EdgeColor','auto','Linewidth',1.5,'DisplayStyle','stairs');%histogram
end
plot(RBCOscEdges,abs(HealthyRBCOscFit),'k--','Linewidth',1.5);%healthy fit
legend(fieldnames(HealthyGasExchangeStruct),'Location','best','Interpreter','none')
    %bins
ylims = ylim;
xlims = [-100 100];
pat = patch([xlims(1) xlims(1) thresholds.RBCOsc(1,1) thresholds.RBCOsc(1,1)],[0 ylims(2) ylims(2) 0],EightBinColor(1,:),'HandleVisibility','off','FaceAlpha',1,'AlphaDataMapping','none');
pat.LineStyle = 'none';
pat.FaceVertexAlphaData = 0.1;
pat.FaceAlpha = 'flat';
uistack(pat,'bottom');
for bin = 2:size(thresholds.RBCOsc,2)
    pat = patch([thresholds.RBCOsc(1,bin-1) thresholds.RBCOsc(1,bin-1) thresholds.RBCOsc(1,bin) thresholds.RBCOsc(1,bin)],[0 ylims(2) ylims(2) 0],EightBinColor(bin,:),'HandleVisibility','off','FaceAlpha',1,'AlphaDataMapping','none');
    pat.LineStyle = 'none';
    pat.FaceVertexAlphaData = 0.1;
    pat.FaceAlpha = 'flat';
    uistack(pat,'bottom');
end
pat = patch([thresholds.RBCOsc(1,7) thresholds.RBCOsc(1,7) xlims(2) xlims(2)],[0 ylims(2) ylims(2) 0],EightBinColor(8,:),'HandleVisibility','off','FaceAlpha',1,'AlphaDataMapping','none');
pat.LineStyle = 'none';
pat.FaceVertexAlphaData = 0.1;
pat.FaceAlpha = 'flat';
uistack(pat,'bottom');
axis([-100 100 -inf inf])
title('Healthy RBC Oscillation Cohort')
hold off


%% Save Figures
%Combined Healthy
saveas(CombVentHistFig,[DataFolderLocation,'\CombVent',CohortTag,'.tif']);
saveas(CombDissHistFig,[DataFolderLocation,'\CombDiss',CohortTag,'.tif']);
saveas(CombBarsHistFig,[DataFolderLocation,'\CombBar',CohortTag,'.tif']);
saveas(CombRBCsHistFig,[DataFolderLocation,'\CombRBC',CohortTag,'.tif']);
saveas(CombRBCBarHistFig,[DataFolderLocation,'\CombRBC_Bar',CohortTag,'.tif']);
saveas(CombRBCOscHistFig,[DataFolderLocation,'\CombRBC_Osc',CohortTag,'.tif']);
%Individual Healthy
saveas(IndVentHistFig,[DataFolderLocation,'\IndVent',CohortTag,'.tif']);
saveas(IndDissHistFig,[DataFolderLocation,'\IndDiss',CohortTag,'.tif']);
saveas(IndBarsHistFig,[DataFolderLocation,'\IndBar',CohortTag,'.tif']);
saveas(IndRBCsHistFig,[DataFolderLocation,'\IndRBC',CohortTag,'.tif']);
saveas(IndRBC_BarsHistFig,[DataFolderLocation,'\IndRBC_Bar',CohortTag,'.tif']);
saveas(IndRBCOscHistFig,[DataFolderLocation,'\IndRBC_Osc',CohortTag,'.tif']);


%% Save Data
close all;%no need to save figures
clearvars HealthyGasExchangeStruct HealthyGasExchangeCell%clear large room for quicker saving/loading
save([DataFolderLocation,'\HealthyCohort']);


%% Display Results
disp(['A total of ', num2str(NumHealthyData), ' healthy subject data were found.'])
disp([num2str(HealthyErrors), ' of the ', num2str(HealthyCohortSize), ' subjects did not have GasExchangeWorkspace.'])


%% Finished
disp('Healthy Cohort Processing Completed!')