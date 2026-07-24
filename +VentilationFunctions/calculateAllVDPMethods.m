function [Ventilation, Outputs, VDPTable] = calculateAllVDPMethods( ...
    Ventilation, Proton, MainInput, subjectID)
% calculateAllVDPMethods
% Calculates VDP using all enabled methods and saves one row per subject
% to an Excel summary file.
%
% Methods: TH60, K-means, Adaptive K-means, LBmean, GLBmean, HyLBm, GLB99p.
%
% Example:
% excelFile = fullfile(Ventilation.parentPath,'All_VDP_Methods.xlsx');
% [Ventilation, Outputs, VDPTable] = calculateAllVDPMethods( ...
%     Ventilation, Proton, MainInput, subjectID, excelFile);

if nargin < 4 || isempty(subjectID)
    error('A subjectID must be provided.');
end
subjectID = string(subjectID);

excelFile = fullfile(Ventilation.parentPath,'All_VDP_Methods.xlsx');

excelFile = char(excelFile);
excelFolder = fileparts(excelFile);
if ~isempty(excelFolder) && ~isfolder(excelFolder)
    mkdir(excelFolder);
end

if isfield(Ventilation,'parentPath')
    parentPath = Ventilation.parentPath;
else
    parentPath = '';
end

Outputs = struct();
Outputs.N4.Image = Ventilation.Image;

VDP_TH60 = NaN;
VDP_Kmeans = NaN;
VDP_AKmeans = NaN;
VDP_LBmean = NaN;
VDP_GLBmean = NaN;
VDP_GLB99p = NaN;
VDP_HyLBm = NaN;

%% Exclude vessel voxels from lung mask
if ~isfield(Ventilation,'LungMask')
    error('Ventilation.LungMask is missing.');
end

if isfield(Ventilation,'VesselMask') && ~isempty(Ventilation.VesselMask)
    maskarray = double(Ventilation.LungMask + Ventilation.VesselMask) ...
        .* double(Ventilation.LungMask);
    maskarray(maskarray > 1) = 0;
else
    warning('Ventilation.VesselMask was not found; using the original lung mask.');
    maskarray = double(Ventilation.LungMask);
end
Ventilation.LungMask = double(maskarray);

%% TH60
if isMethodEnabled(Ventilation,'ThreshAnalysis')
    fprintf('Calculating TH60 VDP for %s...\n',subjectID);
    Ventilation.parentPath = parentPath;
    Ventilation = VentilationFunctions.calculate_VDP_CCHMC( ...
        Ventilation,Proton,MainInput);

    Outputs.N4.TH60.VDP = Ventilation.Threshold.VDP;
    Outputs.N4.TH60.defectArray = Ventilation.Threshold.defectArray;
    Outputs.N4.TH60.BinPrecentValues = Ventilation.Threshold.THBins;
    VDP_TH60 = extractScalarVDP(Ventilation.Threshold.VDP);
end
close all force;

%% K-means
if isMethodEnabled(Ventilation,'Kmeans')
    fprintf('Calculating K-means VDP for %s...\n',subjectID);
    Ventilation.parentPath = parentPath;
    Ventilation = VentilationFunctions.kmeans.VDP_calculationASB( ...
        Ventilation,Proton,MainInput);

    Outputs.N4.Kmeans.VDP = Ventilation.KmeansVDP;
    Outputs.N4.Kmeans.defectArray = Ventilation.Kmeans_segmentation;
    Outputs.N4.Kmeans.BinPrecent = Ventilation.wholelung_VDP;
    VDP_Kmeans = extractScalarVDP(Ventilation.KmeansVDP);
end
close all force;

%% Adaptive K-means
if isMethodEnabled(Ventilation,'AKmeans')
    fprintf('Calculating Adaptive K-means VDP for %s...\n',subjectID);
    Ventilation.parentPath = parentPath;
    Ventilation = VentilationFunctions.Akmeans.VDP_calculationASB( ...
        Ventilation,Proton,MainInput);

    Outputs.N4.AKmeans.VDP = Ventilation.Akmeans_VDP;
    Outputs.N4.AKmeans.defectArray = Ventilation.Aksegmentation;
    Outputs.N4.AKmeans.BinPrecent = Ventilation.Akmeans_ventP_raw_cor;
    VDP_AKmeans = extractScalarVDP(Ventilation.Akmeans_VDP);
end
close all force;

%% Linear-binning methods
if isMethodEnabled(Ventilation,'LB_Analysis')
    [Ventilation,Outputs.N4.LBmean,VDP_LBmean] = runLBMethod( ...
        Ventilation,Proton,MainInput,parentPath,'LBmean');
    close all force;

    [Ventilation,Outputs.N4.GLBmean,VDP_GLBmean] = runLBMethod( ...
        Ventilation,Proton,MainInput,parentPath,'GLBmean');
    close all force;

    [Ventilation,Outputs.N4.HybridGLBm,VDP_HyLBm] = runLBMethod( ...
        Ventilation,Proton,MainInput,parentPath,'HybridLBm');
    close all force;

    [Ventilation,Outputs.N4.GLB99p,VDP_GLB99p] = runLBMethod( ...
        Ventilation,Proton,MainInput,parentPath,'GLBpercentile');
    close all force;
end

Ventilation.parentPath = parentPath;

%% One-row summary table
VDPTable = table( ...
    subjectID,VDP_TH60,VDP_Kmeans,VDP_AKmeans,VDP_LBmean, ...
    VDP_GLBmean,VDP_GLB99p,VDP_HyLBm, ...
    'VariableNames',{'SubjectID','TH60_VDP','Kmeans_VDP', ...
    'AKmeans_VDP','LBmean_VDP','GLBmean_VDP','GLB99p_VDP','HyLBm_VDP'});

writeVDPTableToExcel(VDPTable,excelFile);

fprintf('\nVDP calculations completed for %s.\n',subjectID);
fprintf('Excel summary saved to:\n%s\n',excelFile);
end

function [Ventilation,methodOutput,VDP] = runLBMethod( ...
    Ventilation,Proton,MainInput,parentPath,normalizationMethod)

    fprintf('Calculating %s VDP...\n',normalizationMethod);
    Ventilation.parentPath = parentPath;
    Ventilation = configureLinearBinning(Ventilation,normalizationMethod);
    Ventilation = VentilationFunctions.calculate_LB_VDP( ...
        Ventilation,Proton,MainInput);
    
    methodOutput.VDP = Ventilation.LB_VDP;
    methodOutput.defectArray = Ventilation.VentBinMap2;
    methodOutput.BinPrecent = Ventilation.BinsPercent;
    VDP = extractScalarVDP(Ventilation.LB_VDP);
end

function enabled = isMethodEnabled(S,fieldName)
    enabled = false;
    if ~isfield(S,fieldName)
        return;
    end
    value = S.(fieldName);
    if isstring(value) || ischar(value)
        enabled = strcmpi(strtrim(string(value)),"yes");
    elseif islogical(value) || isnumeric(value)
        enabled = logical(value);
    end
end

function value = extractScalarVDP(inputValue)
    if isempty(inputValue)
        value = NaN;
        return;
    end
    inputValue = double(inputValue);
    if ~isscalar(inputValue)
        warning('VDP output had multiple values; the first value was used.');
        inputValue = inputValue(1);
    end
    value = inputValue;
end

function Ventilation = configureLinearBinning(Ventilation,method)
    Ventilation.LB_Normalization = method;
    
    switch method
        case 'LBmean'
            Ventilation.LBThresholds = [0.33,0.66,1.00,1.33,1.66];
            Ventilation.Hdist = [0.001029,0.995768,0.244090];
    
        case 'GLBmean'
            Ventilation.LBThresholds = ...
                [0.495130,0.752533,1.001607,1.244833,1.483568];
            Ventilation.Hdist = [-0.771445,1.136303,0.282776];
    
        case 'HybridLBm'
            Ventilation.LBThresholds = [0.50,0.75,1.00,1.25,1.50];
            Ventilation.Hdist = [-0.771445,1.136303,0.282776];
    
        case 'GLBmed'
            Ventilation.LBThresholds = ...
                [0.489025,0.744938,0.991227,1.230820,1.465304];
            Ventilation.Hdist = [-0.824047,1.131236,0.282992];
    
        case 'GLBpercentile'
            Ventilation.LBThresholds = ...
                [0.288930,0.462393,0.622368,0.773740,0.918873];
            Ventilation.Hdist = [-1.186698,0.738826,0.198451];
    
        otherwise
            error('Unknown linear-binning method: %s',method);
    end
end

function writeVDPTableToExcel(newRow,excelFile)
    expectedVariables = newRow.Properties.VariableNames;
    
    if isfile(excelFile)
        try
            existingTable = readtable(excelFile,'Sheet','VDP_Summary', ...
                'TextType','string');
        catch readError
            warning(['Could not read the existing Excel summary: %s\n', ...
                'A new VDP summary sheet will be created.'],readError.message);
            existingTable = table();
        end
    else
        existingTable = table();
    end
    
    if isempty(existingTable)
        combinedTable = newRow;
    else
        missingVariables = setdiff(expectedVariables, ...
            existingTable.Properties.VariableNames);
        if ~isempty(missingVariables)
            error(['The existing VDP_Summary sheet is missing: %s. ', ...
                'Use a new Excel filename or correct the sheet.'], ...
                strjoin(missingVariables,', '));
        end
    
        existingTable = existingTable(:,expectedVariables);
        existingTable.SubjectID = string(existingTable.SubjectID);
        newRow.SubjectID = string(newRow.SubjectID);
    
        matchingRows = existingTable.SubjectID == newRow.SubjectID;
        existingTable(matchingRows,:) = [];
        combinedTable = [existingTable;newRow];
    end
    
    combinedTable = sortrows(combinedTable,'SubjectID');
    writetable(combinedTable,excelFile,'Sheet','VDP_Summary', ...
        'WriteMode','overwritesheet');
end
