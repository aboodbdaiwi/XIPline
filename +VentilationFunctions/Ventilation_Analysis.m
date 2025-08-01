
function [Ventilation] = Ventilation_Analysis(Ventilation,Proton,MainInput)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% Matlab script for processing VDP Analysis  %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Inputs:
%      
%   Outputs:
%                   
%   Package: https://github.com/aboodbdaiwi/xipline
%
%   Author: Abdullah S. Bdaiwi
%   Work email: abdullah.bdaiwi@cchmc.org
%   Personal email: abdaiwi89@gmail.com
%   Website: https://www.cincinnatichildrens.org/research/divisions/c/cpir
%
%   Please add updates at the end. Ex: 3/10/24 - ASB: update ....

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1) Apply the settings:
% We shall use a custom 'input_params' function to define the initial
% settings:
close all;
% f = waitbar(0,'Processing ventilation analysis...');
% pause(.1) 
Ventilation.Image = double(Ventilation.UncorrectedImage);
Ventilation.UncorrectedImage = Ventilation.Image;
Ventilation.LungMask = double(Ventilation.LungMask);
Ventilation.AirwayMask = double(Ventilation.AirwayMask);

settings = VentilationFunctions.input_params(...
1,... % Median filter
0,... % RF correction
1,... % Save data
1,... % Calculate SNR
Ventilation.N4Analysis,... % N4 analysis
Ventilation.IncompleteThresh,... % Incomplete threshold
Ventilation.CompleteThresh,... % Complete threshold
Ventilation.HyperventilatedThresh); % Hyperventilated threshold
Ventilation.MedianFilter = 'yes';
% Usage: settings = Functions.input_params(medfilter, RFcorrection,
% savedata, calculateSNR, N4, incomplete, complete, hyper)

MR = Ventilation.Image;
mkdir([MainInput.OutputPath '\Ventilation_Analysis']);
parentPath = [MainInput.OutputPath '\Ventilation_Analysis\'];
Ventilation.parentPath = parentPath;
%parentPath = [MainInput.XeDataLocation,'\'];
FileNames = Ventilation.filename;

try
    Ventilation.AirwayMask = double(Ventilation.AirwayMask);
    airwaymask = Ventilation.AirwayMask;
catch
    Ventilation.AirwayMask = zeros(size(Ventilation.LungMask));
    airwaymask = Ventilation.AirwayMask;
end

if ~isfield(Ventilation, 'VesselMask')
    Ventilation.VesselMask = zeros(size(Ventilation.LungMask));
else
    Ventilation.VesselMask = double(Ventilation.VesselMask);
end

combined_mask = zeros(size(Ventilation.LungMask), 'double');

% Assign labels with priority: vessels > airway > lung
combined_mask(Ventilation.LungMask > 0) = 1;
combined_mask(Ventilation.AirwayMask > 0) = 2;
combined_mask(Ventilation.VesselMask > 0) = 3;
Ventilation.combined_mask = combined_mask;

maskarray = double(Ventilation.LungMask);
maskarray(Ventilation.AirwayMask == 1) = 0;
maskarray(Ventilation.VesselMask == 1) = 0;

% maskarray = double(Ventilation.LungMask - Ventilation.VesselMask);
% maskarray(maskarray > 1) = 0;
% maskarray(maskarray < 1) = 0;

maskarray = double(maskarray);
% imslice(maskarray)
%% Healthy Reference
Ventilation.HealthyRef.Date = '20250627'; 
Ventilation.HealthyRef.CoV = '0.13±0.02';
Ventilation.HealthyRef.VHI = '7.7±1.8';

Age = MainInput.Age;
% TH method
Ventilation.HealthyRef.TH.VDPULN = ['≤',num2str(round((0.49633 + 0.060413* Age + 1.6449 * 1.9806),2))];
Ventilation.HealthyRef.TH.LVV = '2.3±1.8';
Ventilation.HealthyRef.TH.HVV = '0.05±0.1';

% LB method
if isfield(Ventilation, 'LB_Normalization')
    switch Ventilation.LB_Normalization
        case 'LBmean'
            Ventilation.HealthyRef.LB.VDPULN = ['≤', num2str(round((-0.31238 + 0.020369 * Age + 1.6449 * 0.68247),2))];
            Ventilation.HealthyRef.LB.LVV = '6±2.3';
            Ventilation.HealthyRef.LB.HVV = '5.5±3.3';        
        case 'GLBmean'
            Ventilation.HealthyRef.LB.VDPULN = ['≤', num2str(round((0.11008 + 0.051214 * Age + 1.6449 * 1.5833),2))];
            Ventilation.HealthyRef.LB.LVV = '10.6±3';
            Ventilation.HealthyRef.LB.HVV = '10±3.5';  
        case 'GLBpercentile'
            Ventilation.HealthyRef.LB.VDPULN = ['≤', num2str(round((-0.65729 + 0.06534 * Age + 1.6449 * 1.874),2))];
            Ventilation.HealthyRef.LB.LVV = '11±8';
            Ventilation.HealthyRef.LB.HVV = '14±8.7';  
        otherwise
            warning('Unrecognized LB_Normalization method: %s', Ventilation.LB_Normalization);
    end
else
    warning('Field "LB_Normalization" not found in Ventilation structure.');
end

% Kmeans method
Ventilation.HealthyRef.HK.VDPULN = ['≤',num2str(round((-0.15823 + 0.015696* Age + 1.6449 * 0.59544),2))];
Ventilation.HealthyRef.HK.LVV = '0.8±3.4';
Ventilation.HealthyRef.HK.HVV = '15±5';
% AKmeans method
Ventilation.HealthyRef.AK.VDPULN = ['≤',num2str(round((-0.31238 + 0.020369* Age + 1.6449 * 0.68247),2))];
Ventilation.HealthyRef.AK.LVV = '30±11';
Ventilation.HealthyRef.AK.HVV = '19±4.8';

Ventilation.HealthyRef.skewness = '0.04±0.28';
Ventilation.HealthyRef.kurtosis = '4.0±1.4';
Ventilation.HealthyRef.DDI2D = '≤1';
Ventilation.HealthyRef.DDI3D = '≤3';
Ventilation.HealthyRef.AgeCorrected = 'yes';

Ventilation.writereport = 'yes';
% MainInput.ImageQuality = '';
% MainInput.Note = '';
%% Calculate SNR:
switch settings.calculate_SNR
    case "yes"
%         waitbar(.20,f,'Calculating SNR...');
%         pause(.1)        
        disp('Calculating SNR...')
        [Ventilation] = VentilationFunctions.calculate_SNR(Ventilation);       
    case "no"
        disp('SNR calculation has been skipped.')
    otherwise
        disp('Input error. Check inputs/outputs of input_params().')
end
close all;

%% 5) Apply RF correction to images:
% Generally, we choose not to perform RF correction on images if we run N4
% bias correction later on.
switch settings.RF_correction
    case "yes"
        disp('RF correction performed on images.')
        MR = VentilationFunctions.rfCorrection(MR,maskarray);
    case "no"
        disp('RF correction has been skipped.')
    otherwise
        disp('Input error. Check inputs/outputs of input_params().')
end

%% 6) Apply N4 Bias Correction:
% NOTE: There are opportunities to make this function more robust, such as
% making sure it is functional for 2D images as well as 3D (it currently is
% optimized for 3D imagesets). Similarly, we could also provide it with
% additional settings inputs, but I chose to keep it simple and only
% provide the optimal settings as chosen when I was the research assistant.
% -- written by Joey Plummer, 05/10/2021

switch settings.N4_bias_analysis
    case "yes"    
        %waitbar(.10,f,'Performing N4 bias analysis...');
%         pause(.1)
        disp('Running N4 bias field correction...')
        [N4, Bias] = VentilationFunctions.N4_bias_correction(MR, maskarray, parentPath);
        Ventilation.Bias = Bias;       
        % Output Nifti file of N4 corrected images to save for later use.
        if isempty(FileNames)
            FileNames = 'Ventilation_Image.nii.gz';
        end
        [~,Name,ext] = fileparts(FileNames);
        if iscell(ext) == 0 % Deal with multiple dicom vs single dicom/nifti upload.
            % skip
        else
            ext = ext{1};
            Name = Name{1};
        end
        
        if ext == ".gz"
            [~,Name,~] = fileparts(Name); % Deal with double extension (.nii.gz)
        else
        end
        
        % Apply rotations and flip to account for niftiwrite rotation.
        N4_2 = imrotate(N4,90); 
        N4_2 = flip(N4_2,1);
        
        % Save images as Nifti.
        % NameN4 = Name + "N4";
%         niftiwrite(abs(N4_2),[parentPath,char(NameN4)]);
        % niftiwrite(abs(N4_2),parentPath + "Ventilation_ImagesN4"); % Or do this until we have a naming convention        
        MR = double(N4);
        Ventilation.Image_uncorrected = Ventilation.Image;
        Ventilation.Image = MR;
    case "no"
        disp('N4 bias field correction has been skipped.')
    otherwise
        disp('Input error. Check inputs/outputs of input_params().')
end
%% %% save image and mask
% Generate timestamp suffix
timestamp = lower(datestr(now, '_yyyymmddHHMM')); % e.g., '_202507161430'
Ventilation.timestamp = timestamp(2:end);         % removes the leading '_'
% ==== File Deletion Helper Function ====
delete_if_exist = @(pattern) cellfun(@(f) delete(fullfile(parentPath, f)), ...
    {dir(fullfile(parentPath, pattern)).name}, 'UniformOutput', false);

% ==== Delete any older matching files ====
delete_if_exist('uncorrectedimage*.nii*');
delete_if_exist('image*.nii*');
delete_if_exist('combined_mask*.nii*');
delete_if_exist('lungmask*.nii*');
delete_if_exist('airwaymask*.nii*');
delete_if_exist('vesselsmask*.nii*');

% Image
img_name = lower(['uncorrectedimage' timestamp '.nii']);
niftiwrite(abs(fliplr(rot90(Ventilation.UncorrectedImage, -1))), fullfile(parentPath, img_name), 'Compressed', true);
info = niftiinfo(fullfile(parentPath, [img_name '.gz']));
info.Description = 'Package Version: Version1';
niftiwrite(abs(fliplr(rot90(Ventilation.UncorrectedImage, -1))), fullfile(parentPath, img_name), info, 'Compressed', true);


% Image
img_name = lower(['image' timestamp '.nii']);
niftiwrite(abs(fliplr(rot90(Ventilation.Image, -1))), fullfile(parentPath, img_name), 'Compressed', true);
info = niftiinfo(fullfile(parentPath, [img_name '.gz']));
info.Description = 'Package Version: Version1';
niftiwrite(abs(fliplr(rot90(Ventilation.Image, -1))), fullfile(parentPath, img_name), info, 'Compressed', true);

% combined_mask
combined_mask_name = lower(['combined_mask' timestamp '.nii']);
niftiwrite(abs(fliplr(rot90(Ventilation.combined_mask, -1))), fullfile(parentPath, combined_mask_name), 'Compressed', true);
info = niftiinfo(fullfile(parentPath, [combined_mask_name '.gz']));
info.Description = 'Package Version: Version1';
niftiwrite(abs(fliplr(rot90(Ventilation.combined_mask, -1))), fullfile(parentPath, combined_mask_name), info, 'Compressed', true);

% Lung mask
lung_name = lower(['lungmask' timestamp '.nii']);
niftiwrite(abs(fliplr(rot90(Ventilation.LungMask, -1))), fullfile(parentPath, lung_name), 'Compressed', true);
info = niftiinfo(fullfile(parentPath, [lung_name '.gz']));
info.Description = 'Package Version: Version1';
niftiwrite(abs(fliplr(rot90(Ventilation.LungMask, -1))), fullfile(parentPath, lung_name), info, 'Compressed', true);

% Airway mask
airway_name = lower(['airwaymask' timestamp '.nii']);
niftiwrite(abs(fliplr(rot90(Ventilation.AirwayMask, -1))), fullfile(parentPath, airway_name), 'Compressed', true);
info = niftiinfo(fullfile(parentPath, [airway_name '.gz']));
info.Description = 'Package Version: Version1';
niftiwrite(abs(fliplr(rot90(Ventilation.AirwayMask, -1))), fullfile(parentPath, airway_name), info, 'Compressed', true);

% Vessels mask
vessels_name = lower(['vesselsmask' timestamp '.nii']);
niftiwrite(abs(fliplr(rot90(Ventilation.VesselMask, -1))), fullfile(parentPath, vessels_name), 'Compressed', true);
info = niftiinfo(fullfile(parentPath, [vessels_name '.gz']));
info.Description = 'Package Version: Version1';
niftiwrite(abs(fliplr(rot90(Ventilation.VesselMask, -1))), fullfile(parentPath, vessels_name), info, 'Compressed', true);


%% calculate ventilation heterogeneity
if strcmp(Ventilation.HeterogeneityIndex,'yes') 
    [Ventilation] = VentilationFunctions.calculate_VHI(Ventilation);
else
    Ventilation.CV_maps = [];
    Ventilation.sliceMeanCV = [];
    Ventilation.sliceVHI = [];
    Ventilation.overallMeanCV = []; 
    Ventilation.overallVHI = [];
end
%% Proton image
if MainInput.NoProtonImage == 1 
    Proton.Image = zeros(size(MR));
    Proton.Image(end,end,end) = 1; % to avoid error for max 
    Proton.ProtonRegistered = Proton.Image;
    Proton.ProtonRegisteredColored = Proton.Image;
end
%% Skewness
lung_voxels = Ventilation.Image(maskarray == 1);
Ventilation.skewness = skewness(double(lung_voxels));
Ventilation.kurtosis = kurtosis(double(lung_voxels)); 
%% 8) Calculate VDP:

if strcmp(Ventilation.ThreshAnalysis,'yes')   % 'yes'; || 'no'

    [Ventilation] = VentilationFunctions.calculate_VDP_CCHMC(Ventilation,Proton,MainInput);    

    disp(['The overall ventilation defect percentage is: ',num2str(Ventilation.Threshold.VDP),'%'])
end 
close all;
%%  Perform linear binning ventilation analysis:
if strcmp(Ventilation.LB_Analysis,'yes') == 1 
    maskarraytrachea = maskarray;
    maskarraytrachea(airwaymask == 1)=0;    
    Ventilation.maskarraytrachea = maskarraytrachea;
    [Ventilation] = VentilationFunctions.calculate_LB_VDP(Ventilation,Proton,MainInput);    
end
close all;
%% K-means VDP
if strcmp(Ventilation.Kmeans,'yes')
    Ventilation.parentPath = parentPath;
    % Riaz code
    % [Ventilation] = VentilationFunctions.calculate_kMeans_VDP(Ventilation,Proton,MainInput); 
    % Samal's code
    % [Ventilation] = VentilationFunctions.kmeansCalc(Ventilation,Proton,MainInput);
    % Grace's (Rachel) code
    [Ventilation] = VentilationFunctions.kmeans.VDP_calculationASB(Ventilation,Proton,MainInput);
end
close all;
%% Adaptive K-means VDP
if strcmp(Ventilation.AKmeans,'yes')
    Ventilation.parentPath = parentPath;
    [Ventilation] = VentilationFunctions.Akmeans.VDP_calculationASB(Ventilation,Proton,MainInput);
end
close all;
%% DDI
if strcmp(Ventilation.DDI2D,'yes') || strcmp(Ventilation.DDI3D,'yes')
    switch Ventilation.DDIDefectMap
        case 'Threshold'
            defect_mask = Ventilation.Threshold.defectArray;
            defect_mask(defect_mask > 2) = 0;
            defect_mask = double(defect_mask > 0);
            defectMap = double(maskarray) + defect_mask;  
        case 'Linear Binning'
            defect_mask = Ventilation.VentBinMap2;
            defect_mask(defect_mask > 1) = 0;
            defectMap = double(maskarray) + defect_mask; 
        case 'Kmeans'
            defect_mask = Ventilation.Kmeans_segmentation;
            defect_mask(defect_mask > 0 & defect_mask < 2) = 1;
            defect_mask(defect_mask > 1) = 0;
            defectMap = double(maskarray) + defect_mask; 
        case 'AKmeans'
            defect_mask = Ventilation.Akmeans_defect_mask;
            defect_mask(defect_mask > 0 & defect_mask < 2) = 1;
            defect_mask(defect_mask > 1) = 0;
            defectMap = double(maskarray) + defect_mask;             
    end
    Ventilation.defectMap_forDDI = defectMap;
end

if strcmp(Ventilation.DDI2D,'yes')
    [Ventilation] = VentilationFunctions.calculateDDI_2D(Ventilation,Proton,MainInput);
end
if strcmp(Ventilation.DDI3D,'yes')
    [Ventilation] = VentilationFunctions.calculateDDI_3D(Ventilation,Proton,MainInput);
end
close all;
%% GLRLM analysis
if Ventilation.GLRLM_Analysis == "yes" % 'yes'; || 'no'
%     f = waitbar(.75,'Performing Linear Binning Analysis ...');
%     waitbar(.90,f,'Performing Texture Analysis ...');
%     pause(.1)     
    cd(parentPath)    
    switch Ventilation.GLRLMDefectMap
        case 'Threshold'
            defect_mask = Ventilation.Threshold.defectArray;
            defect_mask(defect_mask > 2) = 0;
            defect_mask = double(defect_mask > 0);
            defectMap = double(maskarray) + defect_mask;  
        case 'Linear Binning'
            defect_mask = Ventilation.VentBinMap2;
            defect_mask(defect_mask > 1) = 0;
            defectMap = double(maskarray) + defect_mask; 
        case 'Kmeans'
            defect_mask = Ventilation.Kmeans_segmentation;
            defect_mask(defect_mask > 0 & defect_mask < 2) = 1;
            defect_mask(defect_mask > 1) = 0;
            defectMap = double(maskarray) + defect_mask; 
        case 'AKmeans'
            defect_mask = Ventilation.Akmeans_defect_mask;
            defect_mask(defect_mask > 0 & defect_mask < 2) = 1;
            defect_mask(defect_mask > 1) = 0;
            defectMap = double(maskarray) + defect_mask;             
    end
    Ventilation.defectMap_forGLRLM = defectMap;    
    [Ventilation] = VentilationFunctions.GLRLM_Analysis(Ventilation);
end 
close all;
%% write report
if strcmp(Ventilation.writereport,'yes')
    if strcmp(Ventilation.ThreshAnalysis,'yes')
        files = dir(fullfile(Ventilation.outputpath, 'ThresholdVDP_Report*.pptx')); 
        for i = 1:length(files)
            delete(fullfile(Ventilation.outputpath, files(i).name));
        end
        VentilationFunctions.ThresholdVDP_Report(Ventilation,Proton, MainInput);
    end
    if strcmp(Ventilation.LB_Analysis,'yes')
        files = dir(fullfile(Ventilation.outputpath, 'LBVDP_Report*.pptx'));
        for i = 1:length(files)
            delete(fullfile(Ventilation.outputpath, files(i).name));
        end        
        VentilationFunctions.LBVDP_Report(Ventilation,Proton, MainInput);
    end
    if strcmp(Ventilation.Kmeans,'yes')
        files = dir(fullfile(Ventilation.outputpath, 'KmeansVDP_Report*.pptx'));
        for i = 1:length(files)
            delete(fullfile(Ventilation.outputpath, files(i).name));
        end          
        VentilationFunctions.KmeansVDP_Report(Ventilation,Proton, MainInput);
    end
    if strcmp(Ventilation.AKmeans,'yes')
        files = dir(fullfile(Ventilation.outputpath, 'AKmeansVDP_Report*.pptx'));
        for i = 1:length(files)
            delete(fullfile(Ventilation.outputpath, files(i).name));
        end          
        VentilationFunctions.AKmeansVDP_Report(Ventilation,Proton, MainInput);
    end    
end

%% %% save maps in mat file
save_data=[parentPath,'\','Ventilation_Analysis','.mat'];
save(save_data);  

VentilationExcelFile = fullfile(parentPath, 'Ventilation_workspace.xlsx');
Global.exportStructToExcel(Ventilation, VentilationExcelFile);
MainInputExcelFile = fullfile(parentPath, 'MainInput_workspace.xlsx');
Global.exportStructToExcel(MainInput, MainInputExcelFile);
ProtonExcelFile = fullfile(parentPath, 'Proton_workspace.xlsx');
Global.exportStructToExcel(Proton, ProtonExcelFile);

disp('ventilation analysis completed');
close all;
end %end of function



% Please add updates here
% 
% 
% 
% 

