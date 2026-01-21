
function [Diffusion] = Diffusion_Analysis(Diffusion,MainInput)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%% Diffusion Analysis  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%   Inputs:
%      
%   Outputs:
%                   
%   Package: https://github.com/aboodbdaiwi/XIPline
%
%   Author: Abdullah S. Bdaiwi
%   Work email: abdullah.bdaiwi@cchmc.org
%   Personal email: abdaiwi89@gmail.com
%   Website: https://www.cincinnatichildrens.org/research/divisions/c/cpir
%
%   Please add updates at the end. Ex: 3/10/24 - ASB: updated flip angle caluclation for spiral diffusion
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% extract varabiles
% f = waitbar(0,'Please wait...');
% pause(.5)

WinBUGSPath = mfilename('fullpath');
idcs = strfind(WinBUGSPath,filesep);%determine location of file separators
WinBUGSPath = [WinBUGSPath(1:idcs(end)-1),filesep];%remove file
WinBUGSPath = [WinBUGSPath,'\WinBUGS14'];

% WinBUGSPath = 'P:\OneDrive - cchmc\Lab\WinBUGS14';

Diffusion.WinBUGSPath = WinBUGSPath;  
Diffusion.MA_WinBUGSPath = WinBUGSPath;

try
    PatientAge = eval(MainInput.PatientAge);
catch
    PatientAge = MainInput.Age;
end

% ADC_Analysis = Diffusion.ADC_Analysis;
ADC_Analysis = 'yes';
Diffusion.ADCLB_RefMean = 0.0002*MainInput.Age+0.029; % mean equetion for healthy co. 
Diffusion.ADCLB_RefSD = 5e-5*MainInput.Age+0.0121; 


ADCFittingType = Diffusion.ADCFittingType;

try
    bvalues = eval(Diffusion.bvalues);
catch
    bvalues = Diffusion.b_values;
end
Diffusion.b_values = bvalues;

ADCLB_Analysis = Diffusion.ADCLB_Analysis;
MorphometryAnalysis = Diffusion.MorphometryAnalysis;
% ADCAnalysisType = Diffusion.ADCAnalysisType;  % human | animals
ADCAnalysisType = 'human';

diffimg = Diffusion.Image;
lung_mask = Diffusion.LungMask;
try
    airway_mask = Diffusion.AirwayMask;
catch
    airway_mask = zeros(size(lung_mask));
end
lung_mask = (lung_mask - airway_mask).*lung_mask;
% Diffusion.LungMask = lung_mask;
% figure; imslice(airway_mask)

% create a final mask
final_mask = lung_mask;
final_mask(airway_mask==1)=0;

% create the noise mask
SE = strel('square', 10);
maskarray_dilated = imdilate(final_mask, SE);
airwaymask = Segmentation.SegmentLungthresh(diffimg(:,:,:,1),1,0.3);
zerofillings = double(diffimg(:,:,:,1) == 0);
noise_mask = double(imcomplement(maskarray_dilated).*~airwaymask.*~zerofillings);


%% get rid of any slice that has no signal (un-segmented)
slices = size(diffimg,3);
Nbvalues = length(bvalues);
slice_checker=zeros(1,slices);
for check_slice=1:slices
    sl=lung_mask(:,:,check_slice);
    if sum(sl(:))>1
        slice_checker(check_slice)=1;
    elseif sum(sl(:))== 0
        slice_checker(check_slice) = 0;
    end
end
find_slice_checker = find(slice_checker);
Ndiffimg = zeros(size(lung_mask,1),size(lung_mask,2),length(find_slice_checker), Nbvalues);
Nlung_mask = zeros(size(lung_mask,1),size(lung_mask,2),length(find_slice_checker));
Nairway_mask = zeros(size(lung_mask,1),size(lung_mask,2),length(find_slice_checker));
Nfinal_mask = zeros(size(lung_mask,1),size(lung_mask,2),length(find_slice_checker));
Nnoise_mask = zeros(size(lung_mask,1),size(lung_mask,2),length(find_slice_checker));
for find_sl=1:length(find_slice_checker)
    choose_sl=find_slice_checker(find_sl);
    Ndiffimg(:,:,find_sl,:)=diffimg(:,:,choose_sl,:);
    Nlung_mask(:,:,find_sl)=lung_mask(:,:,choose_sl);
    Nairway_mask(:,:,find_sl)=airway_mask(:,:,choose_sl);
    Nfinal_mask(:,:,find_sl)=final_mask(:,:,choose_sl);
    Nnoise_mask(:,:,find_sl)=noise_mask(:,:,choose_sl);
end
% figure; imslice(Ndiffimg)
% figure; imslice(Nfinal_mask)
diffimg = Ndiffimg;
Diffusion.Image = diffimg;
lung_mask = Nlung_mask;
airway_mask = Nairway_mask;
final_mask = Nfinal_mask;
noise_mask = Nnoise_mask;
slices = length(find_slice_checker);

Diffusion.lung_mask = lung_mask;
Diffusion.airway_mask = airway_mask;
Diffusion.final_mask = final_mask;
Diffusion.noise_mask = noise_mask;

%% create path
Nbvalues = length(bvalues);

if ~isempty(MainInput.OutputPath)
    outputpath = MainInput.OutputPath;
else 
    DataLocation = MainInput.XeDataLocation;
    cd(DataLocation)
    mkdir([DataLocation '\Diffusion_Analysis']);
    outputpath = [DataLocation '\Diffusion_Analysis\'];
end
cd(outputpath)
%% Healthy Reference 
Diffusion.HealthyRef.ADC = '0.00±0.00';
Diffusion.HealthyRef.LDR = '0.00±0.00';
Diffusion.HealthyRef.NDR = '0.00±0.00';
Diffusion.HealthyRef.HDR = '0.00±0.00';
Diffusion.HealthyRef.ADCskewness = '0.00±0.00';
Diffusion.HealthyRef.ADCkurtosis = '0.00±0.00';
Diffusion.HealthyRef.ADCFitR2 = '0.00±0.00';
Diffusion.HealthyRef.DDC = '0.00±0.00';
Diffusion.HealthyRef.LmD = '0.00±0.00';
Diffusion.HealthyRef.alpha = '0.00±0.00';
Diffusion.HealthyRef.R = '0.00±0.00';
Diffusion.HealthyRef.r = '0.00±0.00';
Diffusion.HealthyRef.h = '0.00±0.00';
Diffusion.HealthyRef.Lm = '0.00±0.00';
Diffusion.HealthyRef.SVR = '0.00±0.00';
Diffusion.HealthyRef.Na = '0.00±0.00';
Diffusion.HealthyRef.LmD  = '0.00±0.00';
Diffusion.HealthyRef.AgeCorrected = 'no';

if ~isfield(MainInput, 'ImageQuality') || isempty(MainInput.ImageQuality)
    MainInput.ImageQuality = '';
end

if ~isfield(MainInput, 'Note') || isempty(MainInput.Note)
    MainInput.Note = '';
end
%%  initialize these 
fieldsToInit = {
    'meanADC','stdADC', ...
    'LDR','NDR','HDR', ...
    'ADC_skewness','ADC_kurtosis', ...
    'meanRmap', ...
    'DDC_mean','DDC_std', ...
    'LmD_mean','LmD_std', ...
    'alpha_mean','alpha_std', ...
    'R_mean','R_std', ...
    'r_mean','r_std', ...
    'h_mean','h_std', ...
    'Lm_mean','Lm_std', ...
    'SVR_mean','SVR_std', ...
    'Na_mean','Na_std'
};

for i = 1:numel(fieldsToInit)
    f = fieldsToInit{i};
    Diffusion.(f) = NaN;
end

%% Save images and maskes
% ==== File Deletion Helper Function ====
delete_if_exist = @(pattern) cellfun(@(f) delete(fullfile(outputpath, f)), ...
    {dir(fullfile(outputpath, pattern)).name}, 'UniformOutput', false);

% ==== Delete any older matching files ====
delete_if_exist('DiffusionImage*.nii*');
delete_if_exist('LungMask*.nii*');
delete_if_exist('AirwayMask*.nii*');

% diffusion images
niftiwrite(abs(fliplr(rot90(diffimg,-1))),[outputpath,'\DiffusionImage.nii'],'Compressed',true);
info = niftiinfo([outputpath,'\DiffusionImage.nii.gz']);
info.Description = strcat('Package Version: ', 'Version1');
niftiwrite(abs(fliplr(rot90(diffimg,-1))),[outputpath,'\DiffusionImage.nii'],info,'Compressed',true);
% lung mask
niftiwrite(abs(fliplr(rot90(lung_mask,-1))),[outputpath,'\LungMask.nii'],'Compressed',true);
info = niftiinfo([outputpath,'\LungMask.nii.gz']);
info.Description = strcat('Package Version: ', 'Version1');
niftiwrite(abs(fliplr(rot90(lung_mask,-1))),[outputpath,'\LungMask.nii'],info,'Compressed',true);
% airway mask
niftiwrite(abs(fliplr(rot90(airway_mask,-1))),[outputpath,'\AirwayMask.nii'],'Compressed',true);
info = niftiinfo([outputpath,'\AirwayMask.nii.gz']);
info.Description = strcat('Package Version: ', 'Version1');
niftiwrite(abs(fliplr(rot90(airway_mask,-1))),[outputpath,'\AirwayMask.nii'],info,'Compressed',true);


%% ADC Analysis
% waitbar(.10,f,'Calculating ADC');
% pause(1)
Diffusion.diffimg = diffimg;
Diffusion.lung_mask = lung_mask;
Diffusion.airway_mask = airway_mask;
Diffusion.bvalues = bvalues;
Diffusion.ADCFittingType = ADCFittingType;
Diffusion.ADCAnalysisType = ADCAnalysisType;
Diffusion.WinBUGSPath = WinBUGSPath;
Diffusion.outputpath = outputpath;
Diffusion.PatientAge = PatientAge;
if strcmp(ADC_Analysis, 'yes') 
    [Diffusion] = DiffusionFunctions.ADC_Analysis(Diffusion, MainInput);
end

%% ADC linear binning analysis
% waitbar(.20,f,'Processing Linear Binning Analysis');
% pause(1)

if strcmp(ADCLB_Analysis, 'yes') == 1 && strcmp(ADC_Analysis, 'yes') == 1
    [Diffusion] = DiffusionFunctions.LinearBinningADC(Diffusion, MainInput);   
end
%% Morphometry Analysis
if strcmp(MorphometryAnalysis, 'no') 
    Diffusion.MorphometryAnalysis = 'no';
    Diffusion.SEMMorphometry = 'no';
    Diffusion.CMMorphometry = 'no';
elseif strcmp(MorphometryAnalysis, 'yes')  && Nbvalues < 3
    Diffusion.MorphometryAnalysis = 'no';
    Diffusion.SEMMorphometry = 'no';
    Diffusion.CMMorphometry = 'no';
end
if strcmp(MorphometryAnalysis, 'yes') == 1 && Nbvalues > 3
% CM Morphometry analysis
%     waitbar(.30,f,'Processing CM Morphometry Analysis ');
%     pause(1)
    if strcmp(Diffusion.CMMorphometry, 'yes') 
        cd(outputpath)
        [Diffusion]= DiffusionFunctions.CM_MorphometryFit(Diffusion, MainInput); 

%         [R_map,h_map,r_map,Lm_map,SVR_map,Na_map,DT_map,Dan_map,DL_map,Dmean_map,So_map,MSR_map,Residual_map,LungCMMorphometrySummary]=...
%             DiffusionFunctions.CM_Lung_Morphometry(...
%             diffimg,...
%             final_mask,...
%             noise_mask,...
%             Diffusion.MorphometryAnalysisType,...
%             bvalues,...
%             Diffusion.Do ,...
%             Diffusion.Delta,...
%             outputpath,...
%             WinBUGSPath);
    end


    if strcmp(Diffusion.SEMMorphometry, 'yes') == 1
        cd(outputpath)
        [Diffusion]= DiffusionFunctions.SEM_MorphometryFit(Diffusion, MainInput); 
        
%         [DDC_map,alpha_map,SEMSo_map,LmD_map,LungSEMMorphometrySummary]= DiffusionFunctions.SEM_Morphometry...
%             (diffimg,...
%             final_mask,...
%             noise_mask,...
%             bvalues,...
%             Diffusion.Do ,...
%             Diffusion.Delta,...
%             outputpath,...
%             WinBUGSPath);
    end
elseif strcmp(MorphometryAnalysis, 'yes')  && Nbvalues < 3
    disp('Morphometry cannot be performed with only 2 b values')
end

close all;
%% writere port
if ~isfield(Diffusion, 'writereport') || isempty(Diffusion.writereport)
    Diffusion.writereport = 'yes';
end
if strcmp(Diffusion.writereport,'yes')
    files = dir(fullfile(Diffusion.outputpath, 'DiffusionAnalysis_Report*.pptx')); 
    for i = 1:length(files)
        delete(fullfile(Diffusion.outputpath, files(i).name));
    end
    DiffusionFunctions.DiffusionAnalysis_Report(Diffusion, MainInput);
end
%% save maps in mat file
save_data=[outputpath,'\workspace.mat'];
save(save_data);   

DiffusionExcelFile = fullfile(outputpath, 'Diffusion_workspace.xlsx');
Global.exportStructToExcel(Diffusion, DiffusionExcelFile);
MainInputExcelFile = fullfile(outputpath, 'MainInput_workspace.xlsx');
Global.exportStructToExcel(MainInput, MainInputExcelFile);

% waitbar(1,f,'Finishing');
% pause(1)
% close(f)
end


% Please add updates here
% 
% 
% 
% 

