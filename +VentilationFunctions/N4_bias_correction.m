function [N4, Bias] = N4_bias_correction(MR, maskarray, parentPath)
%% MATLAB script to perform N4 bias correction
% This code is build from the original Python code used by our department
% for N4 bias correction of grayscale imagesets. The function will call in
% a .exe file to run. There are many features to play around with in this
% function. 
%
% Instructions:
% Run: [N4, Bias] = Functions.N4_bias_correction(MR, maskarray, parentPath);
%
% Inputs:
% MR = ventilation image data
% maskarray = corresponding masks for MR 
% parentPath = pathname containing the MR images
%
% Outputs:
% N4 = N4 bias corrected images according to code last updated in August
% 2019
% Bias = bias field correction map
%
% Remarks:
% Feel free to play around with the code, and adjust any settings anywhere.
% Default settings are chosen as optimized during the 2018-2020 period. As
% the field progresses, we may want to change these parameters.
%
% Author: Joseph Plummer.
% Date: 05/05/2021.
%% Potential Inputs
Image = MR;
% Mask = ones(size(Image)); % ones mask
Mask = imbinarize(maskarray); % simple binary mask
% Mask = imdilate(Mask, strel('sphere',2)); % dilated binary mask
% Weight = ones(size(Image)); % ones weights
Weight = double(Image); % weight on image intensity
% Weight = double(Image).^2; % weight on intensity cubed
% Weight(Image>prctile(Image(:),99.9)) = 0; %Not indcluding airway

% %% View Inputs
% figure('Name','Input Image')
% montage(Image, 'DisplayRange', [0 0.5*max(Image(:))])
% figure('Name','Mask')
% montage(Mask, 'DisplayRange', [0 max(Mask(:))])
% figure('Name','Weights')
% montage(Weight, 'DisplayRange', [0 max(Weight(:))])


%% Set location of exe and temp images
% Location of N4BiasFieldCorrection.exe
% N4bias_location = 'C:\Work_LAB\Xe_App\Main App functions\+VentilationFunctions\';
%% Set location of exe and temp images
%location of N4BiasFieldCorrection
N4Path = mfilename('fullpath');
idcs = strfind(N4Path,filesep);%determine location of file separators
N4Path = [N4Path(1:idcs(end)-1),filesep];%remove file
%% Export Image and Mask temporarily for use in N4
niftiwrite(abs(Image),[parentPath,'Image.nii']); % Image
niftiwrite(abs(Mask),[parentPath,'Mask.nii']); % Mask
niftiwrite(abs(Weight),[parentPath,'Weight.nii']); % Weight

%% Run Bias Correction
% Default settings:
% cmd = ['"',parentPath,'N4BiasFieldCorrection.exe" ',...% run bias correction
%     '-d 3 ',... % set to 3 dimensions
%     '-i "',parentPath,'Image.nii" ',... % input image called of Image.nii
%     '-x "',parentPath,'Mask.nii" ',... % import mask called Mask.nii
%     '-w "',parentPath,'Weight.nii" ',... % import mask called Weight.nii
%     '-s 1 ',... % shrink factor
%     '-c [50,0.001] ',... % convergence
%     '-b [5,3] ',... % spline settings
%     '-t [0.15,0.01,200] ',... % histogram settings
%     '-o ["',parentPath,'CorrectedImage.nii","',parentPath,'Bias.nii"]']; % output corrected image and bias field

% Chosen settings:
cmd = ['"',N4Path,'N4BiasFieldCorrection.exe" ',...%run bias correction
    '-d 3 ',... % set to 3 dimensions - MUST MATCH DIMENSIONS OF IMAGE
    '-i "',parentPath,'Image.nii" ',... % input image called of Image.nii
    '-x "',parentPath,'Mask.nii" ',... % import mask called Mask.nii
    '-w "',parentPath,'Weight.nii" ',... % import mask called Weight.nii
    '-s 1 ',... % shrink factor
    '-c [50,0.001] ',... % convergence
    '-b [5,3] ',... % spline settings
    '-t [0.15,0.01,200] ',... % histogram settings
    '-o ["',parentPath,'CorrectedImage.nii","',parentPath,'Bias.nii"]']; % output corrected image and bias field
system(cmd);

%% Import Results
N4 = double(niftiread(strcat(parentPath,"CorrectedImage.nii")));
Bias = double(niftiread(strcat(parentPath,"Bias.nii")));

%% Delete Temp Files
delete([parentPath,'Image.nii']);
delete([parentPath,'Mask.nii']);
delete([parentPath,'Weight.nii']);
delete([parentPath,'CorrectedImage.nii']);
delete([parentPath,'Bias.nii']);

end