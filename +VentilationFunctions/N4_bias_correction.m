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
% Author: Joseph Plummer and Abdullah Bdaiwi
% Date: 05/05/2021.

%% Potential Inputs
% Set location of exe and temp images
%location of N4BiasFieldCorrection
N4Path = mfilename('fullpath');
idcs = strfind(N4Path,filesep);%determine location of file separators
N4Path = [N4Path(1:idcs(end)-1),filesep];%remove file
% maskarray = flip(maskarray,3);
if ~exist('maskarray','var')
    maskarray = Segmentation.SegmentLungthresh(MR,1,1); 
else
    Ostumask = Segmentation.SegmentLungthresh(MR,1,1); 
    maskarray = maskarray + Ostumask;
    maskarray = maskarray > 0;
end

% Dilate the mask arrary according to the structuring element:
SE = strel('square', 10);
maskarray_dilated = imdilate(maskarray, SE);
% Export Image and Mask temporarily for use in N4
niftiwrite(abs(MR+0.001),[parentPath,'Image.nii']); % Image
niftiwrite(abs(maskarray_dilated),[parentPath,'Mask.nii']); % Mask
niftiwrite(abs(maskarray),[parentPath,'Weight.nii']); % Weight

% Run Bias Correction
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

% % Chosen settings:
% cmd = ['"',N4Path,'N4BiasFieldCorrection.exe" ',...%run bias correction
%     '-d 3 ',... % set to 3 dimensions - MUST MATCH DIMENSIONS OF IMAGE
%     '-i "',parentPath,'Image.nii" ',... % input image called of Image.nii
%     '-x "',parentPath,'Mask.nii" ',... % import mask called Mask.nii
%     '-w "',parentPath,'Weight.nii" ',... % import mask called Weight.nii
%     '-s 1 ',... % shrink factor
%     '-c [50,0.001] ',... % convergence
%     '-b [5,3] ',... % spline settings
%     '-t [0.15,0.01,200] ',... % histogram settings
%     '-o ["',parentPath,'CorrectedImage.nii","',parentPath,'Bias.nii"]']; % output corrected image and bias field
% system(cmd);
%% general settings (Abdullah)

% cmd = ['"',N4Path,'N4BiasFieldCorrection.exe"',...%run bias correction
%     ' -d 3 -i "',parentPath,'Image.nii"',... % set to 3 dimensions and input image of Image.nii
%     ' -s 2',... % shrink by factor of 1
%     ' -w "',parentPath,'Weight.nii" ',... % import mask called Weight.nii
%     ' -o ["',parentPath,'CorrectedImage.nii","',parentPath,'Bias.nii"]']; % output corrected image and bias field
% system(cmd);

% cmd = ['"',N4Path,'N4BiasFieldCorrection.exe"',...%run bias correction
%     ' -d 3 -i "',parentPath,'Image.nii"',... % set to 3 dimensions and input image of Image.nii
%     ' -s 1',... % shrink by factor of 1
%     ' -x "',parentPath,'Weight.nii" ',... % import mask called Weight.nii
%     ' -c [100,50,50,50]',... % convergence
%     ' -b [200]',... % spline settings %hf, lr, ap
%     ' -t [1e-10]',... % histogram settings
%     ' -o ["',parentPath,'CorrectedImage.nii","',parentPath,'Bias.nii"]']; % output corrected image and bias field
% system(cmd);

%% new settings (Abdullah)
%Main issue is drop off at top of lungs (dual loop) and RL edges (Polarean)
% Z Correction
cmd = ['"',N4Path,'N4BiasFieldCorrection.exe"',...%run bias correction
    ' -d 3 -i "',parentPath,'Image.nii"',... % set to 3 dimensions and input image of Image.nii
    ' -s 1',... % shrink by factor of 1
    ' -x "',parentPath,'Mask.nii" ',... % import mask called Mask.nii
    ' -w "',parentPath,'Weight.nii" ',... % import mask called Weight.nii
    ' -c [25,0]',... % convergence
    ' -b [1x1x14,3]',... % spline settings %hf, lr, ap
    ' -t [0.75,0.01,100]',... % histogram settings
    ' -o ["',parentPath,'CorrectedImage.nii","',parentPath,'Bias.nii"]']; % output corrected image and bias field
system(cmd);
% X Correction
cmd = ['"',N4Path,'N4BiasFieldCorrection.exe"',...%run bias correction
    ' -d 3 -i "',parentPath,'CorrectedImage.nii"',... % set to 3 dimensions and input image of Image.nii
    ' -s 1',... % shrink by factor of 1
    ' -x "',parentPath,'Mask.nii" ',... % import mask called Mask.nii
    ' -w "',parentPath,'Weight.nii" ',... % import mask called Weight.nii
    ' -c [25,0]',... % convergence
    ' -b [1x14x1,3]',... % spline settings %hf, lr, ap
    ' -t [0.75,0.01,100]',... % histogram settings
    ' -o ["',parentPath,'CorrectedImage.nii","',parentPath,'Bias.nii"]']; % output corrected image and bias field
system(cmd);
% Y Correction
cmd = ['"',N4Path,'N4BiasFieldCorrection.exe"',...%run bias correction
    ' -d 3 -i "',parentPath,'CorrectedImage.nii"',... % set to 3 dimensions and input image of Image.nii
    ' -s 1',... % shrink by factor of 1
    ' -x "',parentPath,'Mask.nii" ',... % import mask called Mask.nii
    ' -w "',parentPath,'Weight.nii" ',... % import mask called Weight.nii
    ' -c [25,0]',... % convergence
    ' -b [14x1x1,3]',... % spline settings %hf, lr, ap
    ' -t [0.75,0.01,100]',... % histogram settings
    ' -o ["',parentPath,'CorrectedImage.nii","',parentPath,'Bias.nii"]']; % output corrected image and bias field
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