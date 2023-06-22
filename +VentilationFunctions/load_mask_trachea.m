function [maskarraytrachea] = load_mask_trachea
%% Data loading function for the ventilation mask AND trachea 
% The goal of this function is to upload the trachea mask data file in *.nii.gz
% format. This is used for the linear binning ventilation analysis.
%
% Instructions:
% Run: Functions.load_mask_trachea
%
% Inputs:
% None.
%
% Outputs:
% maskarraytrachea = double array containing mask array slices
%
% Author: Joseph Plummer.
% Date: 05/27/2021.

% Open file finder to obtain ventilation mask:
[parentFile2,parentPath2] = uigetfile({'*.nii.gz';'*.nii'}, 'Select Trachea Mask Nifti file');

% Read in mask file to double format:
A2 = Functions.load_nii([parentPath2 parentFile2]);
A3 = A2.img;
A3 = double(squeeze(A3));

% Apply rotations if necessary:
Ifv2 = A3;
I290 = imrotate(A3,90);
Ifv2 = flipdim(I290,1); % Original
maskarraytrachea = Ifv2;

end