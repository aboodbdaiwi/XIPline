function [maskarray] = load_mask
%% Data loading function for the ventilation mask
% The goal of this function is to upload the mask data file in *.nii.gz
% format.
%
% Instructions:
% Run: Functions.load_mask
%
% Inputs:
% None.
%
% Outputs:
% maskarray = double array containing mask array slices
%
% Author: Joseph Plummer.
% Date: 05/05/2021.

% Open file finder to obtain ventilation mask:
[parentFile2,parentPath2] = uigetfile('*.nii.gz', 'Select Ventilation Mask Nifti file');

% Read in mask file to double format:
A2 = Functions.load_nii([parentPath2 parentFile2]);
A3 = A2.img;
A3 = double(squeeze(A3));
I290 = imrotate(A3,90);
Ifv2 = flipdim(I290,1); % Original
maskarray = Ifv2;

end