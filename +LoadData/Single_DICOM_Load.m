function [imag_vol, path, FileName, info] = Single_DICOM_Load(DataPath)
%% Function to load a single DICOM file
% Input:
%   DataPath: Full path to the DICOM file
% Output:
%   imag_vol: 2D image volume
%   path:     Path to the DICOM file
%   FileName: Name of the DICOM file
%   info:     DICOM metadata

% Extract path and filename
[path, FileName, ext] = fileparts(DataPath);
FileName = [FileName, ext];

% Move to file directory
start_path = pwd;
cd(path);

% Load image and metadata
imag_vol = double(squeeze(dicomread(DataPath)));
info = dicominfo(DataPath);
% figure; Global.imslice(imag_vol,'imag_vol')
cd(start_path);
end
