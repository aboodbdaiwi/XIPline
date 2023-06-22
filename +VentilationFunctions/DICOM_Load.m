function [imag_vol, path, FileNames] = DICOM_Load(filename_extension);
%% A funtion of read in DICOM Files
% Written by Zackary I. Cleveland 02/26/2018
% 
% Input:    filename_extension: extension of file name, defaults to .dcm
%                               Note, must be a string variable.
%
% Output:   imag_vol:   2D or 3D image output
%           path:       pathway of file containing source data
%           FileNames:  names of files used to generate volume
%%
%% Select files to 2D image or 3D image output

% set file extension
if nargin >=1
    FNE = ['.',filename_extension];
else
    FNE = '.dmc'; % default to DICOM extention for search
end

% clear preixisting file names if needed
if exist('FileNames')~=0
   clear('FileNames');
end

% select file
if nargin >=1 
    [FileNames,path]=uigetfile(['*',FNE],'Select files','multiselect','on');
elseif strcmp(FNE, '.dmc')==1;
    [FileNames,path]=uigetfile('*.dcm','Select DICOM files','multiselect','on');
else        
    [FileNames,path]=uigetfile('multiselect','on'); %load other file extensions
end
 
cd(path) % move to file containing data
% im_info=dicominfo([path FileNames]);
%%
%% Load data and generate output

if iscell(FileNames)==1; % data stored as multiple slices
    
    num_files=length(FileNames);
    all_slices = cell(num_files,1);% load image slices into a cell

    for indx = 1:num_files
        img_slice = dicomread(FileNames{indx});
        all_slices{indx} = img_slice;
    end

    imag_size = size(all_slices{1});
    imag_vol = zeros(imag_size(1),imag_size(2),num_files);

    for slice = 1:num_files
        imag_vol(:,:,slice) = all_slices{slice};
    end
%     im_info=dicominfo([path FileNames]);
else % data stored a single file

    info=dicominfo([path FileNames]);
    imag_vol = dicomread(info);
    imag_vol = double(squeeze(imag_vol));
%     im_info=dicominfo([path FileNames]);
end

end
%%