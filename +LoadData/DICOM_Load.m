function [imag_vol, path, FileNames] = DICOM_Load(DataPath)
%% A funtion of read in lung masks segmented in Amira and exported as DICOM Files
% 
% Written by Z.I. Cleveland 03/25/2015
% DICOM extension commented out 09/27/2017 ZIC
% 
% Input:    num_reps:   number of repetions used in data acquisition
%           save_save:  command to save as a .mat file must be 'save_data'
%                       to work
%
% Output:   imag_vol:   3D image output volume
%           FileNames:  names of all files used to generate volume
%
%% 
%%
save_data = 'no'; % 'save_data';

% Load data
if exist('FileNames')~=0
    clear('FileNames');
end

start_path = pwd;
%[FileNames,path]=uigetfile('*.dcm','Select Dicom files','multiselect','on');
% [FileNames,path]=uigetfile('.dcm','multiselect','on');

DataFiles = dir([DataPath,'\*.dcm']);
DataFiles = struct2cell(DataFiles);
FileNames = DataFiles(1,:);
path = char(DataFiles(2,1));

cd(path);

num_files=length(FileNames);
num_slices= num_files;

if num_files > 1    
    msg = sprintf(['\n' int2str(num_slices) ' slices found.']);
    disp(msg);

    % load image slices into a cell
    all_slices = cell(num_files,1);

    for indx = 1:num_files
        img_slice = squeeze(dicomread(FileNames{indx}));
        all_slices{indx} = img_slice;
    end

    % generate 3D volume
    imag_size = size(all_slices{1});
    imag_vol = zeros(imag_size(1),imag_size(2),num_slices);
    % Check if the variable is empty
    if isempty(imag_vol)
%         nonEmptyCount = sum(~cellfun('isempty', all_slices));
        imag_size = size(all_slices{num_files});
        imag_vol = zeros(imag_size(1),imag_size(2),num_files);       
    end
    for slice = 1:num_slices
        if ~isempty(all_slices{slice})
            imag_vol(:,:,slice) = all_slices{slice};
        end
    end
elseif num_files == 1
    imag_vol = squeeze(dicomread(FileNames{1}));     
end

% Save data
if nargin >=1
    
    save_check = strcmp(save_data,'save_data');
    
    if save_check == 1 
        save_name = FileNames{1};
        save_name(end-8:end)=[];
        save_name = [save_name '.mat'];
        cd ../        
        save(save_name);
        cd(path);
    else
        msg = sprintf('\nUnrecognized save command. Data was not saved.');
        disp(msg);
    end
end

cd(start_path)

%end