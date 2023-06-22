function [imag_vol, imag_vol_norm, path, FileNames] = load_data
%% Data loading function for all ventilation input types
% The goal of this function is to upload any data type into MATLAB
% workspace. Feel free to add more data loading methods to the function.
%
% Instructions:
% Run: Functions.load_data
%
% Inputs:
% None.
%
% Outputs:
% imag_vol = image data for multiple slices.
% imag_vol_norm = normalized image data between 0 and 1.
% path = pathname.
% FileNames = filename.
%
% Author: Joseph Plummer.
% Date: 05/05/2021.

%% Select files for 2D or 3D output:
% Clear preixisting file names if needed:
if exist('FileNames')~=0
    clear(FileNames);
end

% Select file:
[FileNames,path] = uigetfile({'*.dcm;*.nii.gz'},'Select files','multiselect','on');
% NOTE: make sure to edit the filters on uigetfilter if you add other file
% loading types.

% Move to file containing data
cd(path)

%% Use switch statement to select between loading options:
% Different file extensions require different loading functions. Add new
% loading functions to the switch-case statements if you want new
% extensions to be loaded. These could include .raw, .lab, .sin, .data,
% .list, etc...

if iscell(FileNames) == 1
    [~,~,ext] = fileparts(FileNames);
    ext = ext{1};
else 
    [~,~,ext] = fileparts(FileNames);
end


switch ext
    % Case statement for .dcm:
    case '.dcm'
        if iscell(FileNames) == 1; % data stored as multiple slices
            
            num_files = length(FileNames);
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
            % Optional: write out a nifti file conversion for later use.
            imag_vol2 = imrotate(imag_vol,90); % niftiwrite likes to rotate and flip, account for this
            imag_vol2 = flipdim(imag_vol2,1);
            niftiwrite(abs(imag_vol2),[path,'Ventilation_MultipleDicomsToSingleNiftiConversion']);
        else % data stored a single file
            
            info = dicominfo([path FileNames]);
            imag_vol = dicomread(info);
            imag_vol = double(squeeze(imag_vol));
            % Optional: write out a nifti file conversion for later use.
            imag_vol2 = imrotate(imag_vol,90); % niftiwrite likes to rotate and flip, account for this
            imag_vol2 = flipdim(imag_vol2,1);
            niftiwrite(abs(imag_vol2),[path,'Ventilation_SingleDicomToSingleNiftiConversion']);
        end
         
        fprintf('.dcm file loaded.\n')
        
        
        % New case statement for .nii.gz file loading:
    case '.gz'
        fprintf('.nii.gz file loaded.\n')
        % Read in images into double format.
        A1 = Functions.load_nii([path FileNames]); % Original
        A = A1.img;
        A = double(squeeze(A));
        
        % The orientation is typically off for the Philips ventilation images.
        % Correct for this here:
        I90 = imrotate(A, -90);
        Ifv = flipdim(I90, 1);
        imag_vol = Ifv;
        
    otherwise
        % Add more cases for other file type loadings.
        disp('Selected file type not supported for upload to this software version.');
        
end

%% Normalize the intensity of the original image to fall between [0,1].
    imag_vol_norm = imag_vol - min(imag_vol(:));
    imag_vol_norm = imag_vol_norm/max(imag_vol_norm(:));

end
