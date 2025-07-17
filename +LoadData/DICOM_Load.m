function [imag_vol, path, FileNames, info] = DICOM_Load(DataPath)
%% A function to read in lung masks segmented in Amira and exported as DICOM Files
% 
% Written by Z.I. Cleveland 03/25/2015
% DICOM extension commented out 09/27/2017 ZIC
% 
% Input:    DataPath:   path to the folder containing DICOM files
%
% Output:   imag_vol:   3D image output volume
%           FileNames:  names of all files used to generate volume
%
%% 
save_data = 'no'; % 'save_data';

% Load data
if exist('FileNames', 'var')
    clear('FileNames');
end

start_path = pwd;

% Get list of DICOM files
DataFiles = dir(fullfile(DataPath, '*.dcm'));

DataFiles = struct2cell(DataFiles);
FileNames = DataFiles(1,:);
path = char(DataFiles(2,1));
[~, sortedIdx] = sort(FileNames);
cd(path);

num_files = length(FileNames);
num_slices = num_files;

if num_files > 1    
    msg = sprintf(['\n' int2str(num_slices) ' slices found.']);
    disp(msg);

    % Initialize arrays to store slice positions and instance numbers
    slicePositions = zeros(num_files, 3);
    instanceNumbers = zeros(num_files, 1);
    acquisitionTimes = cell(num_files, 1);

    % Load image slices and their metadata
    all_slices = cell(num_files, 1);

    for indx = 1:num_files
        img_slice = squeeze(dicomread(FileNames{indx}));
        all_slices{indx} = img_slice;
        info = dicominfo(FileNames{indx});        

        % Read Image Position Patient (0020,0032)
        if isfield(info, 'ImagePositionPatient')
            slicePositions(indx, :) = info.ImagePositionPatient;
        end
        
        % Read Instance Number (0020,0013)
        if isfield(info, 'InstanceNumber') && ~isempty(info.InstanceNumber)
            instanceNumbers(indx) = info.InstanceNumber;
        else
            instanceNumbers(indx) = indx;
        end
        
        % Read Acquisition Time (0008,0032)
        if isfield(info, 'AcquisitionTime')
            acquisitionTimes{indx} = info.AcquisitionTime;
        end
    end

    % Generate 3D volume
    imag_size = size(all_slices{1});
    imag_vol = zeros(imag_size(1),imag_size(2),num_slices);

    for slice = 1:num_slices
        if ~isempty(all_slices{slice})
            imag_vol(:,:,slice) = all_slices{slice};
        end
    end

    try 
        % Find unique elements and their counts
        [uniqueVec, ~, idx] = unique(instanceNumbers);
        counts = histc(idx, 1:numel(uniqueVec));
        
        % Find the repeated elements
        repeatedElements = uniqueVec(counts > 1);
        
        % Display the results
        if isempty(repeatedElements)
            disp('No repeated elements found.');
            % Define the new range
            new_min = 1;
            new_max =  length(instanceNumbers);
            % Create a vector of unique integers in the range [new_min, new_max]
            unique_integers = linspace(new_min, new_max, length(instanceNumbers));
            % Sort the original vector and get the sorted indices
            [~, sorted_idx] = sort(instanceNumbers);
            % Assign unique integers based on the sorted order
            unique_integers_sorted = unique_integers(sorted_idx);
            % Create a new vector to store the final unique integers in the original order
            final_unique_integers = zeros(size(instanceNumbers));
            % Place the unique integers in the original order
            final_unique_integers(sorted_idx) = unique_integers;
            imag_sorted = zeros(size(imag_vol));
            for i = 1:size(imag_vol,3)
                imag_sorted(:,:,i) = squeeze(imag_vol(:,:,find(final_unique_integers == i)));
            end            
            % if max(instanceNumbers(:)) > size(imag_vol,3)
            %     for i = 1:size(imag_vol,3)
            %         imag_sorted(:,:,i) = squeeze(imag_vol(:,:,final_unique_integers(i)));
            %     end
            % else
            %     for i = 1:size(imag_vol,3)
            %         imag_sorted(:,:,instanceNumbers(i)) = squeeze(imag_vol(:,:,find(instanceNumbers == i)));
            %     end
            % end
            imag_vol = imag_sorted;
        else
            disp('Repeated elements and their counts:');
            for i = 1:length(repeatedElements)
                fprintf('Element %d is repeated %d times\n', repeatedElements(i), counts(uniqueVec == repeatedElements(i)));
            end
            % Define the 3D array and the vector
    
            % Find the unique numbers and the counts of each number
            uniqueNumbers = unique(instanceNumbers);
            numSlices = length(uniqueNumbers);
            counts = histc(instanceNumbers, uniqueNumbers);
            numBValues = counts(1);  % Assuming the number of repetitions is the same for all unique numbers
            
            % Initialize the new reshaped array
            reshapedArray = zeros(size(imag_vol,1), size(imag_vol,2), numSlices, numBValues);
            
            % Loop through the unique numbers in the vector
            for i = 1:numSlices
                % Get the indices of the current unique number
                indices = find(instanceNumbers == uniqueNumbers(i));
                % Assign the corresponding slices to the reshaped array
                reshapedArray(:, :, i, :) = imag_vol(:, :, indices);
            end
            imag_vol = reshapedArray;
            % Display the size of the reshaped array to verify
            disp(size(reshapedArray));
            if ndims(imag_vol) > 3
                numVolumes = size(imag_vol, 4);
                meanVals = zeros(1, numVolumes);
            
                % Compute mean intensity for each 3D volume along the 4th dimension
                for i = 1:numVolumes
                    vol = squeeze(imag_vol(:,:,:,i));
                    meanVals(i) = mean(vol(:));
                end
            
                % Find the volume with the highest mean
                [~, maxIdx] = max(meanVals);
                imag_vol = squeeze(imag_vol(:,:,:,maxIdx));
            end
        end
    catch
        disp('slices are stored without sorting')
    end

elseif num_files == 1
    imag_vol = squeeze(dicomread(FileNames{1}));  
    info = dicominfo(FileNames{1});
end

% Save data if required
if nargin >= 1
    save_check = strcmp(save_data, 'save_data');
    if save_check == 1 
        save_name = FileNames{1};
        save_name(end-8:end) = [];
        save_name = [save_name '.mat'];
        cd ../        
        save(save_name);
        cd(path);
    else
        msg = sprintf('\nUnrecognized save command. Data was not saved.');
        disp(msg);
    end
end
% imslice(imag_vol)

cd(start_path)


%end
