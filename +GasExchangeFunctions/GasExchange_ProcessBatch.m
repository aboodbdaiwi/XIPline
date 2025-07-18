function GasExchange_ProcessBatch(DataFolderLocation,NewMask,NewImages,NewReg,HealthyOnly)
%GasExchange_ProcessBatch - Runs GasExchange_ProcessSingleSubject() for all folders within a folder
%   Requires a folder with the Dissolved_Xe raw/lab/sin files,
%    Dissolved_Xe_Mask raw/lab/sin files, and two (and only two!) sets of 
%    data/list files for the dissolved and proton images
%   If no folder is provided, user will be prompted to select the folder
%
%   Syntax:  GasExchange_ProcessBatch(DataLocation)
%
%   Inputs:
%      DataLocation - (optional)string containing the path of the folder
%      NewMask - (optional) 1 forces new masks to be made, 0 imports masks if present
%      NewImages - (optional) 1 forces new images to be made, 0 imports images if present
%      NewReg - (optional) 1 forces new registration to be made, 0 imports reg tform if present
%      HealthyOnly - (optional) 0 analyzes all files present, 1 analyzes the healthy subjects listed here
%
%   Outputs:
%      none - all data is exported and saved following analysis
%
%   Example: 
%      GasExchange_ProcessBatch()
%      GasExchange_ProcessBatch('C:\Users\mwillmering\Documents\DissolvedData')
%
%   See also: GasExchange_ProcessSingleSubject
% 
%   Package: https://github.com/cchmc-cpir/CCHMC-Gas-Exchange-Processing-Package
% 
%   Author: Matthew Willmering
%   Work email: matthew.willmering@cchmc.org
%   Personal email: matt.willmering@gmail.com
%   Website: https://cpir.cchmc.org/

clearvars -except DataFolderLocation NewMask NewImages NewReg HealthyOnly; clc; close all;
disp('Gas Exchange Batch Processing Starting...')


%% Get Data Folder Location
%Determine folder to use
disp('Looking for necessary files...')
if exist('DataFolderLocation','var') == 0
    DataFolderLocation = uigetdir('C:\', 'Select Folder containing Gas Exchange Data Folders');
end
if exist('NewMask','var') == 0
  NewMask = 1; %1 forces new masks to be made, 0 imports masks if present
end
if exist('NewImages','var') == 0
  NewImages = 1; %1 forces new images to be made, 0 imports images if present
end
if exist('NewReg','var') == 0
  NewReg = 1; %1 forces registration to run
end
if exist('HealthyOnly','var') == 0
  HealthyOnly = 0; %0 forces all subjects to run
end
batch_tic = tic;

%% Determine folders in file
folders = dir(DataFolderLocation); %obtain directory
folders = folders(~ismember({folders.name},{'.','..','All Summary Powerpoints','Scan Info'})); %remove current, parent, Summary folder, Scan Info folder
folders = folders([folders.isdir]'); %remove non-folder files
numFolders = size(folders,1);

%% Modify folders if healthy only
if(HealthyOnly)
    [HealthyCohort, ~] = GasExchange.GetHealthyCohort;

    for ii = numFolders:-1:1%start from end
        folder = strsplit(folders(ii).name); %split file name at spaces

        if(size(folder,2) ~= 2) %only use folder if date and subject, no notes
            folders(ii) = [];%remove from list
            continue;
        else
            subject = folder(1,2);
        end

        Healthy = sum(contains(HealthyCohort,subject)); %determine if matched any healthys
        if (~Healthy)%remove from list
            folders(ii) = [];
        end

    end
    numFolders = size(folders,1);
end


%% Run Process for each Folder
numErrors = 0;
folderErrors = {};
f = waitbar(0,['Running Subject 1 of ', num2str(numFolders)]);
frames = java.awt.Frame.getFrames();%get waitbar fram via Java
frames(end).setAlwaysOnTop(1);%keep on top via Java
movegui(f,[50 600]);
for ii = 1:numFolders
    %update progress bar
    try
        waitbar(ii/numFolders,f,['Running Subject ', num2str(ii), ' of ', num2str(numFolders), ' - (', folders(ii).name,')']);
    catch
        %if waitbar is closed, don't do anything, just continue
    end
    %run single subject
    try
        success = GasExchange_ProcessSingleSubject([folders(ii).folder,'\',folders(ii).name],NewMask,NewImages,NewReg);
        if success == 0 %built in error catch
            numErrors = numErrors + 1;
            folderErrors = [folderErrors, folders(ii).name];
        end
    catch %if function ran into error not accounted for
        numErrors = numErrors + 1;
        folderErrors = [folderErrors, folders(ii).name];
    end
end
%update progress bar
try
    waitbar(1,f,'Finished processing all subjects');
catch
    %if waitbar is closed, don't do anything, just continue
end
close all;


%% Print Errors
disp(['There were ', num2str(numErrors), ' Errors Out Of The ', num2str(numFolders), ' Folders.'])
if numErrors ~= 0
    disp('The following folders had errors occur and were not processed:')
    for ii = 1:numErrors
        disp(folderErrors{1,ii})
    end
end
elapsedTime = toc(batch_tic)/60;
disp(['The processing took ', num2str(elapsedTime), ' minutes to complete for all subjects'])


%% Completed
disp(['Gas Exchange Imaging Process Completed For ', char(DataFolderLocation), '.'])



