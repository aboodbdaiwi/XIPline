function CCHMC_Db_GX_Pipeline(MainInput)
%----------------------------- initial inputs -----------------------------
warning('off', 'all'); % Turn off all warnings 'off' | 'on'

SubjectID = MainInput.SubjectID;
Age = MainInput.Age;
Sex = MainInput.Sex;
Disease = MainInput.Disease;
ScanDate = MainInput.ScanDate;
ReconType = MainInput.ReconType;

cal_file = MainInput.cal_file;
gx_file = MainInput.gx_file;
anat_file = MainInput.anat_file;

xedatapath = MainInput.gx_file;
Hdatapath = MainInput.anat_file;
analysisfolder = MainInput.analysisfolder;
MainInput.OutputPath = analysisfolder;

Ventilation = '';
Diffusion = '';
MainInput.AnalysisType = 'GasExchange';
MainInput.Institute = 'CCHMC'; 
MainInput.CCHMC_DbGxAnalysis = 'yes';

% Check if we have anatomical images or not
if isempty(Hdatapath) || strcmp(char(Hdatapath(58:end)), 'NULL')
    MainInput.NoProtonImage = 'yes';  % There is no proton images
else
    MainInput.NoProtonImage = 'no';    % There is  proton images 
end

analysisSubfolder = fullfile(analysisfolder, 'GasExchange_Analysis');
if ~exist(analysisSubfolder, 'dir')
    mkdir(analysisSubfolder);
end
MainInput.analysissesFolder = analysisSubfolder;
MainInput.OutputPath = analysisSubfolder;

% Define file extensions to delete
fileExtensions = {'*.pdf', '*.ppt', '*.pptx'};
% Loop through each file type and delete
for i = 1:numel(fileExtensions)
    files = dir(fullfile(analysisSubfolder, fileExtensions{i}));
    for k = 1:numel(files)
        filePath = fullfile(analysisSubfolder, files(k).name);
        delete(filePath);
        fprintf('Deleted: %s\n', filePath);
    end
end

MainInput.SequenceType = '3D Radial';
if contains(MainInput.Scanner, 'Siemens', 'IgnoreCase', true) 
    MainInput.Scanner = 'Siemens';
elseif contains(MainInput.Scanner, 'Philips', 'IgnoreCase', true) 
    MainInput.Scanner = 'Philips';
elseif contains(MainInput.Scanner, 'GE', 'IgnoreCase', true) 
    MainInput.Scanner = 'GE';    
end

%----------------------------- copy data over -----------------------------
[srcDir, ~] = fileparts(cal_file); 
destDir = fullfile(analysisSubfolder, 'cal_data'); %destination location
% Create destination if it doesn’t exist
if ~exist(destDir, 'dir')
    mkdir(destDir);
end
% Get all files (non-folders) from source
files = dir(srcDir);
files = files(~[files.isdir]); % remove '.' and '..'

% Copy each file
for k = 1:numel(files)
    srcPath  = fullfile(files(k).folder, files(k).name);
    destPath = fullfile(destDir, files(k).name);
    [~,~,ext] = fileparts(destPath);
    if strcmp(ext, '.data')
        MainInput.cal_file = destPath;
    end    
    copyfile(srcPath, destPath);
    fprintf('Copied: %s\n', files(k).name);
end
disp('All files copied successfully!');

[srcDir, ~] = fileparts(gx_file); 
destDir = fullfile(analysisSubfolder, 'gx_data'); % destination location
% Create destination if it doesn’t exist
if ~exist(destDir, 'dir')
    mkdir(destDir);
end
% Get all files (non-folders) from source
files = dir(srcDir);
files = files(~[files.isdir]); % remove '.' and '..'

% Copy each file
for k = 1:numel(files)
    srcPath  = fullfile(files(k).folder, files(k).name);
    destPath = fullfile(destDir, files(k).name);
    [~,~,ext] = fileparts(destPath);
    if strcmp(ext, '.data')
        MainInput.gx_file = destPath;
    end
    copyfile(srcPath, destPath);
    fprintf('Copied: %s\n', files(k).name);
end
disp('All files copied successfully!');
xedatapath = MainInput.gx_file;

if strcmp(MainInput.NoProtonImage, 'no') 
    [srcDir, ~] = fileparts(anat_file); 
    destDir = fullfile(analysisSubfolder, 'anat_data'); % destination location
    % Create destination if it doesn’t exist
    if ~exist(destDir, 'dir')
        mkdir(destDir);
    end
    % Get all files (non-folders) from source
    files = dir(srcDir);
    files = files(~[files.isdir]); % remove '.' and '..'
    
    % Copy each file
    for k = 1:numel(files)
        srcPath  = fullfile(files(k).folder, files(k).name);
        destPath = fullfile(destDir, files(k).name);
        [~,~,ext] = fileparts(destPath);
        if strcmp(ext, '.data')
            MainInput.anat_file = destPath;
        end        
        copyfile(srcPath, destPath);
        fprintf('Copied: %s\n', files(k).name);
    end
    disp('All files copied successfully!');
end
Hdatapath = MainInput.anat_file;

%----------------------------- load data -----------------------------
% Extract the file names
[path,filename] = fileparts(xedatapath);
XeFullPath = xedatapath;
XeDataLocation = path;
[~,XeFileName,xe_ext] = fileparts(xedatapath);

MainInput.XeFullPath = XeFullPath;
MainInput.XeDataLocation = XeDataLocation;
MainInput.XeFileName = filename;
MainInput.XeDataext = xe_ext;
cd(MainInput.XeDataLocation)

% Outputs.XeFullPath = XeFullPath;
Outputs.XeDataLocation = XeDataLocation(57:end);
Outputs.XeFileName = XeFileName;
Outputs.XeDataext = xe_ext;

% select proton
if strcmp(MainInput.NoProtonImage, 'no') == 1
    [path, filename] = fileparts(Hdatapath);
    HFullPath = Hdatapath;
    HDataLocation = path;
    [~,HFileName,H_ext] = fileparts(HFullPath);

    MainInput.HFullPath = HFullPath;
    MainInput.HDataLocation = HDataLocation;
    MainInput.HFileName = filename;
    MainInput.HDataext = H_ext;

    % Outputs.HFullPath = HFullPath;
    Outputs.HDataLocation = HDataLocation(57:end);
    Outputs.HFileName = HFileName;
    Outputs.XeDataext = H_ext;
end
MainInput.denoiseXe = 'no';
[~, ~, GasExchange, Proton, MainInput] = LoadData.LoadReadData(MainInput);
% figure; Global.imslice(double(abs(GasExchange.VentImage)),'Ventilation')
% figure; Global.imslice(Proton.Image,'Proton')

% -----------------------------registration-----------------------------
MainInput.SkipRegistration = 0;
if strcmp(MainInput.NoProtonImage, 'no') 

        % MainInput.RegistrationType = 'ANTs'; 
        % MainInput.TransformType = 'rigid'; % 'translation' | 'rigid' | 'similarity' | 'affine'
        % [Proton] = Registration.PerformRegistration(Proton,Ventilation,GasExchange,MainInput);        

        MainInput.RegistrationType = 'Multimodal'; 
        MainInput.TransformType = 'affine'; % 'translation' | 'rigid' | 'similarity' | 'affine'
        % check if xenon and proton have the same number of slices
        if size(GasExchange.VentImage,3) == size(Proton.Image,3)
            MainInput.SliceSelection = 0;
        else
            MainInput.SliceSelection = 1;
            MainInput.Xestart = 1;
            MainInput.Xeend = size(GasExchange.VentImage,3);
            MainInput.Hstart = 1;
            MainInput.Hend = size(Proton.Image,3);
        end
    
        Xesize = size(GasExchange.VentImage);
        Hsize = size(Proton.Image);
        sizeRatio = Hsize./Xesize;
        
        MainInput.XeVoxelInfo.PixelSize1 = sizeRatio(1);
        MainInput.XeVoxelInfo.PixelSize2 = sizeRatio(2);
        MainInput.XeVoxelInfo.SliceThickness = sizeRatio(3);
        
        MainInput.ProtonVoxelInfo.PixelSize1 = 1;
        MainInput.ProtonVoxelInfo.PixelSize2 = 1;
        MainInput.ProtonVoxelInfo.SliceThickness = 1;
        
        [Proton] = Registration.PerformRegistration(Proton,Ventilation,GasExchange,MainInput);
         % S = orthosliceViewer(Proton.ProtonRegisteredColored);
else
    MainInput.RegistrationType = 'NaN'; 
    MainInput.TransformType = 'NaN'; % 'translation' | 'rigid' | 'similarity' | 'affine'
    Proton.ProtonRegistered = zeros(size(GasExchange.VentImage));
    Proton.ProtonRegisteredColored = zeros(size(GasExchange.VentImage));
end


% -----------------------------segmentation-----------------------------
SegmentMaskMode = 0;    % 0 = new AI mask, 1 = load exisitng mask 
if SegmentMaskMode == 0
    cd(MainInput.XeDataLocation)
    % % % diary Log.txt
    MainInput.SegmentationMethod = 'Auto'; % 'Threshold' || 'Manual' || 'Auto'
    MainInput.SegmentAnatomy = 'Parenchyma'; % 'Airway'; || 'Parenchyma'
    MainInput.Imagestosegment = 'Proton & Xe Registered';  % 'Proton & Xe Registered' | 'Xenon' | 'Registered Proton'
    
    initial_SNR = Global.calcSNR(double(abs(GasExchange.VentImage)));
    if initial_SNR <= 5
        MainInput.Imagestosegment = 'Registered Proton';   
    end
    MainInput.AIScript = 'Python';
    MainInput.thresholdlevel = 1; % 'threshold' 
    MainInput.SE = 1;

    MainInput.SegmentManual = 'Freehand'; % 'AppSegmenter' || 'Freehand'
    MainInput.SliceOrientation = 'isotropic'; % 'coronal' ||'transversal' || 'sagittal' ||'isotropic'
    [Proton,~,~,GasExchange] = Segmentation.PerformSegmentation(Proton,Ventilation,Diffusion,GasExchange,MainInput);
else
  % skip for now
end


% ----------------------------- start analysis -----------------------------
MainInput.ImportHealthyCohort = 0; % 0: import CCHMC healthy Ref., 1: import .mat H. Ref., 2: enter values manually
MainInput.HealthyReferenceType = 'Default';

[GasExchange] = GasExchangeFunctions.GasExchange_Analysis(GasExchange,Proton,MainInput);

clearvars
% Global.tts('Ventilation Analysis completed');
end

