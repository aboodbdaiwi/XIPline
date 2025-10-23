function CCHMC_Db_Diff_Pipeline(MainInput)

warning('off', 'all'); % Turn off all warnings 'off' | 'on'

SubjectID = MainInput.SubjectID;
Age = MainInput.Age;
Sex = MainInput.Sex;
Disease = MainInput.Disease;
ScanDate = MainInput.ScanDate;
ReconType = MainInput.ReconType;

diff_file = MainInput.diff_file;

analysisFolder = MainInput.analysisfolder;

MainInput.OutputPath = analysisFolder;

xedatapath = diff_file;

Ventilation = '';
GasExchange = '';
MainInput.AnalysisType = 'Diffusion';
MainInput.Institute = 'CCHMC'; 

% get segmentation 
MainInput.SE = 1;
MainInput.thresholdlevel = 0.6;

%%

Outputs = [];
Outputs.SUBJECT_ID = SubjectID;
Outputs.Age = Age;
Outputs.Sex = Sex;
Outputs.Disease = Disease;
Outputs.ScanDate = ScanDate;
Outputs.SW_VER = MainInput.ScannerSoftware;
Outputs.DIFF_FILEPATH_NEW = diff_file(57:end);
Outputs.analysispath = analysisFolder(57:end);
Outputs.maindirectory = diff_file(1:56);
Outputs.ReconType = ReconType;
Outputs.MaskPath = analysisFolder(57:end); 
Outputs.Note = MainInput.Note;
Outputs.sernum = MainInput.sernum;

Outputs.AnalysisCode_path = 'https://github.com/aboodbdaiwi/XIPline';
Outputs.AnalysisDate = str2double(datestr(datetime('today'), 'yyyymmdd'));
AnalysisDate = Outputs.AnalysisDate;

%%
analysisSubfolder = fullfile(analysisFolder, 'Ventilation_Analysis');
if ~exist(analysisSubfolder, 'dir')
    mkdir(analysisSubfolder);
end
MainInput.analysissesFolder = analysisSubfolder;

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


%% Load Data
MainInput.NoProtonImage = 1;
MainInput.denoiseXe = 'no';
[Ventilation, Diffusion, GasExchange, Proton, MainInput] = LoadData.LoadReadData(MainInput);
Diffusion.Image = double(Diffusion.Image);

Diffusion.outputpath = analysisSubfolder;
Outputs.AnalysisCode_hash = MainInput.AnalysisCode_hash; 

%% Reorient Diffusion Images if necessary

%%


%%


Image_to_Segment = squeeze(Diffusion.Image(:,:,:,1));
LungMask = Segmentation.SegmentLungthresh(Image_to_Segment,MainInput.SE,MainInput.thresholdlevel);

try
    airway_mask = Diffusion.AirwayMask;
catch
    airway_mask = zeros(size(lung_mask));
end

%%
cd(MainInput.XeDataLocation)

    % % % diary Log.txt
    MainInput.SegmentationMethod = 'Auto'; % 'Threshold' || 'Manual' || 'Auto'
    MainInput.SegmentAnatomy = 'Parenchyma'; % 'Airway'; || 'Parenchyma'
    MainInput.Imagestosegment = 'Xenon';  % 'Xe & Proton Registered' | 'Xenon' | 'Registered Proton'

   
    MainInput.AIScript = 'Python';
    MainInput.thresholdlevel = 1; % 'threshold' 
    MainInput.SE = 1;

    MainInput.SegmentManual = 'Freehand'; % 'AppSegmenter' || 'Freehand'
    % MainInput.SliceOrientation = SliceOrientation; % 'coronal' ||'transversal' || 'sagittal' ||'isotropic'
    [Proton,Ventilation,Diffusion,GasExchange] = Segmentation.PerformSegmentation(Proton,Ventilation,Diffusion,GasExchange,MainInput);

    if exist('Ventilation.AirwayMask', 'var')
        % skip
    else
        Diffusion.AirwayMask = zeros(size(Ventilation.LungMask));
    end

