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
MainInput.CCHMC_DbDiffAnalysis = 'yes';
MainInput.Scanner = 'Philips'; 

Outputs.Institute = MainInput.Institute;
Outputs.Scanner = MainInput.Scanner;

% get segmentation 
MainInput.SE = 1;
MainInput.thresholdlevel = 0.6;

%%

Outputs.SUBJECT_ID = SubjectID;
Outputs.Age = Age;
Outputs.Sex = Sex;
Outputs.Disease = Disease;
Outputs.ScanDate = string(ScanDate);
Outputs.SW_VER = MainInput.ScannerSoftware;
Outputs.DIFF_FILEPATH_NEW = diff_file(57:end);
Outputs.analysispath = analysisFolder(57:end);
Outputs.maindirectory = diff_file(1:56);
Outputs.ReconType = ReconType;
Outputs.MaskPath = analysisFolder(57:end); 
Outputs.Note = MainInput.Note; 
Outputs.ImageQuality = MainInput.ImageQuality;
Outputs.diff_sernum = MainInput.sernum;
Outputs.AnalysisVersion = 'v100'; %MainInput.analysisversion;

Outputs.AnalysisCode_path = 'https://github.com/aboodbdaiwi/XIPline';
Outputs.AnalysisDate = string(str2double(datestr(datetime('today'), 'yyyymmdd')));
Outputs.timestamp = Outputs.AnalysisDate; 

%%
analysisSubfolder = MainInput.diff_analysis_folder;
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

if strcmp(MainInput.denoiseXe,'yes')
    MainInput.denoisewindow = '[5 5]';
else
    MainInput.denoisewindow = 'NA';
end

Outputs.NoProtonImage = MainInput.NoProtonImage;
Outputs.denoising = MainInput.denoiseXe;
Outputs.DenoiseWindow = MainInput.denoisewindow;

[Ventilation, Diffusion, GasExchange, Proton, MainInput] = LoadData.LoadReadData(MainInput);
Diffusion.Image = double(Diffusion.Image);
Outputs.UncorrectedImage = Diffusion.UncorrectedImage;
% figure; imslice(Diffusion.Image)
Diffusion.outputpath = analysisSubfolder;
Outputs.AnalysisCode_hash = MainInput.AnalysisCode_hash; 

%% Segmentation 
cd(MainInput.XeDataLocation)

% % % diary Log.txt
MainInput.SegmentationMethod = 'Auto'; % 'Threshold' || 'Manual' || 'Auto'
MainInput.SegmentAnatomy = 'Parenchyma'; % 'Airway'; || 'Parenchyma'
MainInput.Imagestosegment = 'Xenon';  % 'Xe & Proton Registered' | 'Xenon' | 'Registered Proton'

MainInput.AIScript = 'Python';
MainInput.thresholdlevel = 0.6; % 'threshold' 
MainInput.SE = 1;

MainInput.SegmentManual = 'Freehand'; % 'AppSegmenter' || 'Freehand'
% MainInput.SliceOrientation = SliceOrientation; % 'coronal' ||'transversal' || 'sagittal' ||'isotropic'
[Proton,Ventilation,Diffusion,GasExchange] = Segmentation.PerformSegmentation(Proton,Ventilation,Diffusion,GasExchange,MainInput);

if ~isfield(Diffusion, 'AirwayMask')
    Diffusion.AirwayMask = zeros(size(Diffusion.LungMask));
end

% figure; imslice(Diffusion.LungMask)
%%

if MainInput.num_b_values == 3
    Diffusion.b_values = [0, 7.5, 15];
    Diffusion.SmallDeltaTime = 3.5; % ms
    Diffusion.BigDeltaTime = 3.5; % ms
    Diffusion.GapTime = 0;
     Diffusion.MorphometryAnalysis = 'no';
elseif MainInput.num_b_values == 4
    Diffusion.b_values = [0, 10, 20, 30];
    Diffusion.SmallDeltaTime = 5.0; % ms
    Diffusion.GapTime = 0; 
    Diffusion.BigDeltaTime = Diffusion.SmallDeltaTime + Diffusion.GapTime; % ms
    Diffusion.MorphometryAnalysis = 'yes';
elseif MainInput.num_b_values == 5
    Diffusion.b_values = [0, 6.25, 12.5, 18.75, 25];
    Diffusion.SmallDeltaTime = 3.5; % ms
    Diffusion.BigDeltaTime = 3.5; % ms
    Diffusion.GapTime = 0;
    Diffusion.MorphometryAnalysis = 'yes';
else
    Diffusion.b_values = NaN;
    Diffusion.SmallDeltaTime = NaN; % ms
    Diffusion.BigDeltaTime = NaN; % ms
    Diffusion.GapTime = NaN;
    disp('Number of bvalues is not either 3, 4, or 5....')
end
% Diffusion.MorphometryAnalysis = 'no';

Outputs.num_b_values = MainInput.num_b_values;
Outputs.b_values = Diffusion.b_values; 
Outputs.SmallDeltaTime = Diffusion.SmallDeltaTime;
Outputs.BigDeltaTime = Diffusion.BigDeltaTime;
Outputs.GapTime = Diffusion.GapTime;
Outputs.MorphometryAnalysis = Diffusion.MorphometryAnalysis;
    
%MainInput.PatientAge = MainInput.Age;
Diffusion.ADCFittingType = 'Log Weighted Linear'; % Log Weighted Linear | Bayesian | Non-Linear | Log Linear
Diffusion.ADCLB_Analysis = 'no';
Diffusion.ADCLB_RefMean = 0.0002*MainInput.Age+0.029; % mean equetion for healthy co. 
Diffusion.ADCLB_RefSD = 5e-5*MainInput.Age+0.0121; 

Outputs.ADCFittingType = Diffusion.ADCFittingType;
Outputs.ADCLB_Analysis = Diffusion.ADCLB_Analysis;
Outputs.ADCLB_RefMean = Diffusion.ADCLB_RefMean; % mean equetion for healthy co. 
Outputs.ADCLB_RefSD = Diffusion.ADCLB_RefSD; 

Diffusion.CMMorphometry = 'yes';
Diffusion.SEMMorphometry = 'yes';
Diffusion.MorphometryAnalysisType = 'human';
Diffusion.Do = 0.14; % 100% Xenon
Diffusion.Delta = Diffusion.BigDeltaTime;  
Diffusion.writereport = 'no';
Outputs.Scanner = '3T-Philips';

Outputs.MorphometryAnalysisType = Diffusion.MorphometryAnalysisType;
Outputs.Do = Diffusion.Do;
% run analysis
[Diffusion] = DiffusionFunctions.Diffusion_Analysis(Diffusion,MainInput);

%% store results
Outputs.Image = Diffusion.Image; 
Outputs.Ndiffimg = Diffusion.Ndiffimg; 
Outputs.final_mask = Diffusion.final_mask; 
Outputs.noise_mask = Diffusion.noise_mask; 
Outputs.ADCmap = Diffusion.ADCmap;
Outputs.ADCcoloredmap = Diffusion.ADCcoloredmap; 

Outputs.SNR_table = Diffusion.SNR_table; 
Outputs.SNR_vec = Diffusion.SNR_vec; 
Outputs.meanADC = Diffusion.meanADC;
Outputs.stdADC = Diffusion.stdADC;
Outputs.ADC_hist = Diffusion.ADC_hist;
Outputs.ADC_cv = Diffusion.ADC_cv;
Outputs.ADC_skewness = Diffusion.ADC_skewness;
Outputs.ADC_kurtosis = Diffusion.ADC_kurtosis;

if strcmp(Diffusion.ADCLB_Analysis, 'yes')
    Outputs.LBADCMean = Diffusion.LBADCMean;
    Outputs.LBADCStd = Diffusion.LBADCStd;
    Outputs.LBADC_hist = Diffusion.LBADC_hist;
    Outputs.LB_BinTable = Diffusion.LB_BinTable;
    Outputs.DiffLow1Percent = Diffusion.DiffLow1Percent;
    Outputs.DiffLow2Percent = Diffusion.DiffLow2Percent;
    Outputs.DiffNormal1Percent = Diffusion.DiffNormal1Percent;
    Outputs.DiffNormal2Percent = Diffusion.DiffNormal2Percent;
    Outputs.DiffHigh1Percent = Diffusion.DiffHigh1Percent;
    Outputs.DiffHigh2Percent = Diffusion.DiffHigh2Percent;
else
    Outputs.LBADCMean = [];
    Outputs.LBADCStd = [];
    Outputs.LBADC_hist = [];
    Outputs.LB_BinTable = [];
    Outputs.DiffLow1Percent = [];
    Outputs.DiffLow2Percent = [];
    Outputs.DiffNormal1Percent = [];
    Outputs.DiffNormal2Percent = [];
    Outputs.DiffHigh1Percent = [];
    Outputs.DiffHigh2Percent = [];
end
if strcmp(Diffusion.MorphometryAnalysis, 'yes') && strcmp(Diffusion.CMMorphometry, 'yes')
    Outputs.R_mean = Diffusion.R_mean;
    Outputs.h_mean = Diffusion.h_mean;
    Outputs.r_mean = Diffusion.r_mean;
    Outputs.Lm_mean = Diffusion.Lm_mean;
    Outputs.SVR_mean = Diffusion.SVR_mean;
    Outputs.Na_mean = Diffusion.Na_mean;
    
    Outputs.R_std = Diffusion.R_std;
    Outputs.h_std = Diffusion.h_std;
    Outputs.r_std = Diffusion.r_std;
    Outputs.Lm_std = Diffusion.Lm_std;
    Outputs.SVR_std = Diffusion.SVR_std;
    Outputs.Na_std = Diffusion.Na_std;

    % store result
    Outputs.R_map = Diffusion.R_map;
    Outputs.h_map = Diffusion.h_map;
    Outputs.r_map = Diffusion.r_map;
    Outputs.Lm_map = Diffusion.Lm_map;
    Outputs.SVR_map = Diffusion.SVR_map;
    Outputs.Na_map = Diffusion.Na_map;
    Outputs.So_map = Diffusion.So_map;
end
if strcmp(Diffusion.MorphometryAnalysis, 'yes') && strcmp(Diffusion.SEMMorphometry, 'yes') 
    Outputs.DDC_mean = Diffusion.DDC_mean;
    Outputs.alpha_mean = Diffusion.alpha_mean;
    Outputs.LmD_mean = Diffusion.LmD_mean;
    Outputs.DDC_std = Diffusion.DDC_std;
    Outputs.alpha_std = Diffusion.alpha_std;
    Outputs.LmD_std = Diffusion.LmD_std;  

    Outputs.DDC_map = Diffusion.DDC_map;
    Outputs.alpha_map = Diffusion.alpha_map; 
    Outputs.SEMSo_map = Diffusion.SEMSo_map; 
    Outputs.LmD_map = Diffusion.LmD_map; 
    Outputs.DDC_mean = Diffusion.DDC_mean; 
    Outputs.alpha_mean = Diffusion.alpha_mean; 
    Outputs.LmD_mean = Diffusion.LmD_mean; 
    Outputs.DDC_std = Diffusion.DDC_std; 
    Outputs.alpha_std = Diffusion.alpha_std; 
    Outputs.LmD_std = Diffusion.LmD_std;
end


%% Save info to JSON file
OutputJSONFile = fullfile(analysisFolder, ['DiffAnalysis_','ser-',num2str(MainInput.sernum),'.json']);
Global.exportStructToJSON(Outputs, OutputJSONFile);

%% Save Outputs to the MAT-file
save(fullfile(analysisSubfolder, 'workspace.mat'));
save(fullfile(analysisSubfolder, 'Diffusion_Analysis_Outputs.mat'), 'Outputs');

clearvars

end


