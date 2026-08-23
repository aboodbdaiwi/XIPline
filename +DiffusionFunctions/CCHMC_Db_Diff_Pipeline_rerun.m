function CCHMC_Db_Diff_Pipeline_rerun(MainInput)
%CCHMC_DB_DIFF_PIPELINE_RERUN Rerun diffusion analysis using cached workspace
% and corrected mask.

warning('off', 'all');

UpdatedImageQuality = MainInput.ImageQuality;
UpdatedNote         = MainInput.Note;
UpdatedMaskPath     = MainInput.MaskPath;

analysisFolder    = MainInput.analysisfolder;
% analysisSubfolder = MainInput.diff_analysis_folder;
analysisSubfolder = MainInput.analysisfolder;


workspaceFile = fullfile(analysisSubfolder, 'workspace.mat');

if ~isfile(workspaceFile)
    error("DiffusionFunctions:RerunMissingWorkspace", ...
        "Cannot rerun diffusion analysis because workspace.mat was not found: %s", ...
        workspaceFile);
end

load(workspaceFile);

analysisSubfolder = MainInput.analysisfolder;

% Restore updated inputs after loading old workspace
MainInput.ImageQuality = UpdatedImageQuality;
MainInput.Note         = UpdatedNote;
MainInput.MaskPath     = UpdatedMaskPath;
MainInput.OutputPath   = analysisFolder;
MainInput.analysissesFolder = analysisSubfolder;

if exist('Outputs','var')
    Outputs.ImageQuality = MainInput.ImageQuality;
    Outputs.Note         = MainInput.Note;
    Outputs.MaskPath     = MainInput.MaskPath;
end

% Load corrected mask
maskFile = fullfile(analysisSubfolder, 'LungMask.nii.gz' );
A = loadMaskAsDouble(maskFile);

A = squeeze(A);

if ndims(A) == 4
    A = rot90(A(:,:,:,1), 2);
end
% Replace mask in cached Diffusion structure
if ~exist('Diffusion','var')
    error("DiffusionFunctions:RerunMissingDiffusionStruct", ...
        "workspace.mat does not contain variable 'Diffusion'.");
end

A = squeeze(A);
Diffusion.LungMask = double(A);

if ~isfield(Diffusion, 'AirwayMask') || isempty(Diffusion.AirwayMask)
    Diffusion.AirwayMask = zeros(size(Diffusion.LungMask));
end

% Optional sanity check
imgSize = size(Diffusion.Image);
maskSize = size(Diffusion.LungMask);

if numel(imgSize) >= 3
    expectedMaskSize = imgSize(1:3);
else
    expectedMaskSize = imgSize;
end

if ~isequal(maskSize, expectedMaskSize)
    warning("DiffusionFunctions:RerunMaskSizeMismatch", ...
        "Mask size [%s] does not match image spatial size [%s]. Check orientation.", ...
        num2str(maskSize), num2str(expectedMaskSize));
end

% Delete old reports so they regenerate
deleteFilesMatching(analysisSubfolder, '*.pdf');
deleteFilesMatching(analysisSubfolder, '*.ppt');
deleteFilesMatching(analysisSubfolder, '*.pptx');

% Rerun diffusion analysis with corrected mask
Diffusion.outputpath = analysisSubfolder;
Diffusion.ADCFittingType = 'Log Weighted Linear'; % Log Weighted Linear | Bayesian | Non-Linear | Log Linear
Diffusion.ADCLB_Analysis = 'yes';
Age = str2double(string(MainInput.Age));
Diffusion.ADCLB_RefMean = 0.0003 * Age + 0.024;
Diffusion.ADCLB_RefSD = 2e-5*Age+0.0073;
Diffusion = DiffusionFunctions.Diffusion_Analysis(Diffusion, MainInput);


% Rebuild Outputs using the same fields as the full pipeline
Outputs.ImageQuality = MainInput.ImageQuality;
Outputs.Note         = MainInput.Note;
Outputs.MaskPath     = MainInput.MaskPath;

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
    Outputs.AcinarRadius_mean = Diffusion.R_mean;
    Outputs.h_mean = Diffusion.h_mean;
    Outputs.AlveolarRadius_mean = Diffusion.r_mean;
    Outputs.Lm_mean = Diffusion.Lm_mean;
    Outputs.SVR_mean = Diffusion.SVR_mean;
    Outputs.Na_mean = Diffusion.Na_mean;
    
    Outputs.AcinarRadius_std = Diffusion.R_std;
    Outputs.h_std = Diffusion.h_std;
    Outputs.AlveolarRadius_std = Diffusion.r_std;
    Outputs.Lm_std = Diffusion.Lm_std;
    Outputs.SVR_std = Diffusion.SVR_std;
    Outputs.Na_std = Diffusion.Na_std;

    Outputs.AcinarRadius_map = Diffusion.R_map;
    Outputs.h_map = Diffusion.h_map;
    Outputs.AlveolarRadius_map = Diffusion.r_map;
    Outputs.Lm_map = Diffusion.Lm_map;
    Outputs.SVR_map = Diffusion.SVR_map;
    Outputs.Na_map = Diffusion.Na_map;
    Outputs.So_map = Diffusion.So_map;
else
    Outputs.AcinarRadius_mean = [];
    Outputs.h_mean = [];
    Outputs.AlveolarRadius_mean = [];
    Outputs.Lm_mean = [];
    Outputs.SVR_mean = [];
    Outputs.Na_mean = [];
    
    Outputs.AcinarRadius_std = [];
    Outputs.h_std = [];
    Outputs.AlveolarRadius_std = [];
    Outputs.Lm_std = [];
    Outputs.SVR_std = [];
    Outputs.Na_std = [];

    Outputs.AcinarRadius_map = [];
    Outputs.h_map = [];
    Outputs.AlveolarRadius_map = [];
    Outputs.Lm_map = [];
    Outputs.SVR_map = [];
    Outputs.Na_map = [];
    Outputs.So_map = [];
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
else
    Outputs.DDC_mean = [];
    Outputs.alpha_mean = [];
    Outputs.LmD_mean = [];
    Outputs.DDC_std = [];
    Outputs.alpha_std = [];
    Outputs.LmD_std = [];

    Outputs.DDC_map = [];
    Outputs.alpha_map = [];
    Outputs.SEMSo_map = [];
    Outputs.LmD_map = [];
end

% Save JSON and MAT outputs
OutputJSONFile = fullfile(analysisFolder, ...
    ['DiffAnalysis_', 'ser-', num2str(MainInput.sernum), '.json']);

Global.exportStructToJSON(Outputs, OutputJSONFile);

save(fullfile(analysisSubfolder, 'workspace.mat'));
save(fullfile(analysisSubfolder, 'Diffusion_Analysis_Outputs.mat'), 'Outputs');

disp('Diffusion rerun with corrected mask completed.');

end



function A = loadMaskAsDouble(maskFile)

[~,~,ext] = fileparts(maskFile);

if strcmpi(ext, '.gz') || strcmpi(ext, '.nii')
    try
        Mask = LoadData.load_nii(maskFile);
    catch
        Mask = LoadData.load_untouch_nii(maskFile);
    end
    A = flipud(rot90(double(Mask.img)));

elseif strcmpi(ext, '.dcm')
    A = double(squeeze(dicomread(maskFile)));

else
    error("DiffusionFunctions:RerunUnsupportedMaskType", ...
        "Unsupported mask file type: %s", maskFile);
end

end


function deleteFilesMatching(folderName, pattern)

files = dir(fullfile(folderName, pattern));

for k = 1:numel(files)
    f = fullfile(files(k).folder, files(k).name);
    try
        delete(f);
        fprintf('Deleted old output: %s\n', f);
    catch ME
        warning("DiffusionFunctions:DeleteFailed", ...
            "Could not delete old output %s\nReason: %s", f, ME.message);
    end
end

end