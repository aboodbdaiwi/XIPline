
function CCHMC_Db_Diff_Pipeline_rerun(MainInput)

warning('off', 'all'); % Turn off all warnings 'off' | 'on'

UpdatedImageQuality = MainInput.ImageQuality;
UpdatedNote = MainInput.Note;
UpdatedProcessingNotes = MainInput.ProcessingNotes;
UpdatedAnalysisStatus  = MainInput.AnalysisStatus;

analysisSubfolder = MainInput.analysisFolder;
load(fullfile(MainInput.analysisFolder,'Diffusion_Analysis', 'workspace.mat'));

MainInput.ImageQuality = UpdatedImageQuality;
MainInput.Note = UpdatedNote;
MainInput.ProcessingNotes = UpdatedProcessingNotes;
MainInput.AnalysisStatus  = UpdatedAnalysisStatus;

Outputs.ImageQuality = MainInput.ImageQuality;
Outputs.Note = MainInput.Note;
Outputs.ProcessingNotes = MainInput.ProcessingNotes;
Outputs.AnalysisStatus  = MainInput.AnalysisStatus;

mask_file_name = dir(fullfile(analysisSubfolder, 'lungmask_*.nii.gz'));
mask_file_name = fullfile(mask_file_name.folder, mask_file_name.name);
[~,~,mask_ext] = fileparts(mask_file_name);
if strcmp(mask_ext, '.gz') || strcmp(mask_ext, '.nii')
    try
        Mask = LoadData.load_nii(mask_file_name);
    catch
        Mask = LoadData.load_untouch_nii(mask_file_name);
    end
    A = double(Mask.img);
elseif strcmp(mask_ext, '.dcm')
    A = double(squeeze(dicomread(mask_file_name)));
end
Diffusion.LungMask = A;
[CorrectedMask, Diffusion] = VentilationFunctions.correct_mask_orientation(Diffusion);

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

    % store result
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

    % store result
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
    Outputs.DDC_mean = Diffusion.DDC_mean; 
    Outputs.alpha_mean = Diffusion.alpha_mean; 
    Outputs.LmD_mean = Diffusion.LmD_mean; 
    Outputs.DDC_std = Diffusion.DDC_std; 
    Outputs.alpha_std = Diffusion.alpha_std; 
    Outputs.LmD_std = Diffusion.LmD_std;
else
    Outputs.DDC_mean = [];
    Outputs.alpha_mean = [];
    Outputs.LmD_mean = [];
    Outputs.DDC_std = [];
    Outputs.alpha_std = [];
    Outputs.LmD_std =  [];

    Outputs.DDC_map = [];
    Outputs.alpha_map = [];
    Outputs.SEMSo_map = [];
    Outputs.LmD_map = [];
    Outputs.DDC_mean = [];
    Outputs.alpha_mean = [];
    Outputs.LmD_mean = [];
    Outputs.DDC_std = [];
    Outputs.alpha_std = [];
    Outputs.LmD_std = [];
end


%% Save info to JSON file
OutputJSONFile = fullfile(analysisFolder, ['DiffAnalysis_','ser-',num2str(MainInput.sernum),'.json']);
Global.exportStructToJSON(Outputs, OutputJSONFile);

%% Save Outputs to the MAT-file
save(fullfile(analysisSubfolder, 'workspace.mat'));
save(fullfile(analysisSubfolder, 'Diffusion_Analysis_Outputs.mat'), 'Outputs');

clearvars

end


