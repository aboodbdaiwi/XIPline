
function CCHMC_Db_Diff_Pipeline_NoData(MainInput)
    
    analysisSubfolder = fullfile(MainInput.analysisfolder, 'Diffusion_Analysis');
    if ~exist(analysisSubfolder, 'dir')
        mkdir(analysisSubfolder);
    end
    
    Outputs = [];
    Outputs.SUBJECT_ID = MainInput.SubjectID;
    Outputs.Age = MainInput.Age;
    Outputs.Sex = MainInput.Sex;
    Outputs.Disease = MainInput.Disease;
    Outputs.ScanDate = string(MainInput.ScanDate);
    Outputs.SW_VER = MainInput.ScannerSoftware;
    Outputs.DIFF_FILEPATH_NEW = MainInput.diff_file(57:end);
    Outputs.analysispath = MainInput.analysisfolder(57:end);
    Outputs.maindirectory = MainInput.diff_file(1:56);
    Outputs.ReconType = MainInput.ReconType;
    Outputs.MaskPath = ''; 
    Outputs.Note = MainInput.Note; 
    Outputs.ImageQuality = MainInput.ImageQuality;
    Outputs.diff_sernum = MainInput.sernum;
    Outputs.AnalysisVersion = 'v100'; %MainInput.analysisversion;
    
    Outputs.AnalysisCode_path = 'https://github.com/aboodbdaiwi/XIPline';
    Outputs.AnalysisDate = string(str2double(datestr(datetime('today'), 'yyyymmdd')));
    Outputs.timestamp = Outputs.AnalysisDate; 
    Outputs.Institute = 'CCHMC';
    Outputs.Scanner = MainInput.Scanner;
    
    Outputs.XeDataLocation = Outputs.analysispath;
    [~, XeFileName, xe_ext] = fileparts(MainInput.diff_file); 
    Outputs.XeFileName = XeFileName;
    Outputs.XeDataext = xe_ext;
    
    Outputs.denoising = MainInput.denoiseXe;
    Outputs.DenoiseWindow = '[5 5]';
    
    Outputs.UncorrectedImage = '';
    Outputs.AnalysisCode_hash = ''; 
    
    Outputs.num_b_values = '';
    Outputs.b_values = '';
    Outputs.SmallDeltaTime = '';
    Outputs.BigDeltaTime = '';
    Outputs.GapTime = '';
    Outputs.MorphometryAnalysis = '';
    
    Outputs.ADCFittingType = '';
    Outputs.ADCLB_Analysis = '';
    Outputs.ADCLB_RefMean = ''; % mean equetion for healthy co. 
    Outputs.ADCLB_RefSD = ''; 
    Outputs.Do = 0.14;
    
    Outputs.Image = '';
    Outputs.Ndiffimg = '';
    Outputs.final_mask = '';
    Outputs.noise_mask = '';
    Outputs.ADCmap = '';
    Outputs.ADCcoloredmap = '';
    
    Outputs.SNR_table = '';
    Outputs.SNR_vec = ''; 
    Outputs.meanADC = '';
    Outputs.stdADC = '';
    Outputs.ADC_hist = '';
    Outputs.ADC_cv = '';
    Outputs.ADC_skewness = '';
    Outputs.ADC_kurtosis = '';
    
    Outputs.LBADCMean = '';
    Outputs.LBADCStd = '';
    Outputs.LBADC_hist = '';
    Outputs.LB_BinTable = '';
    Outputs.DiffLow1Percent = '';
    Outputs.DiffLow2Percent = '';
    Outputs.DiffNormal1Percent = '';
    Outputs.DiffNormal2Percent = '';
    Outputs.DiffHigh1Percent = '';
    Outputs.DiffHigh2Percent = '';
    
    Outputs.AcinarRadius_mean = '';
    Outputs.h_mean ='';
    Outputs.AlveolarRadius_mean = '';
    Outputs.Lm_mean ='';
    Outputs.SVR_mean = '';
    Outputs.Na_mean = '';
    
    Outputs.AcinarRadius_std = '';
    Outputs.h_std = '';
    Outputs.AlveolarRadius_std = '';
    Outputs.Lm_std = '';
    Outputs.SVR_std = '';
    Outputs.Na_std ='';

    Outputs.AcinarRadius_map ='';
    Outputs.h_map = '';
    Outputs.AlveolarRadius_map = '';
    Outputs.Lm_map = '';
    Outputs.SVR_map = '';
    Outputs.Na_map = '';
    Outputs.So_map = '';
    
    Outputs.DDC_mean =  '';
    Outputs.alpha_mean = '';
    Outputs.LmD_mean = '';
    Outputs.DDC_std = '';
    Outputs.alpha_std =  '';
    Outputs.LmD_std =  '';
    
    Outputs.DDC_map = '';
    Outputs.alpha_map = '';
    Outputs.SEMSo_map = '';
    Outputs.LmD_map = '';
    Outputs.DDC_mean = '';
    Outputs.alpha_mean = '';
    Outputs.LmD_mean = '';
    Outputs.DDC_std = '';
    Outputs.alpha_std = '';
    Outputs.LmD_std = '';        
    
    OutputJSONFile = fullfile(MainInput.analysisfolder, ['DiffAnalysis_','ser-',num2str(MainInput.sernum),'.json']);
    Global.exportStructToJSON(Outputs, OutputJSONFile);
    
    clearvars

end

