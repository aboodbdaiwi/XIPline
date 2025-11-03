
function CCHMC_Db_Vent_Pipeline_NoData(MainInput)


analysisSubfolder = fullfile(MainInput.analysisFolder, 'Ventilation_Analysis');
if ~exist(analysisSubfolder, 'dir')
    mkdir(analysisSubfolder);
end

Outputs = [];
Outputs.SUBJECT_ID = MainInput.SubjectID;
Outputs.Age = MainInput.Age;
Outputs.Sex = MainInput.Sex;
Outputs.Disease = MainInput.Disease;
Outputs.ScanDate = MainInput.ScanDate;
Outputs.SW_VER = MainInput.ScannerSoftware;
Outputs.VENT_FILEPATH_NEW = MainInput.vent_file;
Outputs.ANATVENT_FILEPATH_NEW = MainInput.anat_file;
Outputs.analysispath = MainInput.analysisFolder(57:end);
Outputs.maindirectory = MainInput.vent_file;
Outputs.ReconType = MainInput.ReconType;
Outputs.MaskPath = ''; 
Outputs.Offset_Frequency = MainInput.freqoffset;
Outputs.ImageQuality = MainInput.ImageQuality;
Outputs.Note = MainInput.Note;
Outputs.xe_sernum = MainInput.xe_sernum;
Outputs.proton_sernum = MainInput.proton_sernum;
Outputs.AnalysisCode_path = 'https://github.com/aboodbdaiwi/XIPline';
Outputs.AnalysisDate = str2double(datestr(datetime('today'), 'yyyymmdd'));
Outputs.XeDataLocation = MainInput.vent_file;
Outputs.XeFileName = '';
Outputs.XeDataext = '';
Outputs.HDataLocation = MainInput.anat_file;
Outputs.HFileName = '';
Outputs.XeDataext = '';
Outputs.AnalysisCode_hash = ''; 
Outputs.RegistrationType = '';
Outputs.TransformType = '';
Outputs.timestamp = str2double(Outputs.AnalysisDate);
Outputs.SNR_slice = '';
Outputs.SNR = '';
Outputs.SNRvv_slice = '';
Outputs.SNR_vv = '';
Outputs.VENT_SERIES_NUMBER = MainInput.xe_sernum;
Outputs.VENTONLINE_SERIES_NUMBER = '';
Outputs.ANATVENT_SERIES_NUMBER = MainInput.proton_sernum;
Outputs.IMAGE_ORIENT = MainInput.SliceOrientation;
Outputs.PIXEL_SPACING_X = MainInput.PIXEL_SPACING_X;
Outputs.PIXEL_SPACING_Y = MainInput.PIXEL_SPACING_Y;
Outputs.SLICE_THICKNESS = MainInput.SLICE_THICKNESS;
Outputs.SCAN_NUM = MainInput.SCAN_NUM;
Outputs.sliceMeanCV = '';
Outputs.sliceVHI = '';
Outputs.overallMeanCV = '';
Outputs.overallVHI = '';
Outputs.TH60.Thresholds = '';
Outputs.TH60.VDP = '';
Outputs.TH60.BinPrecentValues = ''; 
Outputs.TH60.DDI2D.DDI2D_meanSlice = '';
Outputs.TH60.DDI2D.DDI2D_stdSlice = '';
Outputs.TH60.DDI2D.DDI2D_mean = '';
Outputs.TH60.DDI2D.DDI2D_std = '';
Outputs.TH60.DDI2D.DDI2D_max = '';
Outputs.GLRLM.SRE = '';
Outputs.GLRLM.LRE = '';
Outputs.GLRLM.GLN = '';
Outputs.GLRLM.RLN = '';
Outputs.GLRLM.RP = '';
Outputs.GLRLM.LGRE = '';
Outputs.GLRLM.HGRE = '';
Outputs.GLRLM.SRLGR = '';
Outputs.GLRLM.SRHGE = '';
Outputs.GLRLM.LRLGE = '';
Outputs.GLRLM.LRHGE = '';
Outputs.Image = '';
Outputs.LungMask = '';
Outputs.AirwayMask = '';
Outputs.VesselMask = '';
Outputs.Mask_Vent_Reg = '';
Outputs.AnatImage = '';
Outputs.CV_maps = '';
Outputs.TH60.VDP_hist = '';
Outputs.TH60.DefectArray = '';
Outputs.TH60.VentDefectmap = '';    
Outputs.TH60.vdp_per_slice_local = '';
Outputs.TH60.vdp_per_slice_global = '';
Outputs.TH60.DDI2D.DDI2D_DDImap = '';
Outputs.TH60.DDI2D.DDI2D_Stat = '';
Outputs.GLRLM.output = '';
Outputs.AnalysisDate = num2str(Outputs.AnalysisDate);
Outputs.ScanDate = num2str(Outputs.ScanDate);
Outputs.timestamp = '';
Outputs.AnalysisVersion = 'v100';
Outputs.skewness = '';
Outputs.kurtosis = '';

OutputJSONFile = fullfile(MainInput.analysisFolder, ['VentAnalysis_','ser-',num2str(MainInput.xe_sernum),'.json']);
Global.exportStructToJSON(Outputs, OutputJSONFile);


clearvars

end

