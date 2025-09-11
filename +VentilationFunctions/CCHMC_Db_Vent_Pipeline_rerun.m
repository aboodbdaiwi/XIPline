function CCHMC_Db_Vent_Pipeline_rerun(MainInput)

analysisSubfolder = MainInput.analysisFolder;
load(fullfile(MainInput.analysisFolder,'Ventilation_Analysis', 'workspace.mat'));
% mask_file_name = dir(fullfile(analysisSubfolder, 'lungmask_*.nii.gz'));
% mask_file_name = fullfile(mask_file_name.folder, mask_file_name.name);
% [~,~,mask_ext] = fileparts(mask_file_name);
% if strcmp(mask_ext, '.gz') || strcmp(mask_ext, '.nii')
%     try
%         Mask = LoadData.load_nii(mask_file_name);
%     catch
%         Mask = LoadData.load_untouch_nii(mask_file_name);
%     end
%     A = double(Mask.img);
% elseif strcmp(mask_ext, '.dcm')
%     A = double(squeeze(dicomread(mask_file_name)));
% end
% B = A;
B = Ventilation.LungMask;
% B = rot90(B);
% B = flip(B,1);
% B = flip(B,2);
% B = flip(B,1);
% 
B = flip(B,3);
% % B = flip(B,2);
% % figure; imslice(B)
Ventilation.LungMask = B;
figure; imslice(Ventilation.LungMask)

Ventilation.Image = Ventilation.UncorrectedImage;
% B = Ventilation.Image;
% B = flip(B,1);
% B = rot90(B);
% B = flip(B,1); B = flip(B,2);
% Ventilation.Image = B;
% figure; imslice(Ventilation.Image)
% Ventilation.UncorrectedImage = Ventilation.Image;

% B = A;
% B = flip(B,1);
% % % % B = rot90(B);
% B = flip(B,2);
% B = flip(B,3);
% 
delete_if_exist = @(pattern) cellfun(@(f) delete(fullfile(analysisSubfolder, f)), ...
    {dir(fullfile(analysisSubfolder, pattern)).name}, 'UniformOutput', false);

% ==== Delete any older matching files ====
delete_if_exist('Ventilation_Analysis_*.pptx*');

maskarray = double(Ventilation.LungMask + Ventilation.VesselMask);
maskarray(maskarray > 1) = 0;
Ventilation.LungMask = double(maskarray);

% ----------------------------- VDP -----------------------------

myCluster = parcluster('local');
delete(myCluster.Jobs);

[Ventilation] = VentilationFunctions.Ventilation_Analysis(Ventilation, Proton, MainInput);
close all;
Outputs.timestamp = str2double(Ventilation.timestamp);
Outputs.SNR_slice = Ventilation.SNR_slice;
Outputs.SNR = Ventilation.SNR_lung;
Outputs.SNRvv_slice = Ventilation.SNRvv_slice;
Outputs.SNR_vv = Ventilation.SNR_vv;

Outputs.VENT_SERIES_NUMBER = MainInput.xe_sernum;
try
    Outputs.VENTONLINE_SERIES_NUMBER = MainInput.xe_sernum_online;
catch
end
Outputs.ANATVENT_SERIES_NUMBER = MainInput.proton_sernum;
Outputs.IMAGE_ORIENT = MainInput.SliceOrientation;
Outputs.PIXEL_SPACING_X = MainInput.PIXEL_SPACING_X;
Outputs.PIXEL_SPACING_Y = MainInput.PIXEL_SPACING_Y;
Outputs.SLICE_THICKNESS = MainInput.SLICE_THICKNESS;
Outputs.SCAN_NUM = MainInput.SCAN_NUM;
% N4 - 2024 settings (Abood)

if strcmp(Ventilation.HeterogeneityIndex, 'yes')
    Outputs.sliceMeanCV = Ventilation.sliceMeanCV;
    Outputs.sliceVHI = Ventilation.sliceVHI;
    Outputs.overallMeanCV = Ventilation.overallMeanCV;
    Outputs.overallVHI = Ventilation.overallVHI;
end
if strcmp(Ventilation.ThreshAnalysis, 'yes') == 1
    Outputs.TH60.Thresholds =Ventilation.Threshold.Thresholds;
    Outputs.TH60.VDP = Ventilation.Threshold.VDP;
    Outputs.TH60.BinPrecentValues = Ventilation.Threshold.THBins; 
end

if strcmp(Ventilation.DDI2D, 'yes') == 1
    Outputs.TH60.DDI2D.DDI2D_meanSlice = Ventilation.DDI2D_meanSlice;
    Outputs.TH60.DDI2D.DDI2D_stdSlice = Ventilation.DDI2D_stdSlice;
    Outputs.TH60.DDI2D.DDI2D_mean = Ventilation.DDI2D_mean;
    Outputs.TH60.DDI2D.DDI2D_std = Ventilation.DDI2D_std; 
    Outputs.TH60.DDI2D.DDI2D_max = Ventilation.DDI2D_max;
end
if strcmp(Ventilation.DDI3D, 'yes') == 1
    Outputs.TH60.DDI3D.DDI3D_meanSlice = Ventilation.DDI3D_meanSlice;
    Outputs.TH60.DDI3D.DDI3D_stdSlice = Ventilation.DDI3D_stdSlice;
    Outputs.TH60.DDI3D.DDI3D_mean = Ventilation.DDI3D_mean;
    Outputs.TH60.DDI3D.DDI3D_std = Ventilation.DDI3D_std; 
    Outputs.TH60.DDI3D.DDI3D_max = Ventilation.DDI3D_max;
end

if strcmp(Ventilation.GLRLM_Analysis, 'yes') == 1
    Outputs.GLRLM.SRE = Ventilation.SRE;
    Outputs.GLRLM.LRE = Ventilation.LRE;
    Outputs.GLRLM.GLN = Ventilation.GLN;
    Outputs.GLRLM.RLN = Ventilation.RLN;
    Outputs.GLRLM.RP = Ventilation.RP;
    Outputs.GLRLM.LGRE = Ventilation.LGRE;
    Outputs.GLRLM.HGRE = Ventilation.HGRE;
    Outputs.GLRLM.SRLGR = Ventilation.SRLGR;
    Outputs.GLRLM.SRHGE = Ventilation.SRHGE;
    Outputs.GLRLM.LRLGE = Ventilation.LRLGE;
    Outputs.GLRLM.LRHGE = Ventilation.LRHGE;
end

Outputs.Image = double(Ventilation.Image);
Outputs.LungMask = double(Ventilation.LungMask);
Outputs.AirwayMask = double(Ventilation.AirwayMask);
Outputs.VesselMask = double(Ventilation.VesselMask);
Outputs.Mask_Vent_Reg = double(Ventilation.LungMask);
Outputs.AnatImage = Proton.Image;
if strcmp(Ventilation.HeterogeneityIndex, 'yes')
    Outputs.CV_maps = Ventilation.CV_maps;
end
if strcmp(Ventilation.ThreshAnalysis, 'yes') == 1
    Outputs.TH60.VDP_hist = Ventilation.Threshold.VDP_hist;
    Outputs.TH60.DefectArray = Ventilation.Threshold.defectArray;
    Outputs.TH60.VentDefectmap = Ventilation.Threshold.VentDefectmap;    
    Outputs.TH60.vdp_per_slice_local = Ventilation.Threshold.vdp_per_slice_local;
    Outputs.TH60.vdp_per_slice_global = Ventilation.Threshold.vdp_per_slice_global;
end
if strcmp(Ventilation.DDI2D, 'yes') == 1
    Outputs.TH60.DDI2D.DDI2D_DDImap = Ventilation.DDI2D_DDImap;
    Outputs.TH60.DDI2D.DDI2D_Stat = Ventilation.DDI2D_Stat;
end
if strcmp(Ventilation.DDI3D, 'yes') == 1
    Outputs.TH60.DDI3D.DDI3D_DDImap = Ventilation.DDI3D_DDImap;
    Outputs.TH60.DDI3D.DDI3D_Stat = Ventilation.DDI3D_Stat;
end
if strcmp(Ventilation.GLRLM_Analysis, 'yes') == 1
    Outputs.GLRLM.output = Ventilation.GLRLM;
end

Outputs.AnalysisDate = num2str(Outputs.AnalysisDate);
Outputs.ScanDate = num2str(Outputs.ScanDate);
Outputs.timestamp = Ventilation.timestamp;
Outputs.AnalysisVersion = 'v100';
Outputs.skewness = Ventilation.skewness;
Outputs.kurtosis = Ventilation.kurtosis;

OutputJSONFile = fullfile(analysisFolder, ['VentAnalysis_','ser-',num2str(MainInput.xe_sernum),'.json']);
Global.exportStructToJSON(Outputs, OutputJSONFile);

% Check if the directory specified by analysisSubfolder exists, create it if necessary
if ~exist(analysisSubfolder, 'dir')
    mkdir(analysisSubfolder);
end

% Save Outputs to the MAT-file
save(fullfile(analysisSubfolder, 'workspace.mat'));
save(fullfile(analysisSubfolder, 'Ventilation_Analysis_Outputs.mat'), 'Outputs');
disp('all Ventilation Analysis completed');


clearvars
% Global.tts('Ventilation Analysis completed');
end

