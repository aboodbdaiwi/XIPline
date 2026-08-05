function Outputs = analyzeVentilationImage2(MR, lungMask, dataPath)

MR = double(squeeze(MR));
lungMask = double(squeeze(lungMask) > 0);

if ~isequal(size(MR), size(lungMask))
    error('Image and mask sizes do not match.');
end

if ~isfolder(dataPath)
    mkdir(dataPath);
end

outputFolder = fullfile(dataPath, 'VDP_analysis');

if ~isfolder(outputFolder)
    mkdir(outputFolder);
end

Ventilation.Image = MR;
Ventilation.UncorrectedImage = MR;
Ventilation.Image_uncorrected = MR;
Ventilation.LungMask = lungMask;
Ventilation.AirwayMask = zeros(size(lungMask));
Ventilation.VesselMask = zeros(size(lungMask));
Ventilation.parentPath = outputFolder;
Ventilation.HeterogeneityIndex = 'yes';
Ventilation.ThreshAnalysis = 'yes';
Ventilation.N4Analysis = 'no';
Ventilation.IncompleteThresh = 60;
Ventilation.CompleteThresh = 30;
Ventilation.HyperventilatedThresh = 200;
Ventilation.MedianFilter = 'yes';
Ventilation.ImageResolution = [3, 3, 15];

MainInput.XeFileName = '';
MainInput.XeDataLocation = dataPath;
MainInput.XeFullPath = '';
MainInput.NoProtonImage = 'yes';
MainInput.SliceOrientation = 'coronal';

Proton.Image = zeros(size(lungMask));
Proton.ProtonRegistered = zeros(size(lungMask));
Proton.ProtonRegisteredColored = zeros(size(lungMask));

Ventilation = VentilationFunctions.calculate_SNR(Ventilation);
close all;

Ventilation = VentilationFunctions.calculate_VHI(Ventilation);
close all;

Ventilation = VentilationFunctions.calculate_VDP_CCHMC( ...
    Ventilation, Proton, MainInput);
close all;

Outputs.N4.TH60.VDP = Ventilation.Threshold.VDP;
Outputs.N4.TH60.defectArray = Ventilation.Threshold.defectArray;
Outputs.N4.TH60.BinPrecentValues = Ventilation.Threshold.THBins;

Ventilation.LB_Analysis = 'yes';
Ventilation.LB_Normalization = 'HybridLBm';
Ventilation.LBThresholds = [0.5, 0.75, 1.0, 1.22, 1.45];
Ventilation.Hdist = [-0.893154, 1.133836, 0.285275];

Ventilation = VentilationFunctions.calculate_LB_VDP( ...
    Ventilation, Proton, MainInput);
close all;

Outputs.N4.HybridGLBm.VDP = Ventilation.LB_VDP;
Outputs.N4.HybridGLBm.defectArray = Ventilation.VentBinMap2;
Outputs.N4.HybridGLBm.BinPrecent = Ventilation.BinsPercent;

Outputs.DataPath = dataPath;
Outputs.Image = MR;
Outputs.LungMask = lungMask;

Outputs.SNR_slice = Ventilation.SNR_slice;
Outputs.SNRvv_slice = Ventilation.SNRvv_slice;
Outputs.SNR_lung = Ventilation.SNR_lung;
Outputs.SNR_vv = Ventilation.SNR_vv;

Outputs.CV_maps = Ventilation.CV_maps;
Outputs.sliceMeanCV = Ventilation.sliceMeanCV;
Outputs.sliceVHI = Ventilation.sliceVHI;
Outputs.overallMeanCV = Ventilation.overallMeanCV;
Outputs.overallVHI = Ventilation.overallVHI;

Outputs.VDP = Ventilation.Threshold.VDP;
Outputs.defectArray = Ventilation.Threshold.defectArray;
Outputs.BinPercentValues = Ventilation.Threshold.THBins;

save(fullfile(dataPath, 'Outputs.mat'), 'Outputs', '-v7.3');
save(fullfile(dataPath, 'Analysis_Workspace.mat'), '-v7.3');

niftiwrite(single(Outputs.CV_maps), ...
    fullfile(outputFolder, 'CV_Map.nii'), ...
    'Compressed', true);

niftiwrite(single(Outputs.N4.TH60.defectArray), ...
    fullfile(outputFolder, 'TH60_Defect_Map.nii'), ...
    'Compressed', true);

niftiwrite(single(Outputs.N4.HybridGLBm.defectArray), ...
    fullfile(outputFolder, 'HybridGLBm_Defect_Map.nii'), ...
    'Compressed', true);

Results = table( ...
    scalarValue(Outputs.SNR_lung), ...
    scalarValue(Outputs.SNR_vv), ...
    scalarValue(Outputs.overallMeanCV), ...
    scalarValue(Outputs.overallVHI), ...
    scalarValue(Outputs.N4.TH60.VDP), ...
    scalarValue(Outputs.N4.HybridGLBm.VDP), ...
    'VariableNames', { ...
    'SNR_Lung', ...
    'SNR_VentilatedVolume', ...
    'OverallMeanCV', ...
    'OverallVHI', ...
    'TH60_VDP', ...
    'HybridGLBm_VDP'});

writetable(Results, fullfile(dataPath, 'SNR_VHI_VDP.xlsx'));

fprintf('Analysis complete:\n%s\n', dataPath);

end

function value = scalarValue(x)

if isempty(x)
    value = NaN;
else
    value = mean(double(x(:)), 'omitnan');
end

end