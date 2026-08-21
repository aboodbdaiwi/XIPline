
function Outputs = analyzeVentilationImage(imagePath, maskPath)
    % Analyze one ventilation image and one lung mask.
    % Example:
    % Outputs = analyzeVentilationImage('C:\Data\Subject01\xe.nii.gz', ...
    %                                   'C:\Data\Subject01\mask.nii.gz');
    
    dataFolder = fileparts(imagePath);
    outputFolder = fullfile(dataFolder, 'VDP_analysis');
    
    if ~isfolder(outputFolder)
        mkdir(outputFolder);
    end
    
    [MR, imageInfo] = loadNifti(imagePath);
    [lungMask, ~] = loadNifti_mask(maskPath);
    
    MR = double(MR);
    lungMask = double(lungMask > 0);
    
    if ~isequal(size(MR), size(lungMask))
        error('Image and mask sizes do not match.');
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
    
    MainInput.XeFileName = getFileName(imagePath);
    MainInput.XeDataLocation = dataFolder;
    MainInput.XeFullPath = imagePath;
    MainInput.NoProtonImage = 'yes';
    MainInput.SliceOrientation = 'coronal'; 
    Ventilation.ImageResolution = [3,3,15];
    
    Proton.Image = zeros(size(lungMask));
    Proton.ProtonRegistered = zeros(size(lungMask));
    Proton.ProtonRegisteredColored = zeros(size(lungMask));
    
    Ventilation = VentilationFunctions.calculate_SNR(Ventilation);
    close all;
    
    Ventilation = VentilationFunctions.calculate_VHI(Ventilation);
    close all;
    
    Ventilation = VentilationFunctions.calculate_VDP_CCHMC(Ventilation, Proton, MainInput); 
    Outputs.N4.TH60.VDP = Ventilation.Threshold.VDP;
    Outputs.N4.TH60.defectArray = Ventilation.Threshold.defectArray;
    Outputs.N4.TH60.BinPrecentValues = Ventilation.Threshold.THBins; 

    Ventilation.LB_Analysis = 'yes';
    Ventilation.LB_Normalization = 'HybridLBm';
    Ventilation.LBThresholds = [0.5, 0.75, 1.0, 1.22, 1.44]; % 
    Ventilation.Hdist = [-0.893154, 1.133836, 0.285275];
    Ventilation =  VentilationFunctions.calculate_LB_VDP(Ventilation,Proton,MainInput);
    Outputs.N4.HybridGLBm.VDP = Ventilation.LB_VDP;
    Outputs.N4.HybridGLBm.defectArray = Ventilation.VentBinMap2;
    Outputs.N4.HybridGLBm.BinPrecent = Ventilation.BinsPercent;    
    close all;
    
    Outputs.ImagePath = imagePath;
    Outputs.MaskPath = maskPath;
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
    
    save(fullfile(dataFolder, 'Outputs.mat'), 'Outputs', '-v7.3');
    save(fullfile(dataFolder, 'Analysis_Workspace.mat'), '-v7.3');
    
    saveNifti(Outputs.CV_maps, fullfile(outputFolder, 'CV_Map.nii.gz'), imageInfo);
    saveNifti(Outputs.defectArray, fullfile(outputFolder, 'TH60_Defect_Map.nii.gz'), imageInfo);
    
    Results = table( ...
        scalarValue(Outputs.SNR_lung), ...
        scalarValue(Outputs.SNR_vv), ...
        scalarValue(Outputs.overallMeanCV), ...
        scalarValue(Outputs.overallVHI), ...
        scalarValue(Outputs.VDP), ...
        'VariableNames', {'SNR_Lung','SNR_VentilatedVolume','OverallMeanCV','OverallVHI','TH60_VDP'});
    
    writetable(Results, fullfile(dataFolder, 'SNR_VHI_VDP.xlsx'));
    
    fprintf('Analysis complete:\n%s\n', dataFolder);
end

function [volume, info] = loadNifti(filePath)
    try
        nii = LoadData.load_nii(filePath);
    catch
        nii = LoadData.load_untouch_nii(filePath);
    end
    
    volume = double(squeeze(nii.img));
    % volume = permute(volume, [1 3 2]);
    volume = imrotate(volume, -90);
    volume = flip(volume, 2);
    
    try
        info = niftiinfo(filePath);
    catch
        info = [];
    end
end


function [volume, info] = loadNifti_mask(filePath)
    try
        nii = LoadData.load_nii(filePath);
    catch
        nii = LoadData.load_untouch_nii(filePath);
    end
    
    volume = double(squeeze(nii.img));
    % volume = permute(volume, [1 3 2]);
    volume = imrotate(volume, -90);
    volume = flip(volume, 1);
    
    try
        info = niftiinfo(filePath);
    catch
        info = [];
    end
end

function saveNifti(volume, outputPath, info)
    volume = double(squeeze(volume));
    volume = imrotate(flip(volume, 2), -90);
    
    [folder, name, ext] = fileparts(outputPath);
    
    if strcmpi(ext, '.gz')
        [~, name, ~] = fileparts(name);
    end
    
    basePath = fullfile(folder, name);
    
    try
        info.Datatype = 'single';
        info.ImageSize = size(volume);
        niftiwrite(single(volume), basePath, info, 'Compressed', true);
    catch
        niftiwrite(single(volume), basePath, 'Compressed', true);
    end
end

function name = getFileName(filePath)
    [~, name, ext] = fileparts(filePath);
    
    if strcmpi(ext, '.gz')
        [~, name2, ext2] = fileparts(name);
        name = [name2 ext2 ext];
    else
        name = [name ext];
    end
    end
    
    function value = scalarValue(x)
    if isempty(x)
        value = NaN;
    else
        value = mean(double(x(:)), 'omitnan');
    end
end
