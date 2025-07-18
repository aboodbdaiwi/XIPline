function CCHMC_Db_Vent_Pipeline(MainInput)
%----------------------------- initial inputs -----------------------------
warning('off', 'all'); % Turn off all warnings 'off' | 'on'

SubjectID = MainInput.SubjectID;
Age = MainInput.Age;
Sex = MainInput.Sex;
Disease = MainInput.Disease;
ScanDate = MainInput.ScanDate;
Scanner = MainInput.Scanner;
ScannerSoftware = MainInput.ScannerSoftware;
SequenceType = MainInput.SequenceType;
ReconType = MainInput.ReconType;
denoiseXe = MainInput.denoiseXe;
Analyst = MainInput.Analyst;
VoxelSize = MainInput.VoxelSize;
SliceOrientation = MainInput.SliceOrientation;
vent_file = MainInput.vent_file;
anat_file = MainInput.anat_file;
analysisFolder = MainInput.analysisFolder;
mask_file_name = MainInput.mask_file_name;
MainInput.OutputPath = analysisFolder;

xedatapath = vent_file;
Hdatapath = anat_file;

Diffusion = '';
GasExchange = '';
MainInput.AnalysisType = 'Ventilation';
MainInput.Institute = 'CCHMC'; 
MainInput.analysissesFolder = analysisFolder;
MainInput.CCHMC_DbVentAnalysis = 'yes';

Outputs = [];
Outputs.PatientID = SubjectID;
Outputs.Age = Age;
Outputs.Sex = Sex;
Outputs.Disease = Disease;
Outputs.ScanDate = ScanDate;
Outputs.XeImagepath = vent_file;
Outputs.HImagepath = anat_file;
Outputs.analysispath = analysisFolder;
Outputs.SequenceType = SequenceType;
Outputs.ReconType = ReconType;
Outputs.MaskPath = mask_file_name; 

Outputs.AnalysisCode_path = 'https://github.com/aboodbdaiwi/XIPline';
Outputs.AnalysisDate = str2double(datestr(datetime('today'), 'yyyymmdd'));

% Check if we have anatomical images or not
if isempty(Hdatapath)
    MainInput.NoProtonImage = 'yes';  % There is no proton images
else
    MainInput.NoProtonImage = 'no';    % There is  proton images 
end

analysisSubfolder = fullfile(analysisFolder, 'Ventilation_Analysis');
if ~exist(analysisSubfolder, 'dir')
    mkdir(analysisSubfolder);
end

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
Outputs.XeDataLocation = XeDataLocation;
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
    Outputs.HDataLocation = HDataLocation;
    Outputs.HFileName = HFileName;
    Outputs.XeDataext = H_ext;
end
MainInput.denoiseXe = 'no';
[Ventilation, ~, ~, Proton, MainInput] = LoadData.LoadReadData(MainInput);
Ventilation.Image = double(Ventilation.Image);
Ventilation.outputpath = analysisSubfolder;
Outputs.AnalysisCode_hash = MainInput.AnalysisCode_hash; 
% Ventilation.Image = Ventilation.Image(:,:,1:2:end);
% figure; Global.imslice(Ventilation.Image,'Ventilation')

if strcmp(MainInput.XeDataext, '.nii') || strcmp(MainInput.XeDataext, '.gz') 
    A = Ventilation.Image;
    if (size(A,1) < size(A,3)) 
        A = permute(A,[3 2 1]); 
        A = flipud(A);
        % A = imrotate(A,-90*2);
        % A = flip(A,1);  
        % A = flip(A,2);
        % % A = flip(A,3);
        Ventilation.Image = A;    
        % figure; Global.imslice(A,'Ventilation')
    elseif (size(A,2) < size(A,3))
        A = permute(A,[1 3 2]); 
        A = imrotate(A,-90);
        % A = flip(A,1);  
        A = flip(A,2);
        % % A = flip(A,3);
        Ventilation.Image = A;    
        % figure; Global.imslice(A,'Ventilation')        
    else
        A = imrotate(A,-90*2);
        % A = flip(A,1);  
        % A = flip(A,2);
        % A = flip(A,3);
        Ventilation.Image = A; 
    end
end
% figure; Global.imslice(Ventilation.Image,'Ventilation')

if strcmp(MainInput.NoProtonImage, 'no')
    if strcmp(MainInput.HDataext, '.nii') || strcmp(MainInput.HDataext, '.gz') 
        A = Proton.Image;
        A = permute(A,[2 3 1]); 
        A = imrotate(A,-90);
        A = flip(A,1);  
        A = flip(A,2);
        % A = flip(A,3);
        Proton.Image = A;  
    end
end
% figure; Global.imslice(A,'Ventilation')

  
% figure; Global.imslice(Proton.Image,'Proton')
% -----------------------------registration-----------------------------
MainInput.SkipRegistration = 0;
if strcmp(MainInput.NoProtonImage, 'no') 
    try
        MainInput.RegistrationType = 'ANTs'; 
        MainInput.TransformType = 'rigid'; % 'translation' | 'rigid' | 'similarity' | 'affine'
        [Proton] = Registration.PerformRegistration(Proton,Ventilation,GasExchange,MainInput);        
    catch
        MainInput.RegistrationType = 'Multimodal'; 
        MainInput.TransformType = 'affine'; % 'translation' | 'rigid' | 'similarity' | 'affine'
        % check if xenon and proton have the same number of slices
        if size(Ventilation.Image,3) == size(Proton.Image,3)
            MainInput.SliceSelection = 0;
        else
            MainInput.SliceSelection = 1;
            MainInput.Xestart = 1;
            MainInput.Xeend = size(Ventilation.Image,3);
            MainInput.Hstart = 1;
            MainInput.Hend = size(Proton.Image,3);
        end
    
        Xesize = size(Ventilation.Image);
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
    end
else
    Proton.ProtonRegistered = zeros(size(Ventilation.Image));
    Proton.ProtonRegisteredColored = zeros(size(Ventilation.Image));
end

% -----------------------------segmentation-----------------------------
% Check if mask_file_name is empty
if isempty(mask_file_name)
    SegmentMaskMode = 0; % 0 = new AI mask, 1 = load exisitng mask
    disp('The mask_file_name is empty.');
else
    SegmentMaskMode = 1; % 0 = new AI mask, 1 = load exisitng mask
    disp('The mask_file_name is not empty.');
end

%SegmentMaskMode = 1; % 0 = new AI mask, 1 = load exisitng mask
if SegmentMaskMode == 0
    cd(MainInput.XeDataLocation)
    
    % diary Log.txt
    MainInput.SegmentationMethod = 'Auto'; % 'Threshold' || 'Manual' || 'Auto'
    MainInput.SegmentAnatomy = 'Parenchyma'; % 'Airway'; || 'Parenchyma'
    MainInput.Imagestosegment = 'Xenon';  % 'Xe & Proton Registered' | 'Xenon' | 'Registered Proton'
    
    MainInput.thresholdlevel = 1; % 'threshold' 
    MainInput.SE = 1;
    
    MainInput.SegmentManual = 'Freehand'; % 'AppSegmenter' || 'Freehand'
    MainInput.SliceOrientation = SliceOrientation; % 'coronal' ||'transversal' || 'sagittal' ||'isotropic'
    [Proton,Ventilation,~,~] = Segmentation.PerformSegmentation(Proton,Ventilation,Diffusion,GasExchange,MainInput);

    if exist('Ventilation.AirwayMask', 'var')
        % skip
    else
        Ventilation.AirwayMask = zeros(size(Ventilation.LungMask));
    end
else
    % % List all files in the analysispath folder
    % files = dir(fullfile(analysispath, '*.gz')); %.nii
    
    % % Extract the file names
    % fileNames = {files.name};
    % maskpath = fullfile(analysispath,fileNames{end});
    % figure; Global.imslice(Ventilation.Image,'Image')
    % figure; Global.imslice(Ventilation.Image,'Ventilation.Image')
    testmask = Segmentation.SegmentLungthresh(Ventilation.Image,3,1); 
    structuringElementSize = 2;
    se = strel('disk', structuringElementSize); % Create a disk-shaped structuring element
    testmask = double(imclose(testmask, se)); % Perform morphological closing
    structuringElementSize = 2; % Size of the structuring element
    % Create a structuring element for erosion
    se = strel('cube', structuringElementSize); % Create a cube-shaped structuring element
    % Perform erosion to shrink the binary mask
    testmask = imerode(testmask, se); % Perform erosion
    % figure; Global.imslice(testmask,'testmask')
    
    % maskpath = fullfile(analysispath,mask_file_name);
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
    % figure; Global.imslice(A,'A')
    B = A;
    sizeArray1 = size(Ventilation.Image);
    sizeArray2 = size(A);
    
    % Check if slices in the 3d dimension 
    if ~isequal(sizeArray1, sizeArray2)
        minIndex = find(sizeArray2 == min(sizeArray2));
        if minIndex == 1
            B = permute(B, [2, 3, minIndex]);
        elseif minIndex == 2
            B = permute(B, [1, 3, minIndex]);
        elseif minIndex == 3
            B = permute(B, [1, 2, minIndex]);
        end
    else
        % Display a message indicating that the sizes match
        disp('Array sizes are the same.');
    end
    
    % check if number of mask slices don't match number of vent image slices
    if size(B,3) ~= size(Ventilation.Image,3)
        tempmask = zeros(size(Ventilation.Image));
        sliceIndex = zeros(1,size(testmask, 3));
        for slice = 1:size(testmask, 3)
            % Check if any pixel has a value of 1 in the current slice
            if any(testmask(:, :, slice) == 1, 'all')
                sliceIndex(slice) = 1;
            end
        end
        startslice = find(sliceIndex);
        if length(startslice) == size(B,3)
            tempmask(:,:,startslice(1):startslice(end)) = B;
        else
            try
                tempmask(:,:,startslice(1):startslice(1)+size(B,3)-1) = B;
            catch
                disp('mask has to be fixed manually')
            end
        end
        B = tempmask;
    end
    %figure; Global.imslice(B,'B')
    % check if orientation between mask and vent image is matching
    [DiceScore] = Global.DiceCoeff(testmask, B);
    if DiceScore.DSC < 0.7 % play with this thershold if it's still not working
        DSC = zeros(3,4,3,2);
        for op = 1:3
            if op == 1
                for i = 1:4
                    C = imrotate(B,90*i);
                    [DiceScore] = Global.DiceCoeff(testmask, C);
                    DSC(op,i,1,1) = DiceScore.DSC;
                end
            elseif op == 2
                for i = 1:4
                    C = imrotate(B,90*i);
                    for j = 1:3
                        C = flip(C,j);
                        [DiceScore] = Global.DiceCoeff(testmask, C);
                        DSC(op,i,j,1) = DiceScore.DSC;
                    end
                end
            elseif op == 3
                for i = 1:4
                    C = imrotate(B,90*i);
                    for j = 1:3
                        C = flip(C,j);
                        for k = 1:2
                            C = flip(C,k+1);
                            [DiceScore] = Global.DiceCoeff(testmask, C);
                            DSC(op,i,j,k) = DiceScore.DSC;
                        end
                    end
                end
            end
    
        end
        [~, linearIndex] = max(DSC(:));
        [indx1, indx2, indx3, indx4] = ind2sub(size(DSC), linearIndex);
        if indx1 == 1
            D = imrotate(B,90*indx2); 
        elseif indx1 == 2
            D = imrotate(B,90*indx2); 
            D = flip(D,indx3);
        elseif indx1 == 3
            D = imrotate(B,90*indx2); 
            D = flip(D,indx3);
            D = flip(D,indx4);
        end
    %     figure; Global.imslice(D,'mask')
        B = D;
    end
    Ventilation.LungMask = B;
    Ventilation.AirwayMask = zeros(size(Ventilation.LungMask));
    % figure; Global.imslice(D,'Image')
    % figure; Global.imslice(Ventilation.Image,'Image')
    % % figure; Global.imslice(Proton.Image,'proton')
    % figure; Global.imslice(Ventilation.LungMask,'mask')

end
% B = A;
% B = flip(B,1);
% % % % B = rot90(B);
% B = flip(B,2);
% B = flip(B,3);
% 
% % Define the target size for the new array (224x224x10)
% targetSize = [224 224 size(B, 3)];
% 
% % Create a new array of zeros with the target size
% A_centered = zeros(targetSize);
% 
% % Calculate the starting index to place A in the center of A_centered
% startX = floor((targetSize(1) - size(B, 1)) / 2) + 1;
% startY = floor((targetSize(2) - size(B, 2)) / 2) + 1;
% 
% % Place the original array A into the center of A_centered
% A_centered(startX:startX+size(B,1)-1, startY:startY+size(B,2)-1, :) = B;
% B = double(A_centered > 0);
% shifted_mask = zeros(size(B));
% % 
% % Shift the mask 2 pixels to the left for all slices
% shifted_mask(11:end, 1:end-15, :) = B(1:end-10, 16:end, :);
% % shifted_mask(11:end, 16:end, :) = B(1:end-10, 1:end-15, :);
% Ventilation.LungMask = B;
% figure; Global.imslice(Ventilation.LungMask,'mask')
% figure; Global.imslice(B,'mask')
% figure; Global.imslice(Ventilation.Image,'Image')

% Ventilation.LungMask(:,95:end,14) = testmask(:,95:end,14);
% Ventilation.LungMask = double(Ventilation.LungMask > 0);

if ischar(VoxelSize)
    VoxelSize = eval(VoxelSize);
end
Ventilation.DDI3Dx = VoxelSize(1);
Ventilation.DDI3Dy = VoxelSize(2);
Ventilation.DDI3Dz = VoxelSize(3);

Ventilation.ImageResolution = [VoxelSize(1), VoxelSize(2), VoxelSize(2)];
MainInput.PixelSpacing = Ventilation.ImageResolution(1:2);
MainInput.SliceThickness = Ventilation.ImageResolution(3);

MainInput.SegmentVessels = 0;
MainInput.vesselImageMode = 'xenon'; % xenon || proton
switch MainInput.vesselImageMode
    case 'xenon'
        MainInput.frangi_thresh = 0.25;
    case 'proton'  
        MainInput.frangi_thresh = 0.2; % you can change this threshold to increase or decrease the mask
end

if MainInput.SegmentVessels == 1
    [Ventilation] = Segmentation.Vasculature_filter(Proton, Ventilation, MainInput);
else
    Ventilation.VesselMask = zeros(size(Ventilation.Image));
    Ventilation.vessel_stack = zeros(size(Ventilation.Image));
end
maskarray = double(Ventilation.LungMask + Ventilation.VesselMask);
maskarray(maskarray > 1) = 0;
Ventilation.LungMask = double(maskarray);

% ----------------------------- VDP -----------------------------
Ventilation.cincibiasfieldcorr = 0;
Ventilation.N4Analysis = 0;
Ventilation.IncompleteThresh = 60;
Ventilation.RFCorrect = 0;
Ventilation.CompleteThresh = round(Ventilation.IncompleteThresh/2); % change from 15% to half of the IncompleteThresh
Ventilation.HyperventilatedThresh = 200;
Ventilation.HeterogeneityIndex = 'yes';
Ventilation.ThreshAnalysis = 'yes'; % 'yes'; || 'no'
Ventilation.LB_Analysis = 'no'; % 'yes'; || 'no'
Ventilation.LB_Normalization = 'LBmean'; % 'LBmean'; || 'GLBmean' || 'GLBpercentile'            
Ventilation.Kmeans = 'no';  % 'yes'; || 'no'
Ventilation.AKmeans = 'no';  % 'yes'; || 'no'
Ventilation.DDI2D = 'yes';  % 'yes'; || 'no'
Ventilation.DDI3D = 'no';  % 'yes'; || 'no'
Ventilation.DDIDefectMap = 'Threshold';
Ventilation.GLRLM_Analysis = 'yes'; % 'yes'; || 'no'
Ventilation.GLRLMDefectMap = 'Threshold';

Proton.AnatImage_Reg = Proton.ProtonRegistered;
Proton.AnatImage_RegColored = Proton.ProtonRegisteredColored;

close all;
Ventilation.N4Analysis = 1;
switch Ventilation.LB_Normalization
    case 'LBmean'
        Ventilation.LBThresholds = [0.33, 0.66, 1, 1.33, 1.66];  % mean
        Ventilation.Hdist = [0.001029, 0.995768, 0.24409];
    case 'GLBmean'
        Ventilation.LBThresholds = [0.479751, 0.739117, 0.985809, 1.223886, 1.455482];  % median
        Ventilation.Hdist = [0.000464, 0.981203, 0.243749]; 
    case 'GLBpercentile'
        Ventilation.LBThresholds = [0.28893, 0.462393, 0.622368, 0.77374, 0.918873]; % percentile
        Ventilation.Hdist = [-1.186698, 0.738826, 0.198451]; 
end   
MainInput.MaskFullPath = mask_file_name;

[Ventilation] = VentilationFunctions.Ventilation_Analysis(Ventilation, Proton, MainInput);
close all;
Outputs.timestamp = str2double(Ventilation.timestamp);
Outputs.SNR_slice = Ventilation.SNR_slice;
Outputs.Overall_SNR = Ventilation.SNR_lung;
Outputs.SNRvv_slice = Ventilation.SNRvv_slice;
Outputs.SNR_vv = Ventilation.SNR_vv;

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

OutputJSONFile = fullfile(analysisSubfolder, 'Ventilation_Analysis.json');
Global.exportStructToJSON(Outputs, OutputJSONFile);

% Check if the directory specified by analysisSubfolder exists, create it if necessary
if ~exist(analysisSubfolder, 'dir')
    mkdir(analysisSubfolder);
end

% Save Outputs to the MAT-file
save(fullfile(analysisSubfolder, 'workspace.mat'));
save(fullfile(analysisSubfolder, 'Ventilation_Analysis_Outputs.mat'), 'Outputs');
clearvars
% Global.tts('Ventilation Analysis completed');
end

