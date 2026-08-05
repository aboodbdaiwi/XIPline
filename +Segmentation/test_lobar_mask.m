
clc;clear;

maskfile = '\\rds6.cchmc.org\PulMed-43\CPIR_Share\STiM (TomoStat)\ExampleDataSets\HyPOINT_Pre-Post_ETI\Needs adjusting\IRC186H-1029\V2\lungmask_202509082242.nii.gz';
[maskfolder, filename] = fileparts(maskfile);
try
    A1 = LoadData.load_nii(maskfile);
catch
    A1 = LoadData.load_untouch_nii(maskfile);
end

A = double(squeeze(A1.img));
% A = permute(A, [1, 3, 2]);
% Match your original 3D atlas orientation convention
A = imrotate(A, -90);
A = flip(A, 2);
% A = flip(A, 3);
lungmask = A;

% imslice(A);

% 2D

input_2d_mask = lungmask;
sliceOrientation = 'coronal';      % 'axial' or 'coronal'

if strcmpi(sliceOrientation, 'axial')
    sliceDirection = 'I-S';      % or 'S-I'
elseif strcmpi(sliceOrientation, 'coronal')
    sliceDirection = 'A-P';      % or 'P-A'
else
    error('Invalid slice orientation. Use ''axial'' or ''coronal''.');
end

lobar_mask = Segmentation.lobar_mask_registration_2d(input_2d_mask, sliceOrientation, sliceDirection);

% imslice(lobar_mask);

niftiwrite(abs(fliplr(rot90(lobar_mask, -1))), fullfile(maskfolder, 'lobar_mask'), 'Compressed', true);


%% 3D

input_3d_mask = double(ProtonMaskRegistered);
lobar_mask = Segmentation.lobar_mask_registration_3d(input_3d_mask);
lobar_mask1 = lobar_mask;
lobar_mask1 = permute(lobar_mask1, [3,2,1]);
figure; imslice(lobar_mask1);
figure; montage(lobar_mask,'Size',[11 11],'DisplayRange',[0 5]);
% niftiwrite(abs(fliplr(rot90(lobar_mask, -1))), fullfile(maskfolder, 'lobar_mask'), 'Compressed', true);
%% Binary mask RGB stack
figure('Color','w');

subplot(2,1,1)
montage(logical(lungmask),'Size',[1 size(lungmask,3)]);
title('Binary Lung Mask')

subplot(2,1,2)
montage(uint8(lobar_mask),'Size',[1 size(lobar_mask,3)], ...
    'DisplayRange',[0 5]);
colormap(gca, jet(6))
colorbar
title('Lobar Segmentation')

%% Batch 2D processing 
clc;
clear;
close all;

mainFolder = '\\rds6.cchmc.org\PulMed-43\CPIR_Share\STiM (TomoStat)\ExampleDataSets\HyPOINT_Pre-Post_ETI\Healthy Patients';

subjectFolders = dir(mainFolder);
subjectFolders = subjectFolders([subjectFolders.isdir]);
subjectFolders = subjectFolders(~ismember({subjectFolders.name}, {'.', '..'}));

for subjectIndex = 2:numel(subjectFolders)

    subjectFolder = fullfile(mainFolder, subjectFolders(subjectIndex).name);

    visitFolders = dir(fullfile(subjectFolder, 'V*'));
    visitFolders = visitFolders([visitFolders.isdir]);

    for visitIndex = 1:numel(visitFolders)

        visitFolder = fullfile(subjectFolder, visitFolders(visitIndex).name);

        maskFiles = dir(fullfile(visitFolder, 'lungmask*.nii.gz'));

        if isempty(maskFiles)
            continue;
        end

        maskfile = fullfile(visitFolder, maskFiles(1).name);

        try
            A1 = LoadData.load_nii(maskfile);
        catch
            A1 = LoadData.load_untouch_nii(maskfile);
        end

        A = double(squeeze(A1.img));

        A = imrotate(A, -90);
        A = flip(A, 2);

        lungmask = A;

        input_2d_mask = lungmask;
        sliceOrientation = 'coronal';
        sliceDirection = 'A-P';

        lobar_mask = Segmentation.lobar_mask_registration_2d( ...
            input_2d_mask, ...
            sliceOrientation, ...
            sliceDirection);

        niftiwrite( ...
            abs(fliplr(rot90(lobar_mask, -1))), ...
            fullfile(visitFolder, 'lobar_mask'), ...
            'Compressed', true);
    end
end
%% run VDP



% Batch 2D processing 
clc;
clear;
close all;

mainFolder = '\\rds6.cchmc.org\PulMed-43\CPIR_Share\STiM (TomoStat)\ExampleDataSets\HyPOINT_Pre-Post_ETI\Healthy Patients';

subjectFolders = dir(mainFolder);
subjectFolders = subjectFolders([subjectFolders.isdir]);
subjectFolders = subjectFolders(~ismember({subjectFolders.name}, {'.', '..'}));

for subjectIndex = 21%:numel(subjectFolders)

    subjectFolder = fullfile(mainFolder, subjectFolders(subjectIndex).name);

    visitFolders = dir(fullfile(subjectFolder, 'V*'));
    visitFolders = visitFolders([visitFolders.isdir]);

    for visitIndex = 1:numel(visitFolders)

        visitFolder = fullfile(subjectFolder, visitFolders(visitIndex).name);

        imageFiles = dir(fullfile(visitFolder, 'image*.nii.gz'));

        maskFiles = dir(fullfile(visitFolder, 'lungmask*.nii.gz'));
        vesselsmaskFiles = dir(fullfile(visitFolder, 'VesselMask*.nii.gz'));
        
        imagefile = fullfile(visitFolder, imageFiles(1).name);
        try
            A1 = LoadData.load_nii(imagefile);
        catch
            A1 = LoadData.load_untouch_nii(imagefile);
        end
        A = double(squeeze(A1.img));
        A = imrotate(A, -90);
        A = flip(A, 2);
        Image = A; % figure; imslice(Image)

        maskfile = fullfile(visitFolder, maskFiles(1).name);
        try
            A1 = LoadData.load_nii(maskfile);
        catch
            A1 = LoadData.load_untouch_nii(maskfile);
        end
        A = double(squeeze(A1.img));
        A = imrotate(A, -90);
        A = flip(A, 2);
        lungmask = A;

        vesselsmaskfile = fullfile(visitFolder, vesselsmaskFiles(1).name);
        try
            A1 = LoadData.load_nii(vesselsmaskfile);
        catch
            A1 = LoadData.load_untouch_nii(vesselsmaskfile);
        end
        A = double(squeeze(A1.img));
        % A = imrotate(A, -90);
        % A = flip(A, 2);
        % vesselsmask = A;

        A = imrotate(A, -90);
        A = flip(A, 1);
        vesselsmask = A;
        % imslice(vesselsmask)

        maskarray = double(lungmask - vesselsmask).*lungmask;
        maskarray(maskarray > 0) = 1;
        LungMask = double(maskarray);

        Outputs = VentilationFunctions.analyzeVentilationImage2(Image, LungMask, visitFolder);

        % Ventilation.LungMaskOriginal = lungmask;
        % sliceOrientation = 'coronal';
        % sliceDirection = 'A-P';
        % Proton = [];
        % MainInput.vesselImageMode = 'xenon';
        % MainInput.frangi_thresh = 0.25;
        % MainInput.SliceOrientation = 'coronal';
        % MainInput.XeDataLocation = visitFolder;
        % Ventilation = Segmentation.Vasculature_filter(Proton, Ventilation, MainInput);
        % VesselMask = Ventilation.VesselMask;
        % % VesselMaskfile = fullfile(visitFolder, 'VesselMask.mat');
        % % save(VesselMaskfile, 'VesselMask');
        % niftiwrite( ...
        %     abs(fliplr(rot90(VesselMask, -1))), ...
        %     fullfile(visitFolder, 'VesselMask.nii'), ...
        %     'Compressed', true);







    end
end

