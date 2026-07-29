
clc;clear;

maskfile = '\\rds6.cchmc.org\PulMed-43\CPIR_Share\STiM (TomoStat)\ExampleDataSets\HyPOINT_Pre-Post_ETI\Healthy Patients\IRC168H-251\V1\lungmask_202508030120.nii.gz';
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
sliceOrientation = 'axial'; % coronal ('A-P' or 'P-A') or axial ('S-I' or 'I-S')
sliceDirection = 'I-S';
% sliceOrientation = 'coronal'; % coronal ('A-P' or 'P-A') or axial ('S-I' or 'I-S')
% sliceDirection = 'A-P';
lobar_mask = Segmentation.lobar_mask_registration_2d(input_2d_mask, sliceOrientation, sliceDirection);
lobar_mask = double(VentilationFunctions.medFilter_integer(lobar_mask,3)).*input_2d_mask;
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
