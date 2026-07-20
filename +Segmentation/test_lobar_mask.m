
clc;clear;

maskfile = '\\rds6.chmccorp.cchmc.org\PulMed-54\CPIR_Images_Database\IRC1031\analysis\gx_v100\sub-0220\ses-20260602\ser-123238\GasExchange_Analysis\ProtonMaskRegistered.nii.gz';
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

imslice(A);

%% 2D

input_2d_mask = lungmask;
sliceOrientation = 'coronal';
sliceDirection = 'A-P';
lobar_mask = Segmentation.lobar_mask_registration_2d(input_2d_mask, sliceOrientation, sliceDirection);

imslice(lobar_mask);

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