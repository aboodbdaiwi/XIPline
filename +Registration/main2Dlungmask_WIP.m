
%2D multiresolution image registration usign level set and B-spline compose
%C-shape example
close all
clear all
clc;
tic
% addpaths


% fixedimgpath = 'D:\ML\CT_lobe_segmentation\LLS0003\LLS0003\Segmentations\LLS0003_LungLobes.nii.gz';
% fixedimgpath = 'Y:\Riaz\CCHMC_CT_data\IRC740H-025\mask.nii';
fixedimgpath = 'Y:\Riaz\CCHMC_CT_data\IRC740H-025\images.nii';

fixedimg = LoadData.load_nii(fixedimgpath);
fixedimg = fixedimg.img;
A = imrotate(fixedimg,90);
A = flip(A,2);
A = double(round(imresize3(A, [512, 512, 512], 'nearest')));
A = double((A - min(A(:)))/(max(A(:)) - min(A(:))));
figure; imslice(A);


% Ensure integer output
mask_resized = uint8(mask_resized);

% movingimgpath = 'D:\ML\CT_lobe_segmentation\LLS0004\LLS0004\Segmentations\LLS0004_LungLobes.nii.gz';
% movingimgpath = 'D:\OneDrive - cchmc\Lab\FLORET_vs_SoS_IRC740H\CF_IRC740_UTE\FLORET\img78\mask2.nii.gz';
% movingimgpath = 'X:\IRC740H_CF_Non-CF_Bronchiectasis\IRC740H-025\20230911\Vent Analysis vfinal\img_ventilation_mask.nii.gz';
movingimgpath = 'D:\OneDrive - cchmc\Lab\FLORET_vs_SoS_IRC740H\CF_IRC740_UTE\FLORET\img78\image.nii.gz';

movingimg = LoadData.load_nii(movingimgpath);
movingimg = movingimg.img;
B = imrotate(movingimg,90);
B = double(flip(B,2));
% B = B == 1;
B = double(B);
B = double((B - min(B(:)))/(max(B(:)) - min(B(:))));

% BB = double(flip(B,1));
% BB = double(flip(BB,2));
% analysisFolder = 'D:\OneDrive - cchmc\Lab\FLORET_vs_SoS_IRC740H\CF_IRC740_UTE\FLORET\img78';
% niftiwrite(fliplr(rot90(BB,-1)),[analysisFolder,'\mask2.nii'],'Compressed',true);

figure; imslice(B);

Aslice = 259;
Bslice = 240; % 8 for xenon image and 240 for FLORET

I0 = double(squeeze(A(Aslice,:,:))); 
% I1 = double(squeeze(B(Bslice,:,:)) > 0); 
% I1 = double(round(imresize(I1, [size(I0,1), size(I0,2)], 'nearest')));

I1 = double(squeeze(B(Bslice,:,1:370))); 
I1 = double(imresize(I1, [size(I0,1), size(I0,2)], 'nearest'));

figure; subplot(1,2,1);
imagesc(I0); colormap("gray"); title('Moving Image'); axis tight;
subplot(1,2,2);
imagesc(I1); colormap("gray"); title('Fixed Image'); axis tight;


filename='lungmask.avi';

nlevel=3;
max_Iteration=1000;
TOL=0.05;
gamma=1.54;
PlotAll = 1;
SaveVideo = 0;
UseGaussian=0;

% [I0] = MultiresolutionRegistration2D( I0,I1,nlevel,max_Iteration,PlotAll,SaveVideo,filename,UseGaussian,TOL,gamma);

%% 2D

clc; close all; clear all;

fixedimgpath = 'D:\ML\CT_lobe_segmentation\LLS0003\LLS0003\Segmentations\LLS0003_LungLobes.nii.gz';
fixedimg = LoadData.load_nii(fixedimgpath);
fixedimg = fixedimg.img;
A = imrotate(fixedimg,90);
A = flip(A,2);
% figure; imslice(A);

movingimgpath = 'D:\ML\CT_lobe_segmentation\LLS0004\LLS0004\Segmentations\LLS0004_LungLobes.nii.gz';
movingimg = LoadData.load_nii(movingimgpath);
movingimg = movingimg.img;
B = imrotate(movingimg,90);
B = flip(B,2);
% figure; imslice(B);

Aslice = 131;
Bslice = 129;
I1 = double(A(:,:,Aslice) > 0); 
I0 = double(B(:,:,Bslice)); 
%% 

% Convert to single precision for demons
fixed = single(I1);
moving = single(I0);
moving = imhistmatch(moving,fixed);
% --- Step 2: Deformable registration ---
[dispField, moving_reg] = imregdemons(...
    moving, fixed, ...
    [500 400 300 200 100], ...
    'AccumulatedFieldSmoothing', 0.5, ...
    'PyramidLevels', 5);

% figure; imslice(moving_reg);
% --- Step 3: Warp and relabel lobes ---
registered_lobar = zeros(size(fixed), 'like', I0);

for classID = 1:5
    binary_mask = I0 == classID;
    warped_mask = imwarp(binary_mask, dispField, 'Interp', 'nearest');

    if ~isequal(size(warped_mask), size(fixed))
        warped_mask = imresize(warped_mask, size(fixed), 'nearest');
    end

    % Keep inside lung mask
    warped_mask = warped_mask & logical(fixed);

    registered_lobar(warped_mask > 0) = classID;
end

% --- Step 4: Dilate each class (2D, per slice) ---
se = strel('disk', 5);  % use 2D disk for slice
dilated_lobar = registered_lobar;

for classID = 1:5
    class_mask = registered_lobar == classID;

    dilated_mask = imdilate(class_mask, se);
    dilated_mask = dilated_mask & (dilated_lobar == 0) & (fixed > 0);

    dilated_lobar(dilated_mask) = classID;
end

% Combine with original
registered_lobar = max(registered_lobar, dilated_lobar);

% --- Step 5: Fill holes (2D) ---
for classID = 1:5
    mask = registered_lobar == classID;
    mask_filled = imfill(mask, 'holes');
    registered_lobar(mask_filled & ~mask) = classID;
end

% --- Step 6: Visualize ---
figure;
subplot(1,3,1); imshow(I0, []); title('Original Lobar Mask (Moving)');
subplot(1,3,2); imshow(I1, []); title('Fixed Lung Mask');
subplot(1,3,3); imshow(registered_lobar, []); title('Warped + Dilated Mask');

%% 3D
clc;

% Step 0: Inputs
I1 = double(B(:,:,:) > 0);     % Binary fixed lung mask
I0 = double(A(:,:,:));         % Moving lobar segmentation (1–5)

% Step 1: Compute 3D deformable registration
moving = single(I0);
fixed  = single(I1);

[dispField, ~] = imregdemons(...
    moving, fixed, ...
    [500 400 300 200 100], ...
    'AccumulatedFieldSmoothing', 0.5, ...
    'PyramidLevels', 5);

% Step 2: Warp each lobar label using displacement field
registered_lobar = zeros(size(fixed), 'like', I0);

for classID = 1:5
    binary_mask = I0 == classID;
    warped_mask = imwarp(binary_mask, dispField, 'Interp', 'nearest');

    if ~isequal(size(warped_mask), size(fixed))
        warped_mask = imresize3(warped_mask, size(fixed), 'nearest');
    end

    % Ensure no overlap during assignment
    registered_lobar(warped_mask > 0) = classID;
end

% Step 3: Dilation for each lobe
se = strel('sphere', 5);  % Dilation radius: 2 voxels
dilated_lobar = registered_lobar;

for classID = 1:5
    class_mask = registered_lobar == classID;

    % Dilate the class
    dilated_mask = imdilate(class_mask, se);

    % Keep within lung mask and avoid label overlap
    dilated_mask = dilated_mask & (dilated_lobar == 0) & (I1 > 0);

    % Assign class label
    dilated_lobar(dilated_mask) = classID;
end

% Step 4: Final result
registered_lobar = max(registered_lobar, dilated_lobar);

% Step 5: Visualize middle axial slice
midSlice = round(size(fixed, 3) / 2);

figure;
subplot(1,3,1); imshow(I0(:,:,midSlice), []); title('Original Lobar (Moving)');
subplot(1,3,2); imshow(I1(:,:,midSlice), []); title('Fixed Lung Mask');
subplot(1,3,3); imshow(registered_lobar(:,:,midSlice), []); title('Warped + Dilated Lobes');


figure; imslice(registered_lobar)

%% 

% mask2 = logical(I0);  % Ensure it's binary
I02 = centerMultiClassMask(I0);
I12 = centerMultiClassMask(I1);

% imslice(centered);



function centered_mask = centerMultiClassMask(mask)
    % mask: 2D multi-class mask (integer values 0–5)

    % Get all non-zero indices (ignore background)
    [rowIdx, colIdx] = find(mask > 0);
    if isempty(rowIdx)
        error('The input mask contains only background.');
    end

    % Bounding box of all non-zero regions (any lobe)
    top = min(rowIdx);
    bottom = max(rowIdx);
    left = min(colIdx);
    right = max(colIdx);

    % Crop mask to bounding box
    cropped = mask(top:bottom, left:right);
    [h, w] = size(cropped);

    % Original mask size
    [rows, cols] = size(mask);

    % Initialize centered mask
    centered_mask = zeros(rows, cols, 'like', mask);

    % Compute start positions for centering
    row_start = floor((rows - h)/2) + 1;
    col_start = floor((cols - w)/2) + 1;

    % Paste cropped region into center of original-size mask
    centered_mask(row_start:row_start+h-1, col_start:col_start+w-1) = cropped;
end


