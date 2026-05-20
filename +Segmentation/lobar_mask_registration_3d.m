function registered_mask = lobar_mask_registration_3d(input_3d_mask)
    
    % input_3d_mask = ProtonMaskRegistered;
    input_3d_mask = squeeze(input_3d_mask);
    input_3d_mask = logical(input_3d_mask);
    input_3d_mask = permute(input_3d_mask, [3,2,1]);
    % imslice(input_3d_mask)
    % script_location = 'D:\Github\XIPline\+Segmentation';
    script_location = fileparts(mfilename('fullpath'));
    
    lobar_atlas_file = fullfile(script_location, ...
        'LungAndLobeEstimation', ...
        'ProtonLobeMasks', ...
        '0006Lobes.nii.gz');
    
    try
        A1 = LoadData.load_nii(lobar_atlas_file);
    catch
        A1 = LoadData.load_untouch_nii(lobar_atlas_file);
    end
    
    A = double(squeeze(A1.img));
    A = imrotate(A,90);
    A = flip(A,3);
    A = flip(A,2);
    % imslice(A)
    
    lobar_mask = round(A);
    lobar_mask(lobar_mask < 0) = 0;
    lobar_mask(lobar_mask > 5) = 0;
    
    targetSize = size(input_3d_mask);
    lobar_mask = resizeLabelMask3D(lobar_mask,targetSize);
    
    input_lung = input_3d_mask > 0;
    % imslice(lobar_mask_shifted)
    
    lobar_mask = alignNonzeroSlices(lobar_mask,input_lung);
    lobar_mask = resizeLabelMask3D(lobar_mask,targetSize);
    
    atlas_lung = lobar_mask > 0;
    
    shiftVec = centerOfMass3D(input_lung) - centerOfMass3D(atlas_lung);
    
    lobar_mask_shifted = imtranslate3D_label(lobar_mask,shiftVec);
    atlas_lung_shifted = lobar_mask_shifted > 0;
    
    atlas_lung_shifted = imfill(atlas_lung_shifted,'holes');
    input_lung_fill = imfill(input_lung,'holes');
    
    try
        [D,~] = imregdemons(single(atlas_lung_shifted), ...
                            single(input_lung_fill), ...
                            [500 300 150], ...
                            'AccumulatedFieldSmoothing',2.0, ...
                            'PyramidLevels',3, ...
                            'DisplayWaitbar',false);
    catch
        [D,~] = imregdemons(single(atlas_lung_shifted), ...
                            single(input_lung_fill), ...
                            [300 150 75], ...
                            'AccumulatedFieldSmoothing',1.5, ...
                            'PyramidLevels',3);
    end
    
    registered_mask = zeros(targetSize);
    for label = 1:5
        this_lobe = lobar_mask_shifted == label;
        this_lobe_warped = imwarp(single(this_lobe), D);
        this_lobe_warped = this_lobe_warped > 0.5;
        this_lobe_warped = imresize3(this_lobe_warped, ...
            targetSize, 'nearest');
        registered_mask(this_lobe_warped > 0) = label;
    end
    
    registered_mask = round(registered_mask);
    registered_mask(~input_lung) = 0;
    
    empty_lung = input_lung & registered_mask == 0;
    if any(empty_lung(:))
        registered_mask = fillEmptyLungVoxelsByNearestLobe(registered_mask,input_lung);
    end
    
    registered_mask(~input_lung) = 0;
    registered_mask = uint8(registered_mask);
    imslice(registered_mask)
end

function out = resizeLabelMask3D(mask,targetSize)
    labels = unique(mask(:));
    labels(labels == 0) = [];
    out = zeros(targetSize);
    for i = 1:numel(labels)
        L = labels(i);
        tmp = mask == L;
        tmp2 = zeros(targetSize);
        for z = 1:targetSize(3)
            if size(tmp,3) == 1
                sliceIdx = 1;
            else
                sliceIdx = round((z-1) * (size(tmp,3)-1) / ...
                    (targetSize(3)-1)) + 1;
            end
            slice2D = tmp(:,:,sliceIdx);
            tmp2(:,:,z) = imresize(slice2D, ...
                targetSize(1:2), ...
                'nearest');
        end
        out(tmp2 > 0) = L;
    end
    out = round(out);
    out(out < 0) = 0;
    out(out > 5) = 0;
end


function out = alignNonzeroSlices(lobar_mask,input_lung)
    
    targetSize = size(input_lung);
    
    atlas_lung = lobar_mask > 0;
    
    input_z = squeeze(any(any(input_lung,1),2));
    atlas_z = squeeze(any(any(atlas_lung,1),2));
    
    input_slices = find(input_z);
    atlas_slices = find(atlas_z);
    
    if isempty(input_slices) || isempty(atlas_slices)
        out = resizeLabelMask3D(lobar_mask,targetSize);
        return
    end
    
    input_min = min(input_slices);
    input_max = max(input_slices);
    
    atlas_min = min(atlas_slices);
    atlas_max = max(atlas_slices);
    
    atlas_crop = lobar_mask(:,:,atlas_min:atlas_max);
    
    nInputSlices = input_max - input_min + 1;
    
    atlas_crop_resized = resizeLabelMask3D(atlas_crop, ...
        [targetSize(1), targetSize(2), nInputSlices]);
    
    out = zeros(targetSize);
    out(:,:,input_min:input_max) = atlas_crop_resized;

end


function c = centerOfMass3D(mask)
    
    idx = find(mask);
    
    if isempty(idx)
        c = [0 0 0];
        return
    end
    
    [y,x,z] = ind2sub(size(mask),idx);
    
    c = [mean(x), mean(y), mean(z)];

end


function out = imtranslate3D_label(mask,shiftVec)
    
    targetSize = size(mask);
    out = zeros(targetSize);
    
    for label = 1:5
        tmp = mask == label;
    
        tmpShift = imtranslate(tmp, ...
            [shiftVec(1), shiftVec(2), shiftVec(3)], ...
            'nearest', ...
            'FillValues',0);
    
        out(tmpShift > 0) = label;
    end
    
    out = round(out);

end


function registered_mask = fillEmptyLungVoxelsByNearestLobe(registered_mask,input_lung)
    
    registered_mask = double(registered_mask);
    
    empty_lung = input_lung & registered_mask == 0;
    known_lobes = registered_mask > 0;
    
    if ~any(empty_lung(:)) || ~any(known_lobes(:))
        return
    end
    
    [~,nearestIdx] = bwdist(known_lobes);
    
    registered_mask(empty_lung) = registered_mask(nearestIdx(empty_lung));
    
    registered_mask(~input_lung) = 0;

end