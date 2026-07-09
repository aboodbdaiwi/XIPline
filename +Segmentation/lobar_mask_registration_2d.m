function registered_mask = lobar_mask_registration_2d(input_2d_mask, sliceOrientation, sliceDirection)

    input_lung = logical(squeeze(input_2d_mask));

    if ndims(input_lung) == 2
        input_lung = reshape(input_lung, size(input_lung,1), size(input_lung,2), 1);
    end

    targetSize = size(input_lung);

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

    lobar_atlas = round(A);
    lobar_atlas(lobar_atlas < 0) = 0;
    lobar_atlas(lobar_atlas > 5) = 0;

    atlas_2d = extractAtlasSlices2D(lobar_atlas, sliceOrientation, sliceDirection);
    atlas_2d = keepNonzeroSlices(atlas_2d);

    lobar_mask = resizeLabelMask2DStack(atlas_2d, targetSize);
    lobar_mask = alignNonzeroSlices2D(lobar_mask, input_lung);

    [left_lung, right_lung] = splitLeftRightLungsGlobal(input_lung, lobar_mask);

    atlas_lung = lobar_mask > 0;

    if any(atlas_lung(:)) && any(input_lung(:))
        shiftVec = centerOfMass2DStack(input_lung) - centerOfMass2DStack(atlas_lung);
        lobar_mask = imtranslate2DStack_label(lobar_mask, shiftVec);
    end

    atlas_lung = lobar_mask > 0;

    registered_mask = zeros(targetSize);

    for s = 1:targetSize(3)

        fixed = imfill(input_lung(:,:,s), 'holes');
        moving = imfill(atlas_lung(:,:,s), 'holes');

        if ~any(fixed(:)) || ~any(moving(:))
            continue
        end

        try
            [D,~] = imregdemons(single(moving), ...
                                single(fixed), ...
                                [300 150 75], ...
                                'AccumulatedFieldSmoothing',1.5, ...
                                'PyramidLevels',3, ...
                                'DisplayWaitbar',false);
        catch
            [D,~] = imregdemons(single(moving), ...
                                single(fixed), ...
                                [150 75 50], ...
                                'AccumulatedFieldSmoothing',1.0, ...
                                'PyramidLevels',3);
        end

        for label = 1:5

            this_lobe = lobar_mask(:,:,s) == label;

            if ~any(this_lobe(:))
                continue
            end

            this_lobe_warped = imwarp(single(this_lobe), D) > 0.5;

            sliceOut = registered_mask(:,:,s);
            sliceOut(this_lobe_warped & sliceOut == 0) = label;
            registered_mask(:,:,s) = sliceOut;
        end
    end

    registered_mask = nearestFillExactMask2D(registered_mask, input_lung);

    registered_mask = enforceLeftRightLabelsNearest( ...
        registered_mask, input_lung, left_lung, right_lung);

    registered_mask = nearestFillExactMask2D(registered_mask, input_lung);

    registered_mask(~input_lung) = 0;

    mismatchPixels = nnz((registered_mask > 0) ~= input_lung);
    fprintf('Lobar-mask mismatch pixels = %d\n', mismatchPixels);

    assert(mismatchPixels == 0, ...
        'Final lobar mask does not exactly match input lung mask.');

    registered_mask = uint8(registered_mask);

end


function atlas_2d = extractAtlasSlices2D(lobar_atlas, sliceOrientation, sliceDirection)

    sliceOrientation = lower(strtrim(sliceOrientation));
    sliceDirection = upper(strrep(strtrim(sliceDirection),' ',''));

    switch sliceOrientation

        case 'coronal'
            atlas_2d = permute(lobar_atlas,[3 2 1]);

            switch sliceDirection
                case 'A-P'
                case 'P-A'
                    atlas_2d = flip(atlas_2d,3);
                otherwise
                    error('For coronal, use A-P or P-A.');
            end

        case 'axial'
            atlas_2d = lobar_atlas;

            switch sliceDirection
                case 'S-I'
                case 'I-S'
                    atlas_2d = flip(atlas_2d,3);
                otherwise
                    error('For axial, use S-I or I-S.');
            end

        otherwise
            error('sliceOrientation must be coronal or axial.');
    end

    atlas_2d = round(atlas_2d);
    atlas_2d(atlas_2d < 0) = 0;
    atlas_2d(atlas_2d > 5) = 0;
end


function out = keepNonzeroSlices(mask)

    z_has_lung = squeeze(any(any(mask > 0,1),2));
    idx = find(z_has_lung);

    if isempty(idx)
        out = mask;
    else
        out = mask(:,:,min(idx):max(idx));
    end
end


function out = resizeLabelMask2DStack(mask,targetSize)

    if numel(targetSize) == 2
        targetSize(3) = 1;
    end

    out = zeros(targetSize);

    nSourceSlices = size(mask,3);
    nTargetSlices = targetSize(3);

    for label = 1:5

        tmp = mask == label;
        tmpOut = false(targetSize);

        for s = 1:nTargetSlices

            if nTargetSlices == 1 || nSourceSlices == 1
                sourceIdx = round((nSourceSlices + 1) / 2);
            else
                sourceIdx = round((s-1) * (nSourceSlices-1) / ...
                                  (nTargetSlices-1)) + 1;
            end

            tmpOut(:,:,s) = imresize(tmp(:,:,sourceIdx), targetSize(1:2), 'nearest');
        end

        out(tmpOut) = label;
    end

    out = round(out);
    out(out < 0) = 0;
    out(out > 5) = 0;
end


function out = alignNonzeroSlices2D(lobar_mask,input_lung)

    targetSize = size(input_lung);

    input_z = squeeze(any(any(input_lung,1),2));
    atlas_z = squeeze(any(any(lobar_mask > 0,1),2));

    input_slices = find(input_z);
    atlas_slices = find(atlas_z);

    if isempty(input_slices) || isempty(atlas_slices)
        out = resizeLabelMask2DStack(lobar_mask,targetSize);
        return
    end

    input_min = min(input_slices);
    input_max = max(input_slices);

    atlas_min = min(atlas_slices);
    atlas_max = max(atlas_slices);

    atlas_crop = lobar_mask(:,:,atlas_min:atlas_max);

    nInputSlices = input_max - input_min + 1;

    atlas_crop_resized = resizeLabelMask2DStack(atlas_crop, ...
        [targetSize(1), targetSize(2), nInputSlices]);

    out = zeros(targetSize);
    out(:,:,input_min:input_max) = atlas_crop_resized;
end


function [left_lung,right_lung] = splitLeftRightLungsGlobal(input_lung,lobar_mask)

    left_lung = false(size(input_lung));
    right_lung = false(size(input_lung));

    atlas_left_x = meanX(lobar_mask == 1 | lobar_mask == 2);
    atlas_right_x = meanX(lobar_mask == 3 | lobar_mask == 4 | lobar_mask == 5);

    if isnan(atlas_left_x) || isnan(atlas_right_x)
        atlas_left_is_larger_x = true;
    else
        atlas_left_is_larger_x = atlas_left_x > atlas_right_x;
    end

    CC = bwconncomp(input_lung, 26);

    if CC.NumObjects >= 2

        areas = cellfun(@numel, CC.PixelIdxList);
        [~,ord] = sort(areas,'descend');

        compA = false(size(input_lung));
        compB = false(size(input_lung));

        compA(CC.PixelIdxList{ord(1)}) = true;
        compB(CC.PixelIdxList{ord(2)}) = true;

        xA = meanX(compA);
        xB = meanX(compB);

        if atlas_left_is_larger_x
            if xA > xB
                left_lung = compA;
                right_lung = compB;
            else
                left_lung = compB;
                right_lung = compA;
            end
        else
            if xA < xB
                left_lung = compA;
                right_lung = compB;
            else
                left_lung = compB;
                right_lung = compA;
            end
        end

        remaining = input_lung & ~(left_lung | right_lung);

        if any(remaining(:))
            distL = bwdist(left_lung);
            distR = bwdist(right_lung);

            left_lung(remaining & distL <= distR) = true;
            right_lung(remaining & distR < distL) = true;
        end

    else

        [~,x,~] = ind2sub(size(input_lung),find(input_lung));
        midx = mean(x);

        [Xgrid,~,~] = meshgrid(1:size(input_lung,2), ...
                               1:size(input_lung,1), ...
                               1:size(input_lung,3));

        if atlas_left_is_larger_x
            left_lung = input_lung & Xgrid >= midx;
            right_lung = input_lung & Xgrid < midx;
        else
            left_lung = input_lung & Xgrid < midx;
            right_lung = input_lung & Xgrid >= midx;
        end
    end

    left_lung = left_lung & input_lung;
    right_lung = right_lung & input_lung;
end


function registered_mask = nearestFillExactMask2D(registered_mask,input_lung)

    registered_mask = double(registered_mask);
    registered_mask(~input_lung) = 0;

    known = registered_mask > 0;
    missing = input_lung & registered_mask == 0;

    if any(missing(:)) && any(known(:))
        [~,idx] = bwdist(known);
        registered_mask(missing) = registered_mask(idx(missing));
    end

    missing = input_lung & registered_mask == 0;

    if any(missing(:))
        for s = 1:size(input_lung,3)
            slice = registered_mask(:,:,s);
            miss = missing(:,:,s);

            if any(miss(:))
                known2D = slice > 0;

                if any(known2D(:))
                    [~,idx2D] = bwdist(known2D);
                    slice(miss) = slice(idx2D(miss));
                else
                    slice(miss) = 1;
                end
            end

            registered_mask(:,:,s) = slice;
        end
    end

    registered_mask(~input_lung) = 0;
end


function registered_mask = enforceLeftRightLabelsNearest(registered_mask,input_lung,left_lung,right_lung)

    registered_mask = double(registered_mask);

    registered_mask(~input_lung) = 0;

    bad_left = left_lung & ismember(registered_mask,[3 4 5]);
    bad_right = right_lung & ismember(registered_mask,[1 2]);

    registered_mask(bad_left | bad_right) = 0;

    for s = 1:size(input_lung,3)

        slice = registered_mask(:,:,s);

        Lmask = left_lung(:,:,s);
        Rmask = right_lung(:,:,s);

        left_missing = Lmask & slice == 0;
        left_known = Lmask & ismember(slice,[1 2]);

        if any(left_missing(:)) && any(left_known(:))
            [~,idx] = bwdist(left_known);
            slice(left_missing) = slice(idx(left_missing));
        end

        right_missing = Rmask & slice == 0;
        right_known = Rmask & ismember(slice,[3 4 5]);

        if any(right_missing(:)) && any(right_known(:))
            [~,idx] = bwdist(right_known);
            slice(right_missing) = slice(idx(right_missing));
        end

        registered_mask(:,:,s) = slice;
    end

    registered_mask(~input_lung) = 0;
end


function c = centerOfMass2DStack(mask)

    idx = find(mask);

    if isempty(idx)
        c = [0 0];
        return
    end

    [y,x,~] = ind2sub(size(mask),idx);
    c = [mean(x), mean(y)];
end


function xmean = meanX(mask)

    idx = find(mask);

    if isempty(idx)
        xmean = NaN;
        return
    end

    [~,x,~] = ind2sub(size(mask),idx);
    xmean = mean(x);
end


function out = imtranslate2DStack_label(mask,shiftVec)

    out = zeros(size(mask));

    for s = 1:size(mask,3)

        for label = 1:5

            tmp = mask(:,:,s) == label;

            if ~any(tmp(:))
                continue
            end

            tmpShift = imtranslate(tmp, ...
                [shiftVec(1), shiftVec(2)], ...
                'nearest', ...
                'FillValues',0);

            sliceOut = out(:,:,s);
            sliceOut(tmpShift > 0 & sliceOut == 0) = label;
            out(:,:,s) = sliceOut;
        end
    end

    out = round(out);
end