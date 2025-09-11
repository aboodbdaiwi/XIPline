function [separated_lung, lung] = separate_connected_lung(lung_binary2D)
% Separate two connected lungs by iterative erosion and restore by dilation.
% Returns:
%   separated_lung : logical mask of the separated lungs (right|left)
%   lung.right     : logical mask for the right component
%   lung.left      : logical mask for the left component
%
% Robust to cases where initial erosion deletes all pixels (nCC==0)
% or where the original mask already has >=2 components.

    lung_binary2D = logical(lung_binary2D);
    [nRows, nCols] = size(lung_binary2D);

    lung.right = false(nRows,nCols);
    lung.left  = false(nRows,nCols);
    separated_lung = lung_binary2D;

    if ~any(lung_binary2D(:))
        return; % nothing to do
    end

    % Connectivity & structuring elements
    conn = 8;
    se1  = strel('disk',1);

    % First erosion step
    eroded_lung = imerode(lung_binary2D, se1);
    nErosion = 1;
    [L, nCC] = bwlabel(eroded_lung, conn);

    % If erosion wiped everything, fall back to original labeling
    if nCC == 0
        [L0, n0] = bwlabel(lung_binary2D, conn);
        if n0 < 2
            % Cannot split reliably; return original as 'right', leave 'left' empty
            lung.right = lung_binary2D;
            separated_lung = lung.right;
            return;
        else
            % Already split in the original; use that directly
            L     = L0;
            nCC   = n0;
            nErosion = 0;
        end
    else
        % Keep eroding until we get 2 or more components (or until nothing left)
        while nCC == 1
            eroded_lung = imerode(eroded_lung, se1);
            nErosion = nErosion + 1;
            [L, nCC] = bwlabel(eroded_lung, conn);
            if nCC == 0
                % Gave up all pixels; revert to original labeling
                [L0, n0] = bwlabel(lung_binary2D, conn);
                if n0 < 2
                    lung.right = lung_binary2D;
                    separated_lung = lung.right;
                    return;
                else
                    L     = L0;
                    nCC   = n0;
                    nErosion = 0;
                    break;
                end
            end
        end
    end

    % Now we have nCC >= 2. Take the two largest components on the eroded map.
    stats = regionprops(L, 'Area');
    areas = [stats.Area]';
    [~, order] = sort(areas, 'descend');
    main   = order(1);
    second = order(2);

    rSeed = (L == main);
    lSeed = (L == second);

    % Dilate back the same number of steps; confine to original mask
    if nErosion > 0
        seN = strel('disk', nErosion);
        rDil = imdilate(rSeed, seN);
        lDil = imdilate(lSeed, seN);
    else
        rDil = rSeed;
        lDil = lSeed;
    end

    right = rDil & lung_binary2D;
    left  = lDil & lung_binary2D;

    % Resolve any overlap by assigning voxel to the nearer seed (distance maps)
    overlap = right & left;
    if any(overlap(:))
        Dr = bwdist(~rSeed);
        Dl = bwdist(~lSeed);
        assignToRight = Dr >= Dl; % ties to right
        left(overlap & assignToRight)  = false;
        right(overlap & ~assignToRight) = false;
    end

    lung.right = right;
    lung.left  = left;
    separated_lung = right | left;
end




% function [separated_lung,lung] = separate_connected_lung(lung_binary2D)
% %% Separate connected lung by # of erosion and then restore the
% % original size by the same # of dilation.
% % W. Zha @ 7/25/2016
% [nRows, nCols] = size(lung_binary2D);
% se1 = strel('disk',1);
% eroded_lung = imerode(lung_binary2D,se1);
% nErosion = 1;
% [labeled_lung, nMid] = bwlabel(eroded_lung);
% while nMid == 1 % erode by 1 pixel until separated
%     eroded_lung = imerode(eroded_lung,se1);
%     nErosion = nErosion + 1;
%     [labeled_lung, nMid] = bwlabel(eroded_lung);
%     if sum(labeled_lung(:))==0
% 
%         separated_lung = lung_binary2D;
%         lung.right = lung_binary2D;
%         lung.left = lung_binary2D;
%         return;
%     end
% end
% pixels = zeros(nMid,1);
% 
% for np = 1:nMid
%     pixels(np) = sum(labeled_lung(:) == np);
% end
% [~, sorted_ind] = sort(pixels, 'descend');
% 
% lung.right = zeros(nRows,nCols);
% lung.right(labeled_lung==sorted_ind(1)) = 1;
% seN = strel('disk',nErosion);
% lung.right = imdilate(lung.right,seN);
% lung.left = zeros(nRows,nCols);
% lung.left(labeled_lung==sorted_ind(2)) = 1;
% lung.left = imdilate(lung.left,seN);
% se2 = strel('disk',2);
% lung.left(imdilate(lung.right,se2)==1) = 0;
% separated_lung = lung.right + lung.left;
