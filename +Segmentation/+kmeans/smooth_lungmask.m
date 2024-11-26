function smoothed_mask = smooth_lungmask(lungmask,acquisition)
%% Smooth lung mask by cubic splines.
% Step 1: detect the number of isolated pieces;
% Step 2: separate connected lung if applicable;
% Step 3: morphometric closing to close small open areas in the mask;
% Step 4: convert piece mask into contours and then do a spline fit;
% Step 5: put all piece masks back into the smoothed version of the lung
%         mask;
% W. Zha@07/2016
%%

[nRows, nCols, nSlices] = size(lungmask);
smoothed_mask = zeros(nRows, nCols, nSlices);
contoured_slices = find(squeeze(sum(sum(lungmask,1),2))>0);

% lungmask = morphometric_closing(lungmask);
lungmask = fill_lungsmask(lungmask,'Coronal');
% Spline fit parameters
switch acquisition 
    case 'Coronal'
    spline_coef.breaks = 12;%32;%
    spline_coef.order = 4;%8;%
    case 'Axial'
    spline_coef.breaks = 4;
    spline_coef.order = 8;
    otherwise
         spline_coef.breaks = 20;
    spline_coef.order = 8;
end

for n = contoured_slices.'
    slice_mask = lungmask(:,:,n);
    % slice mask labeling and separation if needed
    [labeled_slice,num] = bwlabel(slice_mask);
    if num ==1
        separated_slice = separate_connected_lung(slice_mask);
        [labeled_slice,num] = bwlabel(separated_slice);
        
    end
    % smoothing 
    smoothed_slice = zeros(nRows, nCols);
    for np = 1:num
        obmask = zeros(nRows, nCols);
        obmask(labeled_slice==np) = 1;
        obmask = morphometric_closing(obmask);
        if sum(obmask(:))>0
        [~,smoothed_piece] = smooth_contours_by_cubic_spline(obmask,spline_coef);
        
        smoothed_slice = smoothed_slice + smoothed_piece;
        end
    end
    % recheck for connected lungs
    [~,num] = bwlabel(smoothed_slice);
    if num ==1
        smoothed_slice = separate_connected_lung(smoothed_slice);
        
    end
    
    smoothed_mask(:,:,n) = smoothed_slice;
end
% %% erode by 1-pixel
% se1 = strel('disk',1);
% for nsl = 1:nSlices
%     smoothed_mask(:,:,nsl) = imerode(smoothed_mask(:,:,nsl),se1);
% end
%% remove sparse voxels
[L,num] = bwlabeln(smoothed_mask);
if num>2
    for np=1:num
        if sum(L(:)==np)<10
            smoothed_mask(L==np)=0;
            
        end
    end
end

