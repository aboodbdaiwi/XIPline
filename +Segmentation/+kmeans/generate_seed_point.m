function [seedPoints, midLungmask,thorax_mask,boundingbox] = generate_seed_point(proton,bodymask)
% Ref: P. Kohlmann et al,"Automatic lung segmentation method for MRI-based
% lung perfusion studies of patients with chronic obstructive pulmonary
% disease". International Journal of Computer Assisted Radiology and
% Surgery, 2015.
%
% Step 1: Segment thorax mask: threshold at the first valley on the histogram;
% Step 2: Within the thorax mask of the mid-slice lung, threshold at the
%         first valley for the coarse lung mask;
% Step 3: Compute the Euclidean distance map to select the seeds that are
%         furthest from the lung boundaries;
%
% W. Zha@07/2016
%%
proton = double(proton);
[nRows, nCols, nSlices] = size(proton);
%% THORAX MASK
[bins, ctrs] = hist(proton(:),200);
valley_index = get_valley_index(bins);

thorax_th = ctrs(valley_index);

thorax_mask = proton;
thorax_mask(thorax_mask<thorax_th) = 0;
thorax_mask(thorax_mask > 0) = 1;
%%
% morphologic closing
se5 = strel('disk',5);
for nsl = 1: nSlices
    thorax_mask(:,:,nsl) = imfill(imclose(imfill(thorax_mask(:,:,nsl),'holes'),se5),'holes');
end
%% leakage
for nsl = 2:(nSlices-1)
    slice_thorax = thorax_mask(:,:,nsl);
    mutual = thorax_mask(:,:,nsl-1) & thorax_mask(:,:,nsl+1);
    diff = mutual - slice_thorax;
    diff(diff<0) = 0;
    [labeled_diff, nDiff] = bwlabel(diff);
    if nDiff>0
        diffPixels = zeros(nDiff,1);
        for np = 1: nDiff
            diffPixels(np) = sum(labeled_diff(:) == np);
        end
        [maxPixel,max_diff] = max(diffPixels);
        if maxPixel>1000
            
            slice_thorax(labeled_diff ==max_diff) = 1;
            thorax_mask(:,:,nsl) = imfill(slice_thorax,'holes');
        end
    end
    
end
% figure,imshow3Dfull(thorax_mask)
[labeled_thorax,num] = bwlabeln(thorax_mask,26);

if num > 1
    pixels = zeros(num,1);
    
    for np = 1:num
        pixels(np) = sum(labeled_thorax(:) == np);
    end
    [~, sorted_ind] = sort(pixels, 'descend');
    thorax_mask(labeled_thorax ~= sorted_ind(1)) = 0;
elseif num ==0
    error('thorax segmentation failed');
    
end
thorax_mask = thorax_mask.*bodymask;

%% addintional bounding box
FOV_box = zeros(nRows,nCols);
scale_factor = 8;
row_list = round(nRows/scale_factor):(nRows-round(nRows/scale_factor));
col_list = round(nCols/scale_factor):(nCols-round(nCols/scale_factor));
FOV_box(row_list,col_list)=1;
%% COARSE MID LUNG MASK
midSlice = round(nSlices/2);
midSlice_thorax = thorax_mask(:,:,midSlice);
midSlice_proton = proton(:,:,midSlice);
% [lungbins, ctrs] = hist(midSlice_proton(midSlice_thorax==1),100);
% lungbins_sm = pkasmooth(lungbins,3);
% valley_index = get_valley_index(lungbins_sm);
% 
% lung_th = ctrs(valley_index);
clustered_mid = kmeans_clustering(midSlice_proton,3);
clustered_mid(clustered_mid>1) = 0;
mid_lung = clustered_mid.*midSlice_thorax.*FOV_box;
% mid_lung = clustered_mid.*FOV_box;
% mid_lung(mid_lung >= lung_th) = 0;
% mid_lung(mid_lung < lung_th & mid_lung>0) = 1;
[labeled_mid,nLungPieces] = bwlabel(mid_lung);
% figure,data_mask_overlay(midSlice_proton,mid_lung)
%%

if nLungPieces > 0 && sum(mid_lung(:))>100
    pixels = zeros(nLungPieces,1);
    
    for np = 1:nLungPieces
        pixels(np) = sum(labeled_mid(:) == np);
    end
    
    [sorted_pixels, sorted_ind] = sort(pixels, 'descend');
    %     if length(sorted_pixels)==1 || max(sorted_pixels)<1000
    %             [maskRight,seeds.right] = regiongrowing(reg_tran_proton/max(reg_tran_proton(:)),RGthresh,seeds.right);
    %
    %         return;
    %     end
    if sorted_pixels(1) - sorted_pixels(2) < 100 || (sorted_pixels(2)>100) % two separated lungs
        %         mid_lung(labeled_mid~=sorted_ind(1) & labeled_mid ~= sorted_ind(2)) = 0;
        lung.right = zeros(nRows,nCols);
        lung.right(labeled_mid==sorted_ind(1)) = 1;
        lung.left = zeros(nRows,nCols);
        lung.left(labeled_mid==sorted_ind(2)) = 1;
    else % coneected lungs
        mid_lung(labeled_mid ~=sorted_ind(1)) = 0;
        [~, lung] = separate_connected_lung(mid_lung);
        
        
    end
else
    seedPoints.right = [];
    seedPoints.left = [];
    boundingbox = zeros(size(proton));
    return
end
lung.left(lung.right==1) = 0;
midLungmask = lung.right*2+lung.left*5;
if sum(midLungmask(:)) == 0
    seedPoints.right = [];
    seedPoints.left = [];
    boundingbox = zeros(size(proton));
    return;
    
end

%% COMPUTE EUCLUDEAN DISTANCE MAP
lung_list = fieldnames(lung);
for nl = 1: numel(lung_list)
    lungname = lung_list{nl};
    lung.(lungname) = imfill(lung.(lungname),'holes');
    lungbound.(lungname) = mask2poly(lung.(lungname),'Exact','CW');
    invalids = (lungbound.(lungname)(:,1)<0 | lungbound.(lungname)(:,2)<0 |...
        lungbound.(lungname)(:,1)>nRows | lungbound.(lungname)(:,2)>nCols);
    lungbound.(lungname)(invalids,:) = [];
    pixel_list = find(lung.(lungname)(:)==1);
    [row_list,col_list] = ind2sub([nRows,nCols],pixel_list);
    dist_list = zeros(length(pixel_list),1);
    for n = 1:length(pixel_list)
        this_row = row_list(n);
        this_col = col_list(n);
        cadindates = ((lungbound.(lungname)(:,2)-this_row).^2)+((lungbound.(lungname)(:,1)-this_col).^2);
        closest_dist = min(cadindates);
        dist_list(n) = closest_dist;
    end
    dist_map.(lungname) = zeros(nRows,nCols);
    dist_map.(lungname)(pixel_list) = dist_list;
    figure,imagesc(dist_map.(lungname));
    [~,max_dist_ind] = max(dist_list);
    hold on;plot(col_list(max_dist_ind),row_list(max_dist_ind),'x')
    seedPoints.(lungname) = [row_list(max_dist_ind), col_list(max_dist_ind),midSlice];
    
end
% lung_range = get_plotting_range(midLungmask,0);

plotting_range = get_plotting_range(midLungmask,8);
boundingbox = zeros(nRows, nCols, nSlices);
boundingbox(plotting_range.y,plotting_range.x,:) = 1;
%% ----------------------------------------
function first_min = get_valley_index(bins)
[~,global_max] = max(bins);
first_max = locate_extreme_index(smooth(bins),'max');
% if global_max < first_max || (global_max<length(bins)/4 && global_max > first_max)
%     first_max = global_max;
% end
local_min = locate_extreme_index(smooth(bins(first_max:end)),'min');
local_min = local_min + first_max - 1;
first_min = local_min;
% second_max = locate_extreme_index(bins(local_min:end),'min');
% [~,secondary_max] = max(bins(local_min:round(length(bins)/2)));
% secondary_max = secondary_max + local_min-1;
% [~, first_min] = min(bins(first_max:secondary_max));
% first_min = first_min + first_max - 1;
second_max = locate_extreme_index(smooth(bins(first_min:end)),'max');
local_min = locate_extreme_index(smooth(bins(first_min+second_max:end)),'min');
local_min = local_min + + first_min + second_max - 1;
%first_min = local_min;