function [diaphragm_sl,liver_mask] = locate_diaphragm_slice(pmask,maxSL)
% estimate the slice# of diaphragm
% W. Zha@2016
if nargin<2
    slice_lungVol = squeeze(sum(sum(pmask,1),2));
    [~, maxSL] = max(slice_lungVol);
end
    
[nRows,nCols,nSlices] = size(pmask);

slice_flag = zeros(nSlices,1);
%% find diaphragm slice - the slice where liver dome appeared

liver_mask = zeros(nRows,nCols,nSlices);
for sl_no = maxSL:nSlices
    slicemask = pmask(:,:,sl_no);
    filledMask = imfill(slicemask,'holes');
    filled_region = filledMask - slicemask;
    [labeled_fill,num] = bwlabel(filled_region);
    for n = 1:num
        this_fill = zeros(nRows,nCols);
        this_fill(labeled_fill==n)=1;
    if sum(this_fill(:))>500
        
        slice_flag(sl_no) = 1;
        liver_mask(:,:,sl_no) = liver_mask(:,:,sl_no) + this_fill;
    end
    end
end
diaphragm_sl = find(slice_flag==1);
if isempty(diaphragm_sl) % If no liver dome appeared, it is the slice that slice lung volume starts to drop
    slice_Vol = squeeze(sum(sum(pmask,1),2));
    reverse_sliceVol = slice_Vol(end:-1:1);
    last_max = find_local_max(reverse_sliceVol);
    diaphragm_sl = nSlices - last_max + 1;
end
%%
