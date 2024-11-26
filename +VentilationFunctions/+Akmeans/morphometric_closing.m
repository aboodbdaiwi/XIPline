function closed_lungmask = morphometric_closing(lungsmask,size_limit)
% morphometric operations:
% 1) fill the holes of size less than 500;
% 2) close the edge by a disk with radius = 2;
% 3) delete small areas of size less than 30;
%
% W. Zha @07/2016
%%
if nargin<2
    size_limit = 500;
end

closed_lungmask = lungsmask;
contoured_list = find(squeeze(sum(sum(lungsmask,1),2))>0);
se2 = strel('disk',3);
if ~isempty(contoured_list)
    for nsl = 1:length(contoured_list)
        sl = contoured_list(nsl);
        slice_lung = lungsmask(:,:,sl);
        if sum(slice_lung(:))<size_limit
            continue;
        end
        
        filled_mask = imfill(slice_lung,'holes');
        [labeled,num] = bwlabel(filled_mask-slice_lung);
        if num==0
            continue;
        end
        pixel_list = zeros(num,1);
        for n=1:num
            pixel_list(n) = sum(labeled(:)==n);
        end
        preserved_region_idx = find(pixel_list>size_limit);
        if ~isempty(preserved_region_idx)
            for n=1:length(preserved_region_idx)
                filled_mask(labeled==preserved_region_idx(n))=0;
            end
        end
        %     if sum(filled_mask(:))-sum(slice_lung(:))>size_limit
        %         filled_mask = imclose(slice_lung,se2);
        %     else
        filled_mask = imfill(imclose(filled_mask,se2),'holes');
        %     end
        closed_lungmask(:,:,sl) = bwareaopen(filled_mask,30);
    end
end