function merged_lobemask=merge_ctlobemask_into_lungmask(lobemask, binarymask)
% Use nearest neighbor to classify binarymask referring to the reference
% lobemask.
% written by Wei Zha @ 2015
%%

valid_values = unique(lobemask);
valid_values(valid_values == 0) = [];
valid_values(valid_values <1 | valid_values>=235)=[];
lobemask(lobemask<1 | lobemask>=235)=0;
[nRows,nCols,nSlices] = size(lobemask);
nSlices_tmp = size(binarymask,3);
if nSlices ~= nSlices_tmp
    error('slice_mistach:merge_ctlobemask_into_lungmask','the length of two masks mismatched..n');
    
end

merged_lobemask = lobemask.*binarymask;
unclassified_pixels = find(binarymask==1 & merged_lobemask==0);

%% reclassify pixels with undesired values
neighbor_6 = [nCols,-nCols,-1,1,nRows*nCols,-nRows*nCols];

index_list_tmp = unclassified_pixels;
nPixels = length(index_list_tmp);
nPixels_pre = 0;
while nPixels ~= nPixels_pre
    nPixels_pre = nPixels;
    for n_p = 1:nPixels
            neighbors = index_list_tmp(n_p)+neighbor_6;
            neighbors(neighbors<1 | neighbors>nRows*nCols*nSlices) = [];
            neighbor_intensities = merged_lobemask(neighbors);
            for n = 1:length(neighbors)
                if sum(valid_values == neighbor_intensities(n)) == 1
                    merged_lobemask(index_list_tmp(n_p)) = neighbor_intensities(n);
                    index_list_tmp(n_p)=0;
                    break;                                    
                end
            end
        
    end
    index_list_tmp(index_list_tmp == 0) = [];
    nPixels = length(index_list_tmp);
end
% toc
if nPixels ~= 0
    warning('pixel_reclassification:clean_up_lobemask','%d isolated pixels detected.\n',nPixels);
end