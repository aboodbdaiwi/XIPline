function lungsmask_corrected = fill_lungsmask(lungsmask,view)
%W. Zha 10/26/2016
lungsmask_corrected=lungsmask;
[nRows,nCols,nSlices] = size(lungsmask);
se2=strel('disk',2);
for nsl = 1:nSlices
    [lungpieces,num] = bwlabel(lungsmask_corrected(:,:,nsl));
    slice_lung = zeros(nRows,nCols);
    for n_p = 1:num
        piece_mask = zeros(nRows,nCols);
       piece_mask(lungpieces==n_p) = 1;
    slicefilled = imerode(imfill(imdilate(piece_mask,se2),'holes'),se2);
    
    slicediff = slicefilled - lungsmask_corrected(:,:,nsl);
    slicediff(slicediff<0)=0;
    filledComp=bwconncomp(slicediff,8);

    if strcmp(view,'axial') && max(filledComp.NumObjects)>300 %'liver dome'
%         fprintf('liver dome on slice %d. holes not filled.',nsl);
                slice_lung = slice_lung | piece_mask;

    else
        slice_lung = slice_lung | slicefilled;
    end
    end
    lungsmask_corrected(:,:,nsl) = slice_lung;
end